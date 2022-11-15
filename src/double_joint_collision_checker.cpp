// Copyright (c) 2015, Robot Control and Pattern Recognition Group,
// Institute of Control and Computation Engineering
// Warsaw University of Technology
//
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Warsaw University of Technology nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYright HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author: Dawid Seredynski
//

#include "planer_utils/double_joint_collision_checker.h"
#include <iostream>
#include <math.h>

    DoubleJointCC::DoubleJointCC(double d0, const std::vector<double >& polygon) :
        d0_(d0)
    {
        int lines_count = polygon.size() / 2;

        for (int line_idx = 0; line_idx < lines_count; line_idx++) {
            int idxA = line_idx * 2;
            int idxB = ((line_idx + 1) % lines_count) * 2;
            Line l;
            l.a(0) = polygon[idxA];
            l.a(1) = polygon[idxA + 1];
            l.b(0) = polygon[idxB];
            l.b(1) = polygon[idxB + 1];
            l.n(0) = l.a(1) - l.b(1);
            l.n(1) = l.b(0) - l.a(0);
            l.n.normalize();
            l.dn = l.a.dot(l.n);
            l.ab = l.b - l.a;
            lines_.push_back(l);
        }

        for (int pt_idx = 0; pt_idx < lines_count; pt_idx++) {
            int l1_idx = (pt_idx - 1 + lines_count) % lines_count;
            int l2_idx = pt_idx;
            if (lines_[l1_idx].n.dot(lines_[l2_idx].b) < lines_[l1_idx].dn) {
                // convex
                PointAngle pa;
                pa.p = lines_[l2_idx].a;
                pa.n1 = lines_[l1_idx].b - lines_[l1_idx].a;
                pa.n1.normalize();
                pa.dn1 = pa.n1.dot(pa.p);
                pa.n2 = lines_[l2_idx].a - lines_[l2_idx].b;
                pa.n2.normalize();
                pa.dn2 = pa.n2.dot(pa.p);
                lines_[l1_idx].db = lines_[l1_idx].ab.dot(pa.p);
                lines_[l2_idx].da = lines_[l2_idx].ab.dot(pa.p);
                convex_points_.push_back(pa);
            }
            else {
                // concave
                DoubleJointCC::Joints n(lines_[l1_idx].n + lines_[l2_idx].n);
                n.normalize();
                PointAngle pa;

                Eigen::Vector3d n1(lines_[l1_idx].n(0), lines_[l1_idx].n(1), 0), n2(lines_[l2_idx].n(0), lines_[l2_idx].n(1), 0);
                double angle = atan2((n1.cross(n2)).norm(), lines_[l1_idx].n.dot(lines_[l2_idx].n));

                pa.p = lines_[l2_idx].a + n * d0_ * 2.0 / cos(angle/2.0);
                pa.n1 = lines_[l1_idx].b - lines_[l1_idx].a;
                pa.n1.normalize();
                pa.dn1 = pa.n1.dot(pa.p);
                pa.n2 = lines_[l2_idx].a - lines_[l2_idx].b;
                pa.n2.normalize();
                pa.dn2 = pa.n2.dot(pa.p);
                lines_[l1_idx].db = lines_[l1_idx].ab.dot(pa.p);
                lines_[l2_idx].da = lines_[l2_idx].ab.dot(pa.p);
                concave_points_.push_back(pa);
            }
        }
    }

    DoubleJointCC::~DoubleJointCC() {
    }

    bool DoubleJointCC::inCollision(const DoubleJointCC::Joints &q) const {
        return !point_inside_polygon(q);
    }

    bool DoubleJointCC::point_inside_polygon(const DoubleJointCC::Joints &p) const {
        bool inside = false;
        for (int line_idx = 0; line_idx < lines_.size(); line_idx++) {
            const DoubleJointCC::Joints &a = lines_[line_idx].a;
            const DoubleJointCC::Joints &b = lines_[line_idx].b;
            if (p(1) > a(1) || p(1) > b(1)) {
                if (p(1) <= a(1) || p(1) <= b(1)) {
                    if (p(0) <= a(0) || p(0) <= b(0)) {
                        double xinters = 0.0;
                        if (a(1) != b(1)) {
                            xinters = (p(1)-a(1))*(b(0)-a(0))/(b(1)-a(1))+a(0);
                        }
                        if (a(0) == b(0) || p(0) <= xinters) {
                            inside = !inside;
                        }
                    }
                }
            }
        }

        return inside;
    }

    bool DoubleJointCC::getMinDistanceIn(const DoubleJointCC::Joints &q, DoubleJointCC::Joints &min_v, double &min_dist, int &min_idx, int &min_type) const {
            bool found = false;
            min_dist = d0_;

            for (int line_idx = 0; line_idx < lines_.size(); line_idx++) {
                if (lines_[line_idx].ab.dot(q) >= lines_[line_idx].da && lines_[line_idx].ab.dot(q) <= lines_[line_idx].db) {
                    double dist = lines_[line_idx].n.dot(q) - lines_[line_idx].dn;
                    if (dist < min_dist && dist > 0) {
                        min_dist = dist;
                        min_v = lines_[line_idx].n * dist;
                        found = true;
                        min_idx = line_idx;
                        min_type = 0;
                    }
                }
            }

            for (int pt_idx = 0; pt_idx < convex_points_.size(); pt_idx++) {
                if (convex_points_[pt_idx].n1.dot(q) >= convex_points_[pt_idx].dn1 && convex_points_[pt_idx].n2.dot(q) >= convex_points_[pt_idx].dn2) {
                    DoubleJointCC::Joints v(q - convex_points_[pt_idx].p);
                    double dist = v.norm();
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_v = v;
                        found = true;
                        min_idx = pt_idx;
                        min_type = 1;
                    }
                }
            }

            for (int pt_idx = 0; pt_idx < concave_points_.size(); pt_idx++) {
                if (concave_points_[pt_idx].n1.dot(q) >= concave_points_[pt_idx].dn1 && concave_points_[pt_idx].n2.dot(q) >= concave_points_[pt_idx].dn2) {
                    DoubleJointCC::Joints v(q - concave_points_[pt_idx].p);
                    double dist = 2.0 * d0_ - v.norm();
                    if (dist < min_dist && dist > 0) {
                        min_dist = dist;
                        v.normalize();
                        min_v = - v * dist;
                        found = true;
                        min_idx = pt_idx;
                        min_type = 2;
                    }
                }
            }
            return found;
    }

    bool DoubleJointCC::getMinDistanceOut(const DoubleJointCC::Joints &q, DoubleJointCC::Joints &min_v, double &min_dist, int &min_idx, int &min_type) const {
            bool found = false;
            min_dist = d0_;

            for (int line_idx = 0; line_idx < lines_.size(); line_idx++) {
                if (lines_[line_idx].ab.dot(q) >= lines_[line_idx].da && lines_[line_idx].ab.dot(q) <= lines_[line_idx].db) {
                    double dist = lines_[line_idx].n.dot(q) - lines_[line_idx].dn;
                    dist *= -1.0;
                    if (dist < min_dist && dist > 0) {
                        min_dist = dist;
                        min_v = lines_[line_idx].n * dist;
                        found = true;
                        min_idx = line_idx;
                        min_type = 0;
                    }
                }
            }

            for (int pt_idx = 0; pt_idx < convex_points_.size(); pt_idx++) {
                if (convex_points_[pt_idx].n1.dot(q) <= convex_points_[pt_idx].dn1 && convex_points_[pt_idx].n2.dot(q) <= convex_points_[pt_idx].dn2) {
                    DoubleJointCC::Joints v(q - convex_points_[pt_idx].p);
                    double dist = v.norm();
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_v = -v;
                        found = true;
                        min_idx = pt_idx;
                        min_type = 1;
                    }
                }
            }

            for (int pt_idx = 0; pt_idx < concave_points_.size(); pt_idx++) {
                if (concave_points_[pt_idx].n1.dot(q) >= concave_points_[pt_idx].dn1 && concave_points_[pt_idx].n2.dot(q) >= concave_points_[pt_idx].dn2) {
                    DoubleJointCC::Joints v(q - concave_points_[pt_idx].p);
                    double dist = v.norm() - 2.0 * d0_;
                    if (dist < min_dist && dist > 0)
                    {
                        min_dist = dist;
                        v.normalize();
                        min_v = -v * dist;
                        found = true;
                        min_idx = pt_idx;
                        min_type = 2;
                    }
                }
            }
            return found;
    }

    double DoubleJointCC::getD0() const {
        return d0_;
    }

    int DoubleJointCC::visualizeBorder(MarkerPublisher *markers_pub, int m_id) const {
                // draw border lines
                for (int line_idx = 0; line_idx < lines_.size(); line_idx++) {
                    KDL::Vector n(lines_[line_idx].n(0), lines_[line_idx].n(1), 0);
                    KDL::Vector lineA(lines_[line_idx].a(0), lines_[line_idx].a(1), 0);
                    KDL::Vector lineB(lines_[line_idx].b(0), lines_[line_idx].b(1), 0);
                    m_id = markers_pub->addVectorMarker(m_id, lineA, lineB, 1, 1, 1, 1, 0.02, "world");
                    m_id = markers_pub->addVectorMarker(m_id, (lineA + lineB)/2, (lineA + lineB)/2 + n*0.2, 0, 0, 1, 1, 0.02, "world");
                }
                // draw convex points
                for (int pt_idx = 0; pt_idx < convex_points_.size(); pt_idx++) {
                    KDL::Vector pt(convex_points_[pt_idx].p(0), convex_points_[pt_idx].p(1), 0);
                    m_id = markers_pub->addSinglePointMarker(m_id, pt, 0, 1, 0, 1, 0.1, "world");
                }
                // draw concave points
                for (int pt_idx = 0; pt_idx < concave_points_.size(); pt_idx++) {
                    KDL::Vector pt(concave_points_[pt_idx].p(0), concave_points_[pt_idx].p(1), 0);
                    m_id = markers_pub->addSinglePointMarker(m_id, pt, 1, 0, 0, 1, 0.1, "world");
                }
                return m_id;
    }

    int DoubleJointCC::visualizeRegion(MarkerPublisher *markers_pub, int m_id, int min_idx, int min_type) const {
                // draw enlarged region (closest to the point)
                if (min_type == 0) {
                    KDL::Vector lineA(lines_[min_idx].a(0), lines_[min_idx].a(1), 0);
                    KDL::Vector lineB(lines_[min_idx].b(0), lines_[min_idx].b(1), 0);
                    m_id = markers_pub->addVectorMarker(m_id, lineA, lineB, 1, 1, 1, 1, 0.04, "world");
                }
                else if (min_type == 1) {
                    KDL::Vector pt(convex_points_[min_idx].p(0), convex_points_[min_idx].p(1), 0);
                    m_id = markers_pub->addSinglePointMarker(m_id, pt, 0, 1, 0, 1, 0.2, "world");
                }
                else if (min_type == 2) {
                    KDL::Vector pt(concave_points_[min_idx].p(0), concave_points_[min_idx].p(1), 0);
                    m_id = markers_pub->addSinglePointMarker(m_id, pt, 1, 0, 0, 1, 0.2, "world");
                }
                return m_id;
    }

