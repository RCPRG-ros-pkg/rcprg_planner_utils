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

#include "planer_utils/task_wcc.h"
#include <iostream>
#include <math.h>

const double Task_WCC::polygon_q5q6_[] = {
-0.397855401039,-2.90307354927,
2.12894010544,-2.90307354927,
1.92475450039,-1.43123674393,
0.77621114254,-1.39720571041,
0.350824713707,-1.00585031509,
0.401871085167,-0.571956157684,
0.810242056847,0.414940297604,
1.34622907639,0.942419290543,
2.11192464828,1.01898884773,
2.12894010544,2.8906891346,
-0.814733862877,2.8906891346,
-1.22310483456,2.27813267708,
-2.21850919724,2.29514837265,
-2.22701668739,-1.32063627243,
-0.814733862877,-1.73751521111,
-0.423378348351,-2.09483933449,
};

    Task_WCC::Task_WCC(int ndof, int q5_idx, int q6_idx) :
        ndof_(ndof),
        q5_idx_(q5_idx),
        q6_idx_(q6_idx),
        d0_(20.0/180.0*3.1415)
    {
        int lines_count = sizeof(polygon_q5q6_)/sizeof(double)/2;

        for (int line_idx = 0; line_idx < lines_count; line_idx++) {
            int idxA = line_idx * 2;
            int idxB = ((line_idx + 1) % lines_count) * 2;
            Line l;
            l.a.resize(2);
            l.a(0) = polygon_q5q6_[idxA];
            l.a(1) = polygon_q5q6_[idxA + 1];
            l.b.resize(2);
            l.b(0) = polygon_q5q6_[idxB];
            l.b(1) = polygon_q5q6_[idxB + 1];
            l.n.resize(2);
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
            std::cout << l1_idx << " " << l2_idx << std::endl;
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
                Eigen::VectorXd n(lines_[l1_idx].n + lines_[l2_idx].n);
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

    Task_WCC::~Task_WCC() {
    }

    void Task_WCC::distanceLinePoint(const KDL::Vector &lineA, const KDL::Vector &lineB, const KDL::Vector &pt, double &distance, KDL::Vector &p_pt1, KDL::Vector &p_pt2) {
        KDL::Vector v = lineB - lineA;
        double ta = KDL::dot(v, lineA);
        double tb = KDL::dot(v, lineB);
        double tpt = KDL::dot(v, pt);
        if (tpt <= ta) {
            distance = (lineA-pt).Norm();
            p_pt1 = lineA;
            p_pt2 = pt;
        }
        else if (tpt >= tb) {
            distance = (lineB-pt).Norm();
            p_pt1 = lineB;
            p_pt2 = pt;
        }
        else {
            KDL::Vector n(v.y(), -v.x(), v.z());
            n.Normalize();
            double diff = KDL::dot(n, lineA) - KDL::dot(n, pt);
            distance = fabs(diff);
            p_pt1 = pt + (diff * n);
            p_pt2 = pt;
        }
    }

    void Task_WCC::compute(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::MatrixXd &I, const Eigen::MatrixXd &invI, Eigen::VectorXd &torque, Eigen::MatrixXd &N, MarkerPublisher &markers_pub) {

            Eigen::VectorXd p(2);
            p(0) = q(q5_idx_);
            p(1) = q(q6_idx_);
            KDL::Vector pt(q(q5_idx_), q(q6_idx_), 0);
            int m_id = 2000;

            int min_idx = -1;
            int min_type = -1;
            double min_dist = d0_;
            Eigen::VectorXd min_v(2);
            bool found = false;

            for (int line_idx = 0; line_idx < lines_.size(); line_idx++) {
                KDL::Vector lineA(lines_[line_idx].a(0), lines_[line_idx].a(1), 0);
                KDL::Vector lineB(lines_[line_idx].b(0), lines_[line_idx].b(1), 0);
                KDL::Vector n(lines_[line_idx].n(0), lines_[line_idx].n(1), 0);
                m_id = markers_pub.addVectorMarker(m_id, lineA, lineB, 1, 1, 1, 1, 0.02, "world");
                m_id = markers_pub.addVectorMarker(m_id, (lineA + lineB)/2, (lineA + lineB)/2 + n*0.2, 0, 0, 1, 1, 0.02, "world");

                if (lines_[line_idx].ab.dot(p) >= lines_[line_idx].da && lines_[line_idx].ab.dot(p) <= lines_[line_idx].db) {
                    double dist = lines_[line_idx].n.dot(p) - lines_[line_idx].dn;
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
                if (convex_points_[pt_idx].n1.dot(p) >= convex_points_[pt_idx].dn1 && convex_points_[pt_idx].n2.dot(p) >= convex_points_[pt_idx].dn2) {
                    Eigen::VectorXd v(p - convex_points_[pt_idx].p);
                    double dist = v.norm();
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_v = v;
                        found = true;
                        min_idx = pt_idx;
                        min_type = 1;
                    }
                }
                KDL::Vector pt(convex_points_[pt_idx].p(0), convex_points_[pt_idx].p(1), 0);
                m_id = markers_pub.addSinglePointMarker(m_id, pt, 0, 1, 0, 1, 0.1, "world");
            }

            for (int pt_idx = 0; pt_idx < concave_points_.size(); pt_idx++) {
                if (concave_points_[pt_idx].n1.dot(p) >= concave_points_[pt_idx].dn1 && concave_points_[pt_idx].n2.dot(p) >= concave_points_[pt_idx].dn2) {
                    Eigen::VectorXd v(p - concave_points_[pt_idx].p);
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
                KDL::Vector pt(concave_points_[pt_idx].p(0), concave_points_[pt_idx].p(1), 0);
                m_id = markers_pub.addSinglePointMarker(m_id, pt, 1, 0, 0, 1, 0.1, "world");
            }

            if (found) {
                std::cout << p(0) << " " << p(1) << " min_dist: " << min_dist << std::endl;
                if (min_type == 0) {
                    KDL::Vector lineA(lines_[min_idx].a(0), lines_[min_idx].a(1), 0);
                    KDL::Vector lineB(lines_[min_idx].b(0), lines_[min_idx].b(1), 0);
                    m_id = markers_pub.addVectorMarker(m_id, lineA, lineB, 1, 1, 1, 1, 0.04, "world");
                }
                else if (min_type == 1) {
                    KDL::Vector pt(convex_points_[min_idx].p(0), convex_points_[min_idx].p(1), 0);
                    m_id = markers_pub.addSinglePointMarker(m_id, pt, 0, 1, 0, 1, 0.2, "world");
                }
                else if (min_type == 2) {
                    KDL::Vector pt(concave_points_[min_idx].p(0), concave_points_[min_idx].p(1), 0);
                    m_id = markers_pub.addSinglePointMarker(m_id, pt, 0, 1, 0, 1, 0.2, "world");
                }
                m_id = markers_pub.addVectorMarker(m_id, KDL::Vector(p(0), p(1), 0), KDL::Vector(p(0)+min_v(0), p(1)+min_v(1), 0), 0, 1, 0, 1, 0.02, "world");
            }
            else {
                m_id = markers_pub.addSinglePointMarker(m_id, KDL::Vector(p(0), p(1), 0), 0, 1, 0, 1, 0.03, "world");
            }

            markers_pub.addEraseMarkers(m_id, m_id+100);

            for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                torque(q_idx) = 0.0;
            }
            N = Eigen::MatrixXd::Identity(ndof_, ndof_);

            if (!found) {
                return;
            }

            min_v.normalize();

            double depth = d0_ - min_dist;
            if (depth > d0_) {
                depth = d0_;
            }
            else if (depth < 0.0) {
                depth = 0.0;
            }

            double Fmax_ = 1.0;
            double f = depth / d0_;
            double Frep = Fmax_ * f * f;
            double K = 2.0 * Fmax_ / (d0_ * d0_);

            Eigen::MatrixXd J(1, ndof_);
            for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                J(0, q_idx) = 0.0;
            }

            J(0, q5_idx_) = min_v(0);
            J(0, q6_idx_) = min_v(1);

            // calculate relative velocity between points (1 dof)
            double ddij = (J * dq)(0,0);

            double activation = 5.0*f;
            if (activation > 1.0) {
                activation = 1.0;
            }
            if (activation < 0.0) {
                activation = 0.0;
            }
//            if (ddij <= 0.0) {
//                activation = 0.0;
//            }

            N = Eigen::MatrixXd::Identity(ndof_, ndof_) - (J.transpose() * activation * J);

            // calculate collision mass (1 dof)
            double Mdij_inv = (J * invI * J.transpose())(0,0);

            double D = 2.0 * 0.7 * sqrt(Mdij_inv * K);  // sqrt(K/M)
            Eigen::VectorXd d_torque = J.transpose() * (Frep - D * ddij);
            torque += d_torque;
    }

