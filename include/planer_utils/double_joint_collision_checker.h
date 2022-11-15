// Copyright (c) 2015-2017, Robot Control and Pattern Recognition Group,
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

#ifndef DOUBLE_JOINT_COLLISION_CHECKER_H__
#define DOUBLE_JOINT_COLLISION_CHECKER_H__

#include "Eigen/Dense"

#include "rcprg_ros_utils/marker_publisher.h"

class DoubleJointCC {
public:

    typedef Eigen::Matrix<double, 2, 1 > Joints;

    DoubleJointCC(double d0, const std::vector<double >& polygon);
    ~DoubleJointCC();
    bool inCollision(const Joints &q) const;
    bool getMinDistanceIn(const Joints &q, Joints &min_v, double &min_dist, int &min_idx, int &min_type) const;
    bool getMinDistanceOut(const Joints &q, Joints &min_v, double &min_dist, int &min_idx, int &min_type) const;
    double getD0() const;
    int visualizeBorder(MarkerPublisher *markers_pub, int m_id) const;
    int visualizeRegion(MarkerPublisher *markers_pub, int m_id, int min_idx, int min_type) const;
protected:
    bool point_inside_polygon(const Joints &p) const;

    double d0_;

    class Line {
    public:
        Joints a, b, n, ab;
        double dn, da, db;
    };

    class PointAngle {
    public:
        Joints p, n1, n2;
        double dn1, dn2;
    };

    std::vector<PointAngle > convex_points_;
    std::vector<PointAngle > concave_points_;
    std::vector<Line > lines_;

};

#endif  // DOUBLE_JOINT_COLLISION_CHECKER_H__

