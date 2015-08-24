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

#ifndef TASK_WCC_H__
#define TASK_WCC_H__

#include "Eigen/Dense"
#include "Eigen/LU"

#include "planer_utils/marker_publisher.h"

class Task_WCC {
public:
    Task_WCC(int ndof, int q5_idx, int q6_idx);

    ~Task_WCC();

    void distanceLinePoint(const KDL::Vector &lineA, const KDL::Vector &lineB, const KDL::Vector &pt, double &distance, KDL::Vector &p_pt1, KDL::Vector &p_pt2);

    void compute(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::MatrixXd &I, const Eigen::MatrixXd &invI, Eigen::VectorXd &torque, Eigen::MatrixXd &N, MarkerPublisher &markers_pub);

protected:
    int ndof_;
    int q5_idx_, q6_idx_;

    static const double polygon_q5q6_[];

    class Line {
    public:
        Eigen::VectorXd a, b, n, ab;
        double dn, da, db;
    };

    class PointAngle {
    public:
        Eigen::VectorXd p, n1, n2;
        double dn1, dn2;
    };

    std::vector<PointAngle > convex_points_;
    std::vector<PointAngle > concave_points_;
    std::vector<Line > lines_;

    double d0_;
};

#endif  // TASK_WCC_H__

