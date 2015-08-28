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

#ifndef VELMA_Q5Q6_COLLISION_H__
#define VELMA_Q5Q6_COLLISION_H__

#include "Eigen/Dense"

#include "planer_utils/marker_publisher.h"

class VelmaQ5Q6CollisionChecker {
public:
    VelmaQ5Q6CollisionChecker(int q5_idx, int q6_idx, double d0, bool inverted);
    ~VelmaQ5Q6CollisionChecker();
    bool inCollision(const Eigen::VectorXd &q) const;
    bool getMinDistance(const Eigen::VectorXd &q, Eigen::VectorXd &min_v, double &min_dist, int &min_idx, int &min_type) const;
    int visualizeBorder(MarkerPublisher *markers_pub, int m_id) const;
    int visualizeRegion(MarkerPublisher *markers_pub, int m_id, int min_idx, int min_type) const;
protected:
    bool point_inside_polygon(const Eigen::VectorXd &p) const;

    int q5_idx_, q6_idx_;
    double d0_;

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

};

#endif  // VELMA_Q5Q6_COLLISION_H__

