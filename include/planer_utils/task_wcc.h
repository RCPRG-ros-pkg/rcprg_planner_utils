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
#include "planer_utils/activation_function.h"
#include "planer_utils/double_joint_collision_checker.h"

class Task_WCC {
public:
    Task_WCC(int ndof, int q5_idx, int q6_idx, const std::vector<double >& polygon);

    ~Task_WCC();

    void compute(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::MatrixXd &invI, Eigen::VectorXd &torque, Eigen::MatrixXd &N, MarkerPublisher *markers_pub=NULL, int *m_id=NULL);

    bool inCollision() const;

    bool getMinDistance(const Eigen::VectorXd &p, double d0, Eigen::VectorXd &min_v, double &min_dist, int &min_idx, int &min_type) const;

    int getActivationCount() const;

    int visualizeBorder(MarkerPublisher *markers_pub, int m_id) const;
    int visualizeRegion(MarkerPublisher *markers_pub, int m_id, int min_idx, int min_type) const;

protected:

    int ndof_;
    int q5_idx_, q6_idx_;

    double d0_;
    double activation_;
    bool in_collision_;

    DoubleJointCC cc_;

    ActivationFunction af_;
};

#endif  // TASK_WCC_H__

