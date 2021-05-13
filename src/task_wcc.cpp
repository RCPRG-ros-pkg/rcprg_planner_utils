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

    Task_WCC::Task_WCC(int ndof, int q5_idx, int q6_idx, const std::vector<double >& polygon) :
        ndof_(ndof),
        q5_idx_(q5_idx),
        q6_idx_(q6_idx),
        d0_(10.0/180.0*3.1415),
        in_collision_(false),
        af_(0.2 * d0_, 4.0 / d0_),
        cc_(d0_, polygon),
        activation_(0.0)
    {
    }

    Task_WCC::~Task_WCC() {
    }

    void Task_WCC::compute(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::MatrixXd &invI, Eigen::VectorXd &torque, Eigen::MatrixXd &N, MarkerPublisher *markers_pub, int *m_id) {

            torque.setZero();
            N.setIdentity();

            DoubleJointCC::Joints q2(q(q5_idx_), q(q6_idx_));

            if (markers_pub != NULL) {
                *m_id = cc_.visualizeBorder(markers_pub, *m_id);
                *m_id = markers_pub->addSinglePointMarker(*m_id, KDL::Vector(q(q5_idx_), q(q6_idx_), 0), 0, 1, 0, 1, 0.03, "world");
            }

            if (cc_.inCollision(q2)) {
                in_collision_ = true;
                return;
            }
            else {
                in_collision_ = false;
            }

            int min_idx;
            int min_type;
            double min_dist;
            DoubleJointCC::Joints min_v;

            bool found = cc_.getMinDistanceIn(q2, min_v, min_dist, min_idx, min_type);

            if (markers_pub != NULL) {
                if (found) {
                    *m_id = cc_.visualizeRegion(markers_pub, *m_id, min_idx, min_type);
                    *m_id = markers_pub->addVectorMarker(*m_id, KDL::Vector(q(q5_idx_), q(q6_idx_), 0), KDL::Vector(q(q5_idx_)+min_v(0), q(q6_idx_)+min_v(1), 0), 0, 1, 0, 1, 0.02, "world");
                }
                markers_pub->addEraseMarkers(*m_id, *m_id+100);
            }

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

            double Fmax_ = 50.0;
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

            activation_ = 1.0 - af_.func_Ndes(min_dist);

            N = Eigen::MatrixXd::Identity(ndof_, ndof_) - (J.transpose() * activation_ * J);

            // calculate collision mass (1 dof)
            double Mdij_inv = (J * invI * J.transpose())(0,0);

            double D = 2.0 * 0.7 * sqrt(Mdij_inv * K);  // sqrt(K/M)
            Eigen::VectorXd d_torque = J.transpose() * (Frep - D * ddij);
            torque += d_torque;
    }

int Task_WCC::getActivationCount() const {
    if (activation_ > 0.001) {
        return 1;
    }
    return 0;
}

    bool Task_WCC::inCollision() const {
        return in_collision_;
    }

