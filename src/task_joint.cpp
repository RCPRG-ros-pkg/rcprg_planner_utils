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

#include "planer_utils/task_joint.h"

Task_JOINT::Task_JOINT(int ndof, const std::set<int > &excluded_q_idx) :
    ndof_(ndof),
    activation_(ndof_),
//    k_(ndof_),
    k0_(ndof_),
    q_(ndof_, ndof_), d_(ndof_, ndof_),
    J_(ndof_, ndof_),
    tmpNN_(ndof_, ndof_),
    excluded_q_idx_(excluded_q_idx)
{
    for (int q_idx = 0; q_idx < ndof_; q_idx++) {
        if (excluded_q_idx_.find(q_idx) != excluded_q_idx_.end()) {
            activation_[q_idx] = 0.0;
        }
        else {
            activation_[q_idx] = 1.0;
        }
    }
}

Task_JOINT::~Task_JOINT() {
}

void Task_JOINT::compute(const Eigen::VectorXd &q_diff, const Eigen::VectorXd &Kc, const Eigen::VectorXd &Dxi, const Eigen::VectorXd &dq, const Eigen::MatrixXd &I, Eigen::VectorXd &torque, Eigen::MatrixXd &N) {
            // code from joint_limit_avoidance.cpp
            for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                if (excluded_q_idx_.find(q_idx) != excluded_q_idx_.end()) {
                    torque(q_idx) = 0.0;
                }
                else {
                    torque[q_idx] = Kc[q_idx] * q_diff[q_idx];
                }
                if (torque[q_idx] != torque[q_idx]) {
                    std::cout << "Task_JOINT::compute: torque: " << torque.transpose() << std::endl;
                    std::cout << "Task_JOINT::compute: q_diff: " << q_diff.transpose() << std::endl;
                    std::cout << "Task_JOINT::compute: Kc: " << Kc.transpose() << std::endl;
                    return;
                }
            }

            tmpNN_ = Kc.asDiagonal();
            es_.compute(tmpNN_, I);
            q_ = es_.eigenvectors().inverse();
            k0_ = es_.eigenvalues();

            tmpNN_ = k0_.cwiseSqrt().asDiagonal();

            d_.noalias() = 2.0 * q_.adjoint() * Dxi * tmpNN_ * q_;

            torque.noalias() -= d_ * dq;

            for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                if (torque[q_idx] != torque[q_idx]) {
                    std::cout << "Task_JOINT::compute2: torque: " << torque.transpose() << std::endl;
                    std::cout << "Task_JOINT::compute2: dq: " << dq.transpose() << std::endl;
                    std::cout << "Task_JOINT::compute2: d_: " << d_ << std::endl;
                    std::cout << "Task_JOINT::compute2: q_: " << q_ << std::endl;
                    std::cout << "Task_JOINT::compute2: Dxi: " << Dxi << std::endl;
                    std::cout << "Task_JOINT::compute2: tmpNN_: " << tmpNN_ << std::endl;
                    return;
                }
            }
            // calculate jacobian (the activation function)
            J_ = activation_.asDiagonal();
            N = Eigen::MatrixXd::Identity(ndof_, ndof_) - (J_.transpose() * J_);
}

int Task_JOINT::visualize(MarkerPublisher *markers_pub, int m_id, const boost::shared_ptr<KinematicModel> &kin_model, const boost::shared_ptr<self_collision::CollisionModel> &col_model, const std::vector<KDL::Frame > &links_fk) const {
    for (int q_idx = 0; q_idx < ndof_; q_idx++) {
        if (excluded_q_idx_.find(q_idx) == excluded_q_idx_.end()) {
            KDL::Vector axis, origin;
            std::string link_name;
            kin_model->getJointAxisAndOrigin(q_idx, axis, origin);
            origin = KDL::Vector();
            kin_model->getJointLinkName(q_idx, link_name);
            const KDL::Frame &T_B_L = links_fk[col_model->getLinkIndex(link_name)];

            if (activation_[q_idx] < 0.001) {
                m_id = markers_pub->addVectorMarker(m_id, T_B_L * origin, T_B_L * (origin + 0.4 * axis), 0, 1, 0, 1, 0.01, "world");
            }
            else {
                m_id = markers_pub->addVectorMarker(m_id, T_B_L * origin, T_B_L * (origin + 0.4 * axis), 1, activation_[q_idx], activation_[q_idx], 1, 0.01, "world");
            }
        }
    }
    return m_id;
}

