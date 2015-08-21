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

#include "planer_utils/simulator.h"

    DynamicsSimulatorHandPose::DynamicsSimulatorHandPose(int ndof, int ndim, const std::string &effector_name, const boost::shared_ptr<self_collision::CollisionModel> &col_model,
                        const boost::shared_ptr<KinematicModel> &kin_model,
                        const boost::shared_ptr<DynamicModel > &dyn_model,
                        const std::vector<std::string > &joint_names) :
        ndof_(ndof),
        ndim_(ndim),
        effector_name_(effector_name),
        q_(ndof_),
        dq_(ndof_),
        ddq_(ndof_),
        torque_(ndof_),
        col_model_(col_model),
        kin_model_(kin_model),
        dyn_model_(dyn_model),
        links_fk_(col_model->getLinksCount()),
        activation_dist_(0.05)
    {
        for (int q_idx = 0; q_idx < ndof; q_idx++) {
            q_[q_idx] = 0.0;
            dq_[q_idx] = 0.0;
            ddq_[q_idx] = 0.0;
            torque_[q_idx] = 0.0;
        }

        effector_idx_ = col_model->getLinkIndex(effector_name_);

        // joint limits
        Eigen::VectorXd lower_limit(ndof), upper_limit(ndof), limit_range(ndof), max_trq(ndof);
        int q_idx = 0;
        for (std::vector<std::string >::const_iterator name_it = joint_names.begin(); name_it != joint_names.end(); name_it++, q_idx++) {

            if (!col_model->getJointLimits( (*name_it), lower_limit[q_idx], upper_limit[q_idx] )) {
                std::cout << "ERROR: could not find joint with name " << (*name_it) << std::endl;
                return;
            }
            limit_range[q_idx] = 10.0 / 180.0 * 3.141592653589793;
            max_trq[q_idx] = 10.0;
        }

        if (ndim_ != 3 && ndim_ != 6) {
            std::cout << "ERROR: DynamicsSimulatorHandPose works for 3 or 6 dimensions" << std::endl;
        }

        task_JLC_.reset( new Task_JLC(lower_limit, upper_limit, limit_range, max_trq) );
        task_COL_.reset( new Task_COL(ndof_, activation_dist_, 10.0, kin_model_, col_model_) );
        task_HAND_.reset( new Task_HAND(ndof_, ndim_) );

        J_r_HAND_6_.resize(6, ndof_);
        J_r_HAND_.resize(ndim_, ndof_);
    }

    void DynamicsSimulatorHandPose::setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq) {
        q_ = q;
        dq_ = dq;
        ddq_ = ddq;
    }

    void DynamicsSimulatorHandPose::getState(Eigen::VectorXd &q, Eigen::VectorXd &dq, Eigen::VectorXd &ddq) {
        q = q_;
        dq = dq_;
        ddq = ddq_;
    }

    void DynamicsSimulatorHandPose::setTarget(const KDL::Frame &r_HAND_target) {
        r_HAND_target_ = r_HAND_target;
    }

    void DynamicsSimulatorHandPose::oneStep(const KDL::Twist &diff) {
                // calculate forward kinematics for all links
                for (int l_idx = 0; l_idx < col_model_->getLinksCount(); l_idx++) {
                    kin_model_->calculateFk(links_fk_[l_idx], col_model_->getLinkName(l_idx), q_);
                }
                // calculate inertia matrix for whole body
                dyn_model_->computeM(q_);

                //
                // joint limit avoidance
                //
                Eigen::VectorXd torque_JLC(ndof_);
                Eigen::MatrixXd N_JLC(ndof_, ndof_);
                task_JLC_->compute(q_, dq_, dyn_model_->getM(), torque_JLC, N_JLC);

                //
                // collision constraints
                //
                std::vector<self_collision::CollisionInfo> link_collisions;
                self_collision::getCollisionPairs(col_model_, links_fk_, activation_dist_, link_collisions);

                Eigen::VectorXd torque_COL(ndof_);
                for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                    torque_COL[q_idx] = 0.0;
                }
                Eigen::MatrixXd N_COL(Eigen::MatrixXd::Identity(ndof_, ndof_));

                task_COL_->compute(q_, dq_, dyn_model_->getInvM(), links_fk_, link_collisions, torque_COL, N_COL);

                //
                // effector task
                //

                Eigen::VectorXd torque_HAND(ndof_);
                Eigen::MatrixXd N_HAND(ndof_, ndof_);

                Eigen::VectorXd r_HAND_diff(ndim_);
                if (ndim_ == 3) {
                    r_HAND_diff[0] = diff[0];
                    r_HAND_diff[1] = diff[1];
                    r_HAND_diff[2] = diff[5];
                }
                else {
                    for (int dim_idx = 0; dim_idx < 6; dim_idx++) {
                        r_HAND_diff[dim_idx] = diff[dim_idx];
                    }
                }

                double Kc_lin = 20.0;
                double Kc_rot = 2.0;
                Eigen::VectorXd Kc(ndim_);
                if (ndim_ == 3) {
                    Kc[0] = Kc_lin;
                    Kc[1] = Kc_lin;
                    Kc[2] = Kc_rot;
                }
                else {
                    for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
                        Kc[dim_idx] = Kc_lin;
                        Kc[dim_idx+3] = Kc_rot;
                    }
                }

                Eigen::VectorXd Dxi(ndim_);
                for (int dim_idx = 0; dim_idx < ndim_; dim_idx++) {
                    Dxi[dim_idx] = 0.7;
                }

                kin_model_->getJacobian(J_r_HAND_6_, effector_name_, q_);

                if (ndim_ == 3) {
                    for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                        J_r_HAND_(0, q_idx) = J_r_HAND_6_(0, q_idx);
                        J_r_HAND_(1, q_idx) = J_r_HAND_6_(1, q_idx);
                        J_r_HAND_(2, q_idx) = J_r_HAND_6_(5, q_idx);
                    }
                    task_HAND_->compute(r_HAND_diff, Kc, Dxi, J_r_HAND_, dq_, dyn_model_->getInvM(), torque_HAND, N_HAND);
                }
                else {
                    task_HAND_->compute(r_HAND_diff, Kc, Dxi, J_r_HAND_6_, dq_, dyn_model_->getInvM(), torque_HAND, N_HAND);
                }

                torque_ = torque_JLC + N_JLC.transpose() * (torque_COL + (N_COL.transpose() * (torque_HAND)));


                // simulate one step
                Eigen::VectorXd prev_ddq(ddq_), prev_dq(dq_);
                dyn_model_->accel(ddq_, q_, dq_, torque_);
                float time_d = 0.005;
                for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                    dq_[q_idx] += (prev_ddq[q_idx] + ddq_[q_idx]) / 2.0 * time_d;
                    q_[q_idx] += (prev_dq[q_idx] + dq_[q_idx]) / 2.0 * time_d;
                }
    }

    void DynamicsSimulatorHandPose::oneStep() {
        KDL::Frame r_HAND_current;
        kin_model_->calculateFk(r_HAND_current, effector_name_, q_);
        KDL::Twist diff = KDL::diff(r_HAND_current, r_HAND_target_, 1.0);
        oneStep(diff);
    }

