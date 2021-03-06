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
                        const std::vector<std::string > &joint_names, const Eigen::VectorXd &q_eq, const Eigen::VectorXd &max_q,
                        boost::function<KDL::Twist(const KDL::Frame &x, const KDL::Frame &y)> pose_diff_function) :
        ndof_(ndof),
        ndim_(ndim),
        effector_name_(effector_name),
        q_( Eigen::VectorXd::Zero(ndof_) ),
        dq_( Eigen::VectorXd::Zero(ndof_) ),
        ddq_( Eigen::VectorXd::Zero(ndof_) ),
        torque_( Eigen::VectorXd::Zero(ndof_) ),
        col_model_(col_model),
        kin_model_(kin_model),
        dyn_model_(dyn_model),
        links_fk_(col_model->getLinksCount()),
        activation_dist_(0.05),
        torque_JLC_( Eigen::VectorXd::Zero(ndof_) ),
        N_JLC_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        torque_WCCr_( Eigen::VectorXd::Zero(ndof_) ),
        N_WCCr_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        torque_WCCl_( Eigen::VectorXd::Zero(ndof_) ),
        N_WCCl_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        torque_COL_( Eigen::VectorXd::Zero(ndof_) ),
        N_COL_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        torque_HAND_( Eigen::VectorXd::Zero(ndof_) ),
        N_HAND_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        torque_JOINT_( Eigen::VectorXd::Zero(ndof_) ),
        N_JOINT_( Eigen::MatrixXd::Identity(ndof_, ndof_) ),
        Kc_(ndim_),
        Dxi_(ndim_),
        Kc_JOINT_(ndof_),
        Dxi_JOINT_(ndof_),
        r_HAND_diff_(ndim_),
        in_collision_(false),
        max_q_(max_q),
        q_eq_(q_eq),
        pose_diff_function_(pose_diff_function)
    {
        for (int q_idx = 0; q_idx < ndof; q_idx++) {
            Kc_JOINT_[q_idx] = 1.0;
            Dxi_JOINT_[q_idx] = 0.1;
        }

        double Kc_lin = 200.0;
        double Kc_rot = 20.0;
        if (ndim_ == 3) {
            Kc_[0] = Kc_lin;
            Kc_[1] = Kc_lin;
            Kc_[2] = Kc_rot;
        }
        else {
            for (int dim_idx = 0; dim_idx < 3; dim_idx++) {
                Kc_[dim_idx] = Kc_lin;
                Kc_[dim_idx+3] = Kc_rot;
            }
        }

        for (int dim_idx = 0; dim_idx < ndim_; dim_idx++) {
            Dxi_[dim_idx] = 0.7;
        }

        effector_idx_ = col_model->getLinkIndex(effector_name_);

        int wcc_r_q5_idx=-1, wcc_r_q6_idx=-1;
        int wcc_l_q5_idx=-1, wcc_l_q6_idx=-1;

        // joint limits
        Eigen::VectorXd lower_limit(ndof), upper_limit(ndof), limit_range(ndof), max_trq(ndof);
        int q_idx = 0;
        for (std::vector<std::string >::const_iterator name_it = joint_names.begin(); name_it != joint_names.end(); name_it++, q_idx++) {

            if ( (*name_it) == "right_arm_5_joint" ) {
                wcc_r_q5_idx = q_idx;
            }
            else if ( (*name_it) == "right_arm_6_joint" ) {
                wcc_r_q6_idx = q_idx;
            }
            else if ( (*name_it) == "left_arm_5_joint" ) {
                wcc_l_q5_idx = q_idx;
            }
            else if ( (*name_it) == "left_arm_6_joint" ) {
                wcc_l_q6_idx = q_idx;
            }

            lower_limit[q_idx] = kin_model->getLowerLimit(q_idx);
            upper_limit[q_idx] = kin_model->getUpperLimit(q_idx);
            limit_range[q_idx] = 10.0 / 180.0 * 3.141592653589793;
            max_trq[q_idx] = 10.0;
        }

        if (ndim_ != 3 && ndim_ != 6) {
            std::cout << "ERROR: DynamicsSimulatorHandPose works for 3 or 6 dimensions" << std::endl;
        }

        std::vector<double > wcc_polygon({
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
});
        std::vector<double > wcc_polygon_inv;
        for (int i = 0; i < wcc_polygon.size(); ++i) {
            wcc_polygon_inv.push_back(-1.0 * wcc_polygon[i]);
        }

        std::set<int > jlc_excluded_q_idx;
        if (wcc_r_q5_idx != -1 && wcc_r_q6_idx != -1) {
            jlc_excluded_q_idx.insert(wcc_r_q5_idx);
            jlc_excluded_q_idx.insert(wcc_r_q6_idx);
            task_WCCr_.reset( new Task_WCC(ndof, wcc_r_q5_idx, wcc_r_q6_idx, wcc_polygon) );
        }
        if (wcc_l_q5_idx != -1 && wcc_l_q6_idx != -1) {
            jlc_excluded_q_idx.insert(wcc_l_q5_idx);
            jlc_excluded_q_idx.insert(wcc_l_q6_idx);
            task_WCCl_.reset( new Task_WCC(ndof, wcc_l_q5_idx, wcc_l_q6_idx, wcc_polygon_inv) );
        }
        task_JLC_.reset( new Task_JLC(lower_limit, upper_limit, limit_range, max_trq, jlc_excluded_q_idx) );

        task_COL_.reset( new Task_COL(ndof_, activation_dist_, 80.0, kin_model_, col_model_) );
        task_HAND_.reset( new Task_HAND(ndof_, ndim_) );

        J_r_HAND_6_.resize(6, ndof_);
        J_r_HAND_.resize(ndim_, ndof_);

        task_JOINT_.reset( new Task_JOINT(ndof_) );
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

    void DynamicsSimulatorHandPose::oneStep(const KDL::Twist &diff, MarkerPublisher *markers_pub, int m_id) {
                // calculate forward kinematics for all links
                for (int l_idx = 0; l_idx < col_model_->getLinksCount(); l_idx++) {
                    kin_model_->calculateFk(links_fk_[l_idx], col_model_->getLinkName(l_idx), q_);
                }
                // calculate inertia matrix for whole body
                dyn_model_->computeM(q_);

                //
                // joint limit avoidance
                //
                task_JLC_->compute(q_, dq_, dyn_model_->getM(), torque_JLC_, N_JLC_);
                if (markers_pub != NULL) {
                    m_id = task_JLC_->visualize(markers_pub, m_id, kin_model_, col_model_, links_fk_);
                }

                in_collision_ = false;

                //
                // wrist collision constraint
                //
                if (task_WCCr_ != NULL) {
                    task_WCCr_->compute(q_, dq_, dyn_model_->getInvM(), torque_WCCr_, N_WCCr_);
                    //task_WCCr_->compute(q_, dq_, dyn_model_->getM(), dyn_model_->getInvM(), torque_WCCr_, N_WCCr_, markers_pub, &m_id);   // visual debug
                    if (task_WCCr_->inCollision()) {
                        //std::cout << "collision task_WCCr " << q_(6) << " " << q_(7) << std::endl;
                        in_collision_ = true;
                    }
                }

                if (task_WCCl_ != NULL) {
                    task_WCCl_->compute(q_, dq_, dyn_model_->getInvM(), torque_WCCl_, N_WCCl_);
                    if (task_WCCl_->inCollision()) {
                        //std::cout << "collision task_WCCl" << std::endl;
                        in_collision_ = true;
                    }
                }

                //
                // collision constraints
                //
                std::vector<self_collision::CollisionInfo> link_collisions;
                self_collision::getCollisionPairs(col_model_, links_fk_, activation_dist_, link_collisions);
                for (std::vector<self_collision::CollisionInfo>::const_iterator it = link_collisions.begin(); it != link_collisions.end(); it++) {
                    if ( it->dist <= 0.0 ) {
                        std::cout << "collision task_COL " << col_model_->getLinkName(it->link1_idx) << " " << col_model_->getLinkName(it->link2_idx) << std::endl;
                        in_collision_ = true;
                        break;
                    }
                }

                task_COL_->compute(q_, dq_, dyn_model_->getInvM(), links_fk_, N_JLC_ * N_WCCr_ * N_WCCl_, link_collisions, torque_COL_, N_COL_, markers_pub, m_id);

                //
                // effector task
                //
                if (ndim_ == 3) {
                    r_HAND_diff_[0] = diff[0];
                    r_HAND_diff_[1] = diff[1];
                    r_HAND_diff_[2] = diff[5];
                }
                else {
                    for (int dim_idx = 0; dim_idx < 6; dim_idx++) {
                        r_HAND_diff_[dim_idx] = diff[dim_idx];
                    }
                }

                kin_model_->getJacobian(J_r_HAND_6_, effector_name_, q_);

                if (ndim_ == 3) {
                    for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                        J_r_HAND_(0, q_idx) = J_r_HAND_6_(0, q_idx);
                        J_r_HAND_(1, q_idx) = J_r_HAND_6_(1, q_idx);
                        J_r_HAND_(2, q_idx) = J_r_HAND_6_(5, q_idx);
                    }
                    task_HAND_->compute(r_HAND_diff_, Kc_, Dxi_, J_r_HAND_, dq_, dyn_model_->getInvM(), torque_HAND_, N_HAND_);
                }
                else {
//                    J_r_HAND_6_ = J_r_HAND_6_ * N_JLC_ * N_WCCr_ * N_WCCl_ * N_COL_;
//                    Eigen::MatrixXd j = J_r_HAND_6_ * N_JLC_ * N_WCCr_ * N_WCCl_ * N_COL_;
//                    Eigen::FullPivLU< Eigen::MatrixXd > lu(j);
//                    lu.setThreshold(0.1);
//                    if (lu.rank() < 6) {
//                        in_collision_ = true;
//                    }
//                    std::cout << lu.rank() << std::endl;
                    task_HAND_->compute(r_HAND_diff_, Kc_, Dxi_, J_r_HAND_6_, dq_, dyn_model_->getInvM(), torque_HAND_, N_HAND_);
                }

                //
                // joint impedance
                //
                Eigen::VectorXd q_diff = q_eq_ - q_;
                task_JOINT_->compute(q_diff, Kc_JOINT_, Dxi_JOINT_, dq_, dyn_model_->getM(), torque_JOINT_, N_JOINT_);

//                torque_.noalias() = torque_JLC_;// + N_JLC_.transpose() * (torque_COL_ + (N_COL_.transpose() * (torque_HAND_)));
//                torque_.noalias() = (torque_JLC_ + torque_WCCr_ + torque_WCCl_) + (N_JLC_ * N_WCCr_ * N_WCCl_).transpose() * (torque_COL_ + (N_COL_.transpose() * (torque_HAND_)));
                torque_.noalias() = (torque_JLC_ + torque_WCCr_ + torque_WCCl_) + (N_JLC_ * N_WCCr_ * N_WCCl_).transpose() * (torque_COL_ + (N_COL_.transpose() * (torque_HAND_ + N_HAND_.transpose() * (torque_JOINT_))));
//                torque_.noalias() = (torque_JLC_ + torque_WCCr_ + torque_WCCl_) + (N_JLC_ * N_WCCr_ * N_WCCl_) * (torque_COL_ + (N_COL_ * (torque_HAND_ + N_HAND_ * (torque_JOINT_))));
//                torque_.noalias() = torque_JLC_ + N_JLC_.transpose() * (torque_WCCr_ + N_WCCr_.transpose() * (torque_WCCl_ + N_WCCl_.transpose() * (torque_COL_ + (N_COL_.transpose() * (torque_HAND_)))));
//                torque_.noalias() = (torque_JLC_ + torque_WCCr_ + torque_WCCl_) + (N_JLC_ * N_WCCr_ * N_WCCl_).transpose() * (torque_HAND_ + N_HAND_ * (torque_JOINT_));


//                std::cout << torque_COL_.transpose() << std::endl;
//                std::cout << "JLC: " << task_JLC_->getActivationCount() << " WCCr: " << task_WCCr_->getActivationCount() << " WWCl: " << task_WCCl_->getActivationCount() << " COL: " << task_COL_->getActivationCount() << std::endl;

                // simulate one step
                Eigen::VectorXd prev_ddq(ddq_), prev_dq(dq_);
                dyn_model_->accel(ddq_, q_, dq_, torque_);
//                float time_d = 0.005;
//                float time_d = 0.005;
                float time_d = 0.001;

                double max_f = 1.0;
                for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                    dq_[q_idx] += (prev_ddq[q_idx] + ddq_[q_idx]) / 2.0 * time_d;
                    double f = fabs(dq_(q_idx) / max_q_(q_idx));
                    if (max_f < f) {
                        max_f = f;
                    }
                }

                KDL::Vector p(r_HAND_target_.p);
                double qx, qy, qz, qw;
                r_HAND_target_.M.GetQuaternion(qx,qy,qz,qw);

                if (max_f > 1.0) {
                    dq_ /= max_f;
                }

                q_ += (prev_dq + dq_) / 2.0 * time_d;
                for (int q_idx = 0; q_idx < ndof_; q_idx++) {
                    if (q_[q_idx] != q_[q_idx]) {
                        std::cout << "q_: " << q_.transpose() << std::endl;
                        std::cout << "dq_: " << dq_.transpose() << std::endl;
                        std::cout << "torque_: " << torque_.transpose() << std::endl;
                        std::cout << "torque_JLC_: " << torque_JLC_.transpose() << std::endl;
                        std::cout << "torque_WCCr_: " << torque_WCCr_.transpose() << std::endl;
                        std::cout << "torque_WCCl_: " << torque_WCCl_.transpose() << std::endl;
                        std::cout << "torque_COL_: " << torque_COL_.transpose() << std::endl;
                        std::cout << "torque_HAND_: " << torque_HAND_.transpose() << std::endl;
                        std::cout << "torque_JOINT_: " << torque_JOINT_.transpose() << std::endl;
                        std::cout << "N_JLC_: " << std::endl << N_JLC_ << std::endl;
                        std::cout << "N_WCCr_: " << std::endl << N_WCCr_ << std::endl;
                        std::cout << "N_WCCl_: " << std::endl << N_WCCl_ << std::endl;
                        std::cout << "N_COL_: " << std::endl << N_COL_ << std::endl;
                        std::cout << "N_HAND_: " << std::endl << N_HAND_ << std::endl;

                        break;
                    }
                }
    }

    void DynamicsSimulatorHandPose::oneStep(MarkerPublisher *markers_pub, int m_id) {
        KDL::Frame r_HAND_current;
        kin_model_->calculateFk(r_HAND_current, effector_name_, q_);
//        KDL::Twist diff = KDL::diff(r_HAND_current, r_HAND_target_, 1.0);
        KDL::Twist diff = pose_diff_function_(r_HAND_current, r_HAND_target_);
        oneStep(diff, markers_pub, m_id);
    }

    bool DynamicsSimulatorHandPose::inCollision() const {
        return in_collision_;
    }

    void DynamicsSimulatorHandPose::updateMetric(boost::function<KDL::Twist(const KDL::Frame &, const KDL::Frame &)> pose_diff_function) {
        pose_diff_function_ = pose_diff_function;
    }

