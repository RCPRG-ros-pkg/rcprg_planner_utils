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

#ifndef SIMULATOR_H__
#define SIMULATOR_H__

#include <boost/bind.hpp>

#include <string>
#include <stdlib.h>
#include <stdio.h>

#include "Eigen/Dense"

#include "kin_dyn_model/dyn_model.h"
#include <collision_convex_model/collision_convex_model.h>
#include "kin_dyn_model/kin_model.h"
#include "planer_utils/task_col.h"
#include "planer_utils/task_hand.h"
#include "planer_utils/task_jlc.h"
#include "planer_utils/task_wcc.h"
#include "planer_utils/task_joint.h"
#include "planer_utils/marker_publisher.h"

class DynamicsSimulatorHandPose {
private:
    DynamicsSimulatorHandPose(const DynamicsSimulatorHandPose&);
    DynamicsSimulatorHandPose &operator=(const DynamicsSimulatorHandPose&);

protected:
    int ndof_;
    int ndim_;
    int effector_idx_;
    Eigen::VectorXd q_, dq_, ddq_, torque_;
    const boost::shared_ptr<self_collision::CollisionModel> &col_model_;
    const boost::shared_ptr<KinematicModel> &kin_model_;
    const boost::shared_ptr<DynamicModel > &dyn_model_;
    double activation_dist_;
    boost::shared_ptr<Task_JLC> task_JLC_;
    boost::shared_ptr<Task_COL> task_COL_;
    boost::shared_ptr<Task_HAND> task_HAND_;
    boost::shared_ptr<Task_WCC> task_WCCr_;
    boost::shared_ptr<Task_WCC> task_WCCl_;
    boost::shared_ptr<Task_JOINT> task_JOINT_;

    Eigen::VectorXd torque_JLC_;
    Eigen::MatrixXd N_JLC_;
    Eigen::VectorXd torque_WCCr_;
    Eigen::MatrixXd N_WCCr_;
    Eigen::VectorXd torque_WCCl_;
    Eigen::MatrixXd N_WCCl_;
    Eigen::VectorXd torque_COL_;
    Eigen::MatrixXd N_COL_;
    Eigen::VectorXd torque_HAND_;
    Eigen::MatrixXd N_HAND_;
    Eigen::VectorXd torque_JOINT_;
    Eigen::MatrixXd N_JOINT_;

    Eigen::VectorXd Kc_;
    Eigen::VectorXd Dxi_;
    Eigen::VectorXd Kc_JOINT_;
    Eigen::VectorXd Dxi_JOINT_;
    Eigen::VectorXd r_HAND_diff_;

    bool in_collision_;
    Eigen::VectorXd max_q_;
    Eigen::VectorXd q_eq_;

    std::vector<KDL::Frame > links_fk_;
    KDL::Frame r_HAND_target_;
    KinematicModel::Jacobian J_r_HAND_6_, J_r_HAND_;
    std::string effector_name_;

    boost::function<KDL::Twist(const KDL::Frame &, const KDL::Frame &)> pose_diff_function_;

public:

    DynamicsSimulatorHandPose(int ndof, int ndim, const std::string &effector_name, const boost::shared_ptr<self_collision::CollisionModel> &col_model,
                        const boost::shared_ptr<KinematicModel> &kin_model,
                        const boost::shared_ptr<DynamicModel > &dyn_model,
                        const std::vector<std::string > &joint_names, const Eigen::VectorXd &q_eq, const Eigen::VectorXd &max_q,
                        boost::function<KDL::Twist(const KDL::Frame &, const KDL::Frame &)> pose_diff_function = boost::bind(static_cast<KDL::Twist (*)(const KDL::Frame &, const KDL::Frame &, double)>(&KDL::diff), _1, _2, 1.0));

    void setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq);
    void getState(Eigen::VectorXd &q, Eigen::VectorXd &dq, Eigen::VectorXd &ddq);
    void setTarget(const KDL::Frame &r_HAND_target);
    void oneStep(const KDL::Twist &diff, MarkerPublisher *markers_pub=NULL, int m_id=0);
    void oneStep(MarkerPublisher *markers_pub=NULL, int m_id=0);
    bool inCollision() const;

    void updateMetric(boost::function<KDL::Twist(const KDL::Frame &, const KDL::Frame &)> pose_diff_function);
};

#endif  // SIMULATOR_H__

