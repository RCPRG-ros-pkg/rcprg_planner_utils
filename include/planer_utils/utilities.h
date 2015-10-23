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

#ifndef UTILITIES_H__
#define UTILITIES_H__

#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include <boost/bind.hpp>

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "Eigen/Dense"

#include <collision_convex_model/collision_convex_model.h>
#include "marker_publisher.h"

void publishJointState(ros::Publisher &joint_state_pub, const Eigen::VectorXd &q, const std::vector<std::string > &joint_names);
void publishJointState(ros::Publisher &joint_state_pub, const Eigen::VectorXd &q, const std::vector<std::string > &joint_names, const Eigen::VectorXd &ign_q, const std::vector<std::string > &ign_joint_names);
void publishTransform(tf::TransformBroadcaster &br, const KDL::Frame &T_B_F, const std::string &frame_id, const std::string &base_frame_id);
int addRobotModelVis(MarkerPublisher &markers_pub, int m_id, const boost::shared_ptr<self_collision::CollisionModel> &col_model, const std::vector<KDL::Frame > &T);
void getPointOnPath(const std::list<Eigen::VectorXd > &path, double f, Eigen::VectorXd &x);
double getPathLength(const std::list<Eigen::VectorXd > &path);
void printFrameKDL(const KDL::Frame &f);
double getAngle(const KDL::Vector &v1, const KDL::Vector &v2);
void EigenTfToKDL(const Eigen::Isometry3d &tf, KDL::Frame &kdlT);
void KDLToEigenTf(const KDL::Frame &kdlT, Eigen::Isometry3d &tf);
std::ostream& operator<< (std::ostream& stream, const KDL::Frame& f);
std::istream& operator>> (std::istream& stream, KDL::Frame& f);
std::string frameKdl2string(const KDL::Frame &f);
KDL::Frame string2frameKdl(const std::string &str);
std::string double2string(double d);
double string2double(const std::string &str);

double triVariateIsotropicGaussianKernel(const Eigen::Vector3d &x, const Eigen::Vector3d &mean, double sigma);
double biVariateIsotropicGaussianKernel(const Eigen::Vector2d &x, const Eigen::Vector2d &mean, double sigma);
double uniVariateIsotropicGaussianKernel(double x, double mean, double sigma);
double misesFisherKernelConstant(double sigma, int dimensions);
double misesFisherKernel(const Eigen::Vector4d &q, const Eigen::Vector4d &mean, double sigma, double Cp);
double misesFisherKernel(const Eigen::Vector3d &q, const Eigen::Vector3d &mean, double sigma, double Cp);
int vonMisesFisherSample(const Eigen::Vector3d &mean, double pdf_mean, double sigma, double Cp, Eigen::Vector3d &x);
int vonMisesFisherSample(const Eigen::Vector4d &mean, double pdf_mean, double sigma, double Cp, Eigen::Vector4d &x);

#endif  // UTILITIES_H__

