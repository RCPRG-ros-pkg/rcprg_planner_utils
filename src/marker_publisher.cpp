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

#include "planer_utils/marker_publisher.h"
#include <iostream>
#include "geometry_msgs/Point.h"


MarkerPublisher::MarkerPublisher(ros::NodeHandle &nh) :
    nh_(nh)
{
    pub_ = nh_.advertise<visualization_msgs::MarkerArray>("/velma_markers", 1000);
}

MarkerPublisher::~MarkerPublisher() {
}

void MarkerPublisher::publish() {
    pub_.publish( marker_array_ );
    marker_array_.markers.clear();
}

void MarkerPublisher::clear() {
    marker_array_.markers.clear();
}
/*
int MarkerPublisher::addMarkerArray(int m_id, const visualization_msgs::MarkerArray &ma, const std::string &frame_id) {

    for (int i=0; i<ma.markers.size(); i++) {
        marker_array_.markers.push_back( ma.markers[i] );
        marker_array_.markers.back();
    }

}
*/

int MarkerPublisher::addLineListMarker(int m_id, const std::vector<KDL::Vector > &pts, const KDL::Frame &fr, double r, double g, double b, double a, double size, const std::string &frame_id) {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "default";
    marker.id = m_id;
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;
	marker.pose.position.x = fr.p.x();
	marker.pose.position.y = fr.p.y();
	marker.pose.position.z = fr.p.z();
	double qx, qy, qz, qw;
	fr.M.GetQuaternion(qx, qy, qz, qw);
	marker.pose.orientation.x = qx;
	marker.pose.orientation.y = qy;
	marker.pose.orientation.z = qz;
	marker.pose.orientation.w = qw;
    for (int i=0; i<pts.size(); i+=2) {
        geometry_msgs::Point p1, p2;
        p1.x = pts[i].x();
        p1.y = pts[i].y();
        p1.z = pts[i].z();
        p2.x = pts[i+1].x();
        p2.y = pts[i+1].y();
        p2.z = pts[i+1].z();
        marker.points.push_back(p1);
        marker.points.push_back(p2);
    }
    marker.scale.x = size;
    marker.scale.y = 0.0;
    marker.scale.z = 0.0;
    marker.color.a = a;
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker_array_.markers.push_back(marker);
    return m_id + 1;
}

int MarkerPublisher::addSinglePointMarker(int m_id, const KDL::Vector &pos, double r, double g, double b, double a, double size, const std::string &frame_id) {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "default";
    marker.id = m_id;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = pos.x();
    marker.pose.position.y = pos.y();
    marker.pose.position.z = pos.z();
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.pose.orientation.w = 1.0;
    marker.scale.x = size;
    marker.scale.y = size;
    marker.scale.z = size;
    marker.color.a = a;
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker_array_.markers.push_back(marker);
    return m_id + 1;
}

int MarkerPublisher::addSinglePointMarkerCube(int m_id, const KDL::Vector &pos, double r, double g, double b, double a, double size_x, double size_y, double size_z, const std::string &frame_id) {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "default";
    marker.id = m_id;
    marker.type = visualization_msgs::Marker::CUBE;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = pos.x();
    marker.pose.position.y = pos.y();
    marker.pose.position.z = pos.z();
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.pose.orientation.w = 1.0;
    marker.scale.x = size_x;
    marker.scale.y = size_y;
    marker.scale.z = size_z;
    marker.color.a = a;
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker_array_.markers.push_back(marker);
    return m_id + 1;
}

int MarkerPublisher::addMeshMarker(int m_id, const KDL::Vector &pos, double r, double g, double b, double a, double size_x, double size_y, double size_z, const std::string &mesh_path, const std::string &frame_id) {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "default";
    marker.id = m_id;
    marker.type = visualization_msgs::Marker::MESH_RESOURCE;
    marker.action = visualization_msgs::Marker::ADD;
    marker.mesh_resource = mesh_path;
    marker.pose.position.x = pos.x();
    marker.pose.position.y = pos.y();
    marker.pose.position.z = pos.z();
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.pose.orientation.w = 1.0;
    marker.scale.x = size_x;
    marker.scale.y = size_y;
    marker.scale.z = size_z;
    marker.color.a = a;
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker_array_.markers.push_back(marker);
    return m_id + 1;
}

int MarkerPublisher::addVectorMarker(int m_id, const KDL::Vector &v1, const KDL::Vector &v2, double r, double g, double b, double a, double size, const std::string &frame_id) {
    visualization_msgs::Marker marker;
    marker.header.frame_id = frame_id;
    marker.header.stamp = ros::Time();
    marker.ns = "default";
    marker.id = m_id;
    marker.type = visualization_msgs::Marker::ARROW;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.position.x = 0;
    marker.pose.position.y = 0;
    marker.pose.position.z = 0;
    marker.pose.orientation.x = 0.0;
    marker.pose.orientation.y = 0.0;
    marker.pose.orientation.z = 0.0;
    marker.pose.orientation.w = 1.0;
    geometry_msgs::Point p1, p2;
    p1.x = v1.x();
    p1.y = v1.y();
    p1.z = v1.z();
    p2.x = v2.x();
    p2.y = v2.y();
    p2.z = v2.z();
    marker.points.push_back(p1);
    marker.points.push_back(p2);
    marker.scale.x = size;
    marker.scale.y = 2.0 * size;
    marker.scale.z = 0.0;
    marker.color.a = a;
    marker.color.r = r;
    marker.color.g = g;
    marker.color.b = b;
    marker_array_.markers.push_back(marker);
    return m_id + 1;
}

int MarkerPublisher::addCapsule(int m_id, const KDL::Frame &fr, double r, double g, double b, double a, double length, double radius, const std::string &frame_id) {

	KDL::Vector zero;
	KDL::Vector v(0,0,length);
	KDL::Vector v2 = (fr * v) - (fr * zero);

	visualization_msgs::Marker marker;
	marker.header.frame_id = frame_id;
	marker.header.stamp = ros::Time();
	marker.ns = "default";
	marker.id = m_id;
	marker.type = visualization_msgs::Marker::SPHERE;
	marker.action = visualization_msgs::Marker::ADD;
	marker.pose.position.x = fr.p.x() - v2.x()/2.0;
	marker.pose.position.y = fr.p.y() - v2.y()/2.0;
	marker.pose.position.z = fr.p.z() - v2.z()/2.0;
	double qx, qy, qz, qw;
	fr.M.GetQuaternion(qx, qy, qz, qw);
	marker.pose.orientation.x = qx;
	marker.pose.orientation.y = qy;
	marker.pose.orientation.z = qz;
	marker.pose.orientation.w = qw;
	marker.scale.x = radius * 2.0;
	marker.scale.y = radius * 2.0;
	marker.scale.z = radius * 2.0;
	marker.color.a = a;
	marker.color.r = r;
	marker.color.g = g;
	marker.color.b = b;
	marker_array_.markers.push_back(marker);

	visualization_msgs::Marker marker2(marker);
	marker2.id = m_id+1;
	marker2.pose.position.x = fr.p.x() + v2.x()/2.0;
	marker2.pose.position.y = fr.p.y() + v2.y()/2.0;
	marker2.pose.position.z = fr.p.z() + v2.z()/2.0;
	if (length > 0.0001)
		marker_array_.markers.push_back(marker2);

	visualization_msgs::Marker marker3(marker);
	marker3.id = m_id+2;
	marker3.type = visualization_msgs::Marker::CYLINDER;
	marker3.pose.position.x = fr.p.x();
	marker3.pose.position.y = fr.p.y();
	marker3.pose.position.z = fr.p.z();
	marker3.scale.z = length;
//    KDL::Frame T_O_C(KDL::Rotation::RotX(90.0/180.0*3.1415));
    KDL::Frame T_B_C = fr;// * T_O_C;
	T_B_C.M.GetQuaternion(qx, qy, qz, qw);
	marker3.pose.orientation.x = qx;
	marker3.pose.orientation.y = qy;
	marker3.pose.orientation.z = qz;
	marker3.pose.orientation.w = qw;
	if (length > 0.0001)
		marker_array_.markers.push_back(marker3);

/*	visualization_msgs::Marker marker4(marker);
	marker4.id = m_id+3;
	marker4.type = visualization_msgs::Marker::SPHERE;
	marker4.pose.position.x = fr.p.x();
	marker4.pose.position.y = fr.p.y();
	marker4.pose.position.z = fr.p.z();
	marker4.scale.x = 0.01;
	marker4.scale.y = 0.01;
	marker4.scale.z = 0.01;
	marker4.color.a = 1.0;
	marker4.color.r = 0.0;
	marker4.color.g = 1.0;
	marker4.color.b = 0.0;
	marker_array_.markers.push_back(marker4);
*/
	return m_id + 3;
}

int MarkerPublisher::addOctomap(int m_id, const octomap::OcTree &om, const std::string &frame_id) {
//    visualization_msgs::MarkerArray occupiedNodesVis;
    // each array stores all cubes of a different size, one for each depth level:
    int base_idx = marker_array_.markers.size();
    marker_array_.markers.resize(marker_array_.markers.size() + om.getTreeDepth()+1);
//    occupiedNodesVis.markers.resize(om.getTreeDepth()+1);


  // now, traverse all leafs in the tree:
  for (octomap::OcTree::iterator it = om.begin(), end = om.end(); it != end; ++it) {
    if (om.isNodeOccupied(*it)){
      double z = it.getZ();
      double size = it.getSize();
      double x = it.getX();
      double y = it.getY();

      unsigned idx = it.getDepth();

      geometry_msgs::Point cubeCenter;
      cubeCenter.x = x;
      cubeCenter.y = y;
      cubeCenter.z = z;

      marker_array_.markers[idx + base_idx].points.push_back(cubeCenter);
      //occupiedNodesVis.markers[idx].points.push_back(cubeCenter);
      double minX, minY, minZ, maxX, maxY, maxZ;
      om.getMetricMin(minX, minY, minZ);
      om.getMetricMax(maxX, maxY, maxZ);

      double h = (1.0 - std::min(std::max((cubeCenter.z-minZ)/ (maxZ - minZ), 0.0), 1.0)) * 0.8;//m_colorFactor;
      marker_array_.markers[idx + base_idx].colors.push_back(heightMapColor(h));
      //occupiedNodesVis.markers[idx].colors.push_back(heightMapColor(h));
    }
  }

  // finish MarkerArray:
  for (unsigned i= 0; i < om.getTreeDepth()+1; ++i, ++m_id){
    double size = om.getNodeSize(i);

    marker_array_.markers[i + base_idx].header.frame_id = frame_id;
    marker_array_.markers[i + base_idx].header.stamp = ros::Time();
    marker_array_.markers[i + base_idx].ns = "default";
    marker_array_.markers[i + base_idx].id = m_id;
    marker_array_.markers[i + base_idx].type = visualization_msgs::Marker::CUBE_LIST;
    marker_array_.markers[i + base_idx].scale.x = size;
    marker_array_.markers[i + base_idx].scale.y = size;
    marker_array_.markers[i + base_idx].scale.z = size;
    if (marker_array_.markers[i + base_idx].points.size() > 0)
      marker_array_.markers[i + base_idx].action = visualization_msgs::Marker::ADD;
    else
      marker_array_.markers[i + base_idx].action = visualization_msgs::Marker::DELETE;

//    occupiedNodesVis.markers[i].header.frame_id = frame_id;
//    occupiedNodesVis.markers[i].header.stamp = ros::Time();
//    occupiedNodesVis.markers[i].ns = "default";
//    occupiedNodesVis.markers[i].id = m_id;
//    occupiedNodesVis.markers[i].type = visualization_msgs::Marker::CUBE_LIST;
//    occupiedNodesVis.markers[i].scale.x = size;
//    occupiedNodesVis.markers[i].scale.y = size;
//    occupiedNodesVis.markers[i].scale.z = size;
//    if (!m_useColoredMap)
//      occupiedNodesVis.markers[i].color = m_color;
//    if (occupiedNodesVis.markers[i].points.size() > 0)
//      occupiedNodesVis.markers[i].action = visualization_msgs::Marker::ADD;
//    else
//      occupiedNodesVis.markers[i].action = visualization_msgs::Marker::DELETE;
    }
    return m_id;
}

std_msgs::ColorRGBA MarkerPublisher::heightMapColor(double h) {

  std_msgs::ColorRGBA color;
  color.a = 1.0;
  // blend over HSV-values (more colors)

  double s = 1.0;
  double v = 1.0;

  h -= floor(h);
  h *= 6;
  int i;
  double m, n, f;

  i = floor(h);
  f = h - i;
  if (!(i & 1))
    f = 1 - f; // if i is even
  m = v * (1 - s);
  n = v * (1 - s * f);

  switch (i) {
    case 6:
    case 0:
      color.r = v; color.g = n; color.b = m;
      break;
    case 1:
      color.r = n; color.g = v; color.b = m;
      break;
    case 2:
      color.r = m; color.g = v; color.b = n;
      break;
    case 3:
      color.r = m; color.g = n; color.b = v;
      break;
    case 4:
      color.r = n; color.g = m; color.b = v;
      break;
    case 5:
      color.r = v; color.g = m; color.b = n;
      break;
    default:
      color.r = 1; color.g = 0.5; color.b = 0.5;
      break;
  }

  return color;
}

void MarkerPublisher::addEraseMarkers(int from, int to)
{
	for (int i=from; i<to; i++)
	{
		visualization_msgs::Marker marker;
		marker.header.frame_id = "world";
		marker.header.stamp = ros::Time();
		marker.ns = "default";
		marker.id = i;
		marker.type = visualization_msgs::Marker::SPHERE;
		marker.action = visualization_msgs::Marker::DELETE;
		marker_array_.markers.push_back(marker);
	}
}

