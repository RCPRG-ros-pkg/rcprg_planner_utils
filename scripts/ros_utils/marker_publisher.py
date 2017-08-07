#!/usr/bin/env python

# Copyright (c) 2015, Robot Control and Pattern Recognition Group,
# Institute of Control and Computation Engineering
# Warsaw University of Technology
#
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Warsaw University of Technology nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYright HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Author: Dawid Seredynski
#

import rospy
#import tf

from std_msgs.msg import *
from sensor_msgs.msg import *
from geometry_msgs.msg import *
from visualization_msgs.msg import *

#from tf import *
#from tf.transformations import * 
#from tf2_msgs.msg import *

import PyKDL
import copy

class MarkerPublisher:
    def __init__(self, namespace):
        self.pub_marker = rospy.Publisher(namespace, MarkerArray, queue_size=1000)

    def publishSinglePointMarker(self, pt, i, r=1, g=0, b=0, a=1, namespace='default', frame_id='torso_base', m_type=Marker.CUBE, scale=Vector3(0.005, 0.005, 0.005), T=None):
        m = MarkerArray()
        marker = Marker()
        marker.header.frame_id = frame_id
        marker.header.stamp = rospy.Time.now()
        marker.ns = namespace
        marker.id = i
        marker.type = m_type
        marker.action = Marker.ADD
        if T != None:
            point = T*pt
            q = T.M.GetQuaternion()
            marker.pose = Pose( Point(point.x(),point.y(),point.z()), Quaternion(q[0],q[1],q[2],q[3]) )
        else:
            marker.pose = Pose( Point(pt.x(),pt.y(),pt.z()), Quaternion(0,0,0,1) )
        marker.scale = scale
        marker.color = ColorRGBA(r,g,b,a)
        m.markers.append(marker)
        self.pub_marker.publish(m)
        return i+1

    def eraseMarkers(self, idx_from, idx_to, frame_id='torso_base', namespace='default'):
        m = MarkerArray()
        for idx in range(idx_from, idx_to):
            marker = Marker()
            marker.header.frame_id = frame_id
            marker.header.stamp = rospy.Time.now()
            marker.ns = namespace
            marker.id = idx
            marker.action = Marker.DELETE
            m.markers.append(marker)
        if len(m.markers) > 0:
            self.pub_marker.publish(m)


    def publishMultiPointsMarker(self, pt, base_id, r=1, g=0, b=0, namespace='default', frame_id='torso_base', m_type=Marker.CUBE, scale=Vector3(0.002, 0.002, 0.002), T=None):
        m = MarkerArray()
        ret_id = copy.copy(base_id)
        for i in range(0, len(pt)):
            marker = Marker()
            marker.header.frame_id = frame_id
            marker.header.stamp = rospy.Time.now()
            marker.ns = namespace
            marker.id = ret_id
            ret_id += 1
            marker.type = m_type
            marker.action = Marker.ADD
            if T != None:
                point = T*pt[i]
                marker.pose = Pose( Point(point.x(),point.y(),point.z()), Quaternion(0,0,0,1) )
            else:
                marker.pose = Pose( Point(pt[i].x(),pt[i].y(),pt[i].z()), Quaternion(0,0,0,1) )
            marker.scale = scale
            marker.color = ColorRGBA(r,g,b,0.5)
            m.markers.append(marker)
        self.pub_marker.publish(m)
        return ret_id

    def publishMultiPointsMarkerWithSize(self, pt, base_id, r=1, g=0, b=0, namespace='default', frame_id='torso_base', m_type=Marker.SPHERE, T=None):
        m = MarkerArray()
        ret_id = copy.copy(base_id)
        if T == None:
            T = PyKDL.Frame()
        for i in range(0, len(pt)):
            marker = Marker()
            marker.header.frame_id = frame_id
            marker.header.stamp = rospy.Time.now()
            marker.ns = namespace
            marker.id = ret_id
            ret_id += 1
            marker.type = m_type
            marker.action = Marker.ADD
            point = T*pt[i][0]
            marker.pose = Pose( Point(point.x(),point.y(),point.z()), Quaternion(0,0,0,1) )
            marker.scale = Vector3(pt[i][1], pt[i][1], pt[i][1])
            marker.color = ColorRGBA(r,g,b,0.5)
            m.markers.append(marker)
        self.pub_marker.publish(m)
        return ret_id

    def publishVectorMarker(self, v1, v2, i, r, g, b, frame='torso_base', namespace='default', scale=0.001):
        m = MarkerArray()
        marker = Marker()
        marker.header.frame_id = frame
        marker.header.stamp = rospy.Time.now()
        marker.ns = namespace
        marker.id = i
        marker.type = Marker.ARROW
        marker.action = Marker.ADD
        marker.points.append(Point(v1.x(), v1.y(), v1.z()))
        marker.points.append(Point(v2.x(), v2.y(), v2.z()))
        marker.pose = Pose( Point(0,0,0), Quaternion(0,0,0,1) )
        marker.scale = Vector3(scale, 2.0*scale, 0)
        marker.color = ColorRGBA(r,g,b,0.5)
        m.markers.append(marker)
        self.pub_marker.publish(m)
        return i+1

    def publishFrameMarker(self, T, base_id, scale=0.1, frame='torso_base', namespace='default'):
        self.publishVectorMarker(T*PyKDL.Vector(), T*PyKDL.Vector(scale,0,0), base_id, 1, 0, 0, frame, namespace)
        self.publishVectorMarker(T*PyKDL.Vector(), T*PyKDL.Vector(0,scale,0), base_id+1, 0, 1, 0, frame, namespace)
        self.publishVectorMarker(T*PyKDL.Vector(), T*PyKDL.Vector(0,0,scale), base_id+2, 0, 0, 1, frame, namespace)
        return base_id+3

    def publishConstantMeshMarker(self, uri, base_id, r=1, g=0, b=0, scale=0.1, frame_id='torso_base', namespace='default', T=None):
        if T == None:
            T = PyKDL.Frame()
        m = MarkerArray()
        marker = Marker()
        marker.header.frame_id = frame_id
        marker.header.stamp = rospy.Time.now()
        marker.ns = namespace
        marker.id = base_id
        marker.type = Marker.MESH_RESOURCE
        marker.mesh_resource = uri          # e.g. "package://pr2_description/meshes/base_v0/base.dae"
        marker.action = Marker.ADD
        point = T.p
        q = T.M.GetQuaternion()
        marker.pose = Pose( Point(point.x(),point.y(),point.z()), Quaternion(q[0],q[1],q[2],q[3]) )
        marker.scale = Vector3(scale, scale, scale)
        marker.color = ColorRGBA(r,g,b,0.5)
        m.markers.append(marker)
        self.pub_marker.publish(m)
        return base_id+1

    def publishMeshMarker(self, mesh, base_id, r=1, g=0, b=0, scale=0.1, frame_id='torso_base', namespace='default', T=None):
        m = MarkerArray()
        marker = Marker()
        marker.header.frame_id = frame_id
        marker.header.stamp = rospy.Time.now()
        marker.ns = namespace
        marker.id = base_id
        marker.type = Marker.TRIANGLE_LIST
        marker.action = Marker.ADD
        for f in mesh[1]:
            marker.points.append(Point(mesh[0][f[0]][0], mesh[0][f[0]][1], mesh[0][f[0]][2]))
            marker.points.append(Point(mesh[0][f[1]][0], mesh[0][f[1]][1], mesh[0][f[1]][2]))
            marker.points.append(Point(mesh[0][f[2]][0], mesh[0][f[2]][1], mesh[0][f[2]][2]))
        if T != None:
            point = T.p
            q = T.M.GetQuaternion()
            marker.pose = Pose( Point(point.x(),point.y(),point.z()), Quaternion(q[0],q[1],q[2],q[3]) )
        else:
            marker.pose = Pose( Point(0,0,0), Quaternion(0,0,0,1) )
        marker.scale = Vector3(1.0, 1.0, 1.0)
        marker.color = ColorRGBA(r,g,b,0.5)
        m.markers.append(marker)
        self.pub_marker.publish(m)
        return base_id+1

