#!/usr/bin/env python
import roslib; roslib.load_manifest('planer_utils')

import sys
import rospy
import math
import copy

from std_msgs.msg import ColorRGBA
from interactive_markers.interactive_marker_server import *
from interactive_markers.menu_handler import *
from visualization_msgs.msg import *
from geometry_msgs.msg import *
from tf.transformations import * 
import tf_conversions.posemath as pm
import PyKDL

class PoseIntMarker:
    def __init__(self):
        # create an interactive marker server on the topic namespace simple_marker
        self.server = InteractiveMarkerServer('pose_int_marker')

        self.insert6DofGlobalMarker()

       
        self.server.applyChanges();

    def insert6DofGlobalMarker(self):
        int_position_marker = InteractiveMarker()
        int_position_marker.header.frame_id = 'world'
        int_position_marker.name = 'pose_int_marker'
        int_position_marker.scale = 0.2
        int_position_marker.pose = pm.toMsg(PyKDL.Frame(PyKDL.Vector(1,0,1.5)))

        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.ROTATE_AXIS,'x'));
        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.ROTATE_AXIS,'y'));
        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.ROTATE_AXIS,'z'));
        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.MOVE_AXIS,'x'));
        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.MOVE_AXIS,'y'));
        int_position_marker.controls.append(self.createInteractiveMarkerControl6DOF(InteractiveMarkerControl.MOVE_AXIS,'z'));

        box = self.createAxisMarkerControl(Point(0.15,0.015,0.015), Point(0.0, 0.0, 0.0) )
        box.interaction_mode = InteractiveMarkerControl.BUTTON
        box.name = 'button'
        int_position_marker.controls.append( box )
        self.server.insert(int_position_marker, self.processFeedback);
        self.server.applyChanges();

    def processFeedback(self, feedback):
        #print "feedback", feedback.marker_name, feedback.control_name, feedback.event_type

        if ( feedback.marker_name == 'pose_int_marker' ):
            T_W_M = pm.fromMsg(feedback.pose)
            qx,qy,qz,qw = T_W_M.M.GetQuaternion()
            #print "position:(", T_W_M.p.x(), ",", T_W_M.p.y(), ",", T_W_M.p.z(), ")  orientation: (", qx, ",", qy, ",", qz, ",", qw, ")"
            print "PyKDL.Frame(PyKDL.Rotation.Quaternion(", qx, ",", qy, ",", qz, ",", qw, "), PyKDL.Vector(", T_W_M.p.x(), ",", T_W_M.p.y(), ",", T_W_M.p.z(), "))"

    def createAxisMarkerControl(self, scale, position):
        markerX = Marker()
        markerX.type = Marker.ARROW
        markerX.scale = scale
        markerX.pose.position = position
        ori = quaternion_about_axis(0, [0, 1 ,0])
        markerX.pose.orientation = Quaternion(ori[0], ori[1], ori[2], ori[3])
        markerX.color = ColorRGBA(1,0,0,1)
        markerY = Marker()
        markerY.type = Marker.ARROW
        markerY.scale = scale
        markerY.pose.position = position
        ori = quaternion_about_axis(math.pi/2.0, [0, 0 ,1])
        markerY.pose.orientation = Quaternion(ori[0], ori[1], ori[2], ori[3])
        markerY.color = ColorRGBA(0,1,0,1)
        markerZ = Marker()
        markerZ.type = Marker.ARROW
        markerZ.scale = scale
        markerZ.pose.position = position
        ori = quaternion_about_axis(-math.pi/2.0, [0, 1 ,0])
        markerZ.pose.orientation = Quaternion(ori[0], ori[1], ori[2], ori[3])
        markerZ.color = ColorRGBA(0,0,1,1)
        control = InteractiveMarkerControl()
        control.always_visible = True;
        control.markers.append( markerX );
        control.markers.append( markerY );
        control.markers.append( markerZ );
        return control

    def createButtoMarkerControl(self, scale, position, color):
        marker = Marker()
        marker.type = Marker.SPHERE
        marker.scale = scale
        marker.pose.position = position
        marker.pose.orientation = Quaternion(0,0,0,1)
        marker.color = color
        control = InteractiveMarkerControl()
        control.always_visible = True;
        control.markers.append( marker );
        return control

    def createInteractiveMarkerControl6DOF(self, mode, axis):
        control = InteractiveMarkerControl()
        control.orientation_mode = InteractiveMarkerControl.FIXED
        if mode == InteractiveMarkerControl.ROTATE_AXIS:
            control.name = 'rotate_';
        if mode == InteractiveMarkerControl.MOVE_AXIS:
            control.name = 'move_';
        if axis == 'x':
            control.orientation = Quaternion(1,0,0,1)
            control.name = control.name+'x';
        if axis == 'y':
            control.orientation = Quaternion(0,1,0,1)
            control.name = control.name+'x';
        if axis == 'z':
            control.orientation = Quaternion(0,0,1,1)
            control.name = control.name+'x';
        control.interaction_mode = mode
        return control

if __name__ == "__main__":

    rospy.init_node('pose_int_marker', anonymous=True)

    rospy.sleep(1)

    int_marker = PoseIntMarker()

    rospy.spin()

