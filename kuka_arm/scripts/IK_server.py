#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np
import csv




def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        data = [] #Initilaize list where error data is stored
        # Define DH param symbols
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')



        # Joint angle symbols
        # DH Parameters
        s = {alpha0: 0,      a0:   0,    d1: 0.75, 
             alpha1: -pi/2,  a1: 0.35,   d2: 0,     q2: q2- pi/2,
             alpha2: 0,      a2: 1.25,   d3: 0,
             alpha3: -pi/2,  a3: -0.054, d4: 1.5,
             alpha4: pi/2,   a4: 0   ,   d5: 0,
             alpha5: -pi/2,  a5: 0   ,   d6: 0,
             alpha6: 0,      a6: 0   ,   d7: 0.303, q7: 0}


        # Create individual transformation matrices
        T0_1 = Matrix([[             cos(q1),            -sin(q1),            0,              a0],
                       [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
                       [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),  cos(alpha0),  cos(alpha0)*d1],
                       [                   0,                   0,            0,               1]])
        T0_1 = T0_1.subs(s)

        T1_2 = Matrix([[             cos(q2),            -sin(q2),            0,              a1],
                       [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
                       [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),  cos(alpha1),  cos(alpha1)*d2],
                       [                   0,                   0,            0,               1]])
        T1_2 = T1_2.subs(s)

        T2_3 = Matrix([[             cos(q3),            -sin(q3),            0,              a2],
                       [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
                       [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),  cos(alpha2),  cos(alpha2)*d3],
                       [                   0,                   0,            0,               1]])
        T2_3 = T2_3.subs(s)

        T3_4 = Matrix([[             cos(q4),            -sin(q4),            0,              a3],
                       [ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
                       [ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3),  cos(alpha3),  cos(alpha3)*d4],
                       [                   0,                   0,            0,               1]])
        T3_4 = T3_4.subs(s)


        T4_5 = Matrix([[             cos(q5),            -sin(q5),            0,              a4],
                       [ sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
                       [ sin(q5)*sin(alpha4), cos(q5)*sin(alpha4),  cos(alpha4),  cos(alpha4)*d5],
                       [                   0,                   0,            0,               1]])
        T4_5 = T4_5.subs(s)


        T5_6 = Matrix([[             cos(q6),            -sin(q6),            0,              a5],
                       [ sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
                       [ sin(q6)*sin(alpha5), cos(q6)*sin(alpha5),  cos(alpha5),  cos(alpha5)*d6],
                       [                   0,                   0,            0,               1]])
        T5_6 = T5_6.subs(s)

        T6_G = Matrix([[             cos(q7),            -sin(q7),            0,              a6],
                       [ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
                       [ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6)*d7],
                       [                   0,                   0,            0,               1]])
        T6_G = T6_G.subs(s)


        

        
        for x in xrange(0, len(req.poses)):
            joint_trajectory_point = JointTrajectoryPoint()
            
            # Extract end-effector position and orientation from request
	        # px,py,pz = end-effector position
	        # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z
            
            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
            
            P_EE = Matrix([[px],[py],[pz]])
            
            # Get rotation matrix from roll , pitch and yaw angles      
            R_roll = Matrix([[ 1,              0,        0],
                                  [ 0,        cos(roll), -sin(roll)],
                                  [ 0,        sin(roll),  cos(roll)]])

            R_pitch = Matrix([[ cos(pitch),        0,  sin(pitch)],
                                      [       0,        1,        0],
                                      [-sin(pitch),        0,  cos(pitch)]])

            R_yaw = Matrix([[ cos(yaw), -sin(yaw),        0],
                                      [ sin(yaw),  cos(yaw),        0],
                                      [ 0,              0,        1]])

            Rrpy = simplify( R_yaw*R_pitch *R_roll)
            
            # Calculate wrist center from EE position
            P_WC = simplify(P_EE - 0.303*Rrpy * Matrix([[1],[0],[0]]))
            
            # Calculate theta1
            theta1 = atan2(P_WC[1],P_WC[0])
            
            # Calculate theta 3
            T0_2 = T0_1 * T1_2
            
            j2 = (T0_2.evalf(subs={q1: theta1}))[0:3,3]
            
            dist_j2_wc = sqrt((j2[0]-P_WC[0])**2 + (j2[1]-P_WC[1])**2 + (j2[2]-P_WC[2])**2)
            dist_j3_wc = sqrt(s[a3]**2 + s[d4]**2)
            dist_j2_j3 = s[a2]
            
            j3_5_angle_offset = atan2(-0.054,1.5)
            
            cos_theta3_2 = (dist_j2_wc**2 - dist_j2_j3**2 - dist_j3_wc**2)/(2*dist_j2_j3*dist_j3_wc)
            theta3_2 = atan2(sqrt(1-(cos_theta3_2)**2),cos_theta3_2)
            
            theta3 = -(pi/2) + j3_5_angle_offset + theta3_2
            
            #Calculate theta 2
            theta2_a = atan2(P_WC[2]-j2[2],sqrt((j2[0]-P_WC[0])**2 + (j2[1]-P_WC[1])**2))
            cos_theta2_b = (-dist_j3_wc**2 + dist_j2_j3**2 + dist_j2_wc**2)/(2*dist_j2_j3*dist_j2_wc)
            theta2_b = atan2(sqrt(1-(cos_theta2_b)**2),cos_theta2_b)
            theta2 = pi/2-theta2_a-theta2_b
            
            # Obtain transform from 3 to 6 after correction to solve for theta 4, 5 ,6.
            T0_3 = (T0_1 * T1_2 * T2_3).evalf(subs={q1: theta1,q2: theta2,q3:theta3})
            
            R_x = Matrix([[1, 0, 0],
                    [0, cos(pi/2), -sin(pi/2)],
                    [0, sin(pi/2), cos(pi/2)]])
            R_z = Matrix([[cos(pi/2), -sin(pi/2), 0],
                    [sin(pi/2), cos(pi/2), 0],
                    [0, 0, 1]])
            R3_6 = T0_3.T[0:3,0:3] * Rrpy*(R_z * R_x).T
            (theta4,theta5,theta6) = tf.transformations.euler_from_matrix(np.array(R3_6).astype(np.float64), "ryzy")
            # Populate response for the IK request
            
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
            
            # Compute error in pose by forward kinematics calculations
            R_y_total = Matrix([[cos(-pi/2),   0, sin(-pi/2)],[      0,      1,         0],[-sin(-pi/2),  0, cos(-pi/2)]])
            R_z_total = Matrix([[cos(pi), -sin(pi), 0],[sin(pi), cos(pi), 0],[0, 0, 1]])
            R_correction = R_z_total*R_y_total
            T0_G = (T0_1* T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G ).evalf(subs={q1: theta1,q2: theta2,q3:theta3,q4: theta4,q5: theta5,q6:theta6})
            (fk_roll, fk_pitch, fk_yaw) =tf.transformations.euler_from_matrix(np.array(T0_G[:3,:3]*R_correction))
            (fk_x,fk_y,fk_z) = (T0_G[0,3],T0_G[1,3],T0_G[2,3])
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
            data.append([P_EE[0]-T0_G[0,3],P_EE[1]-T0_G[1,3],P_EE[2]-T0_G[2,3],roll-fk_roll,pitch-fk_pitch,yaw-fk_yaw])

            print("Thetas:",theta1, theta2, theta3, theta4, theta5, theta6)
        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        # Dump errors in csv file for accuracy checks or debug
        with open('data.csv', "a") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(data)
        
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
