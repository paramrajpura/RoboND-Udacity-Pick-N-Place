## Project: Kinematics Pick & Place
### Objective: Compute joint angles from given end-effector pose (position and orientation) by solving Inverse Kinematics using Geometry. 

---


**Approach considered for the project:**  


1. Using the kr210.urdf.xacro file, analyze Kuka KR210 and derive DH parameters.
2. Define the transform matrices for each joint starting from base to gripper link using DH parameters.  
3. Understanding the concept of spherical wrist, divide the problem into inverse position and inverse orientation.
4. Compute first 3 joint angles using geometry.
5. Computer last 3 joint angles using known joint angles and transform matrices.
6. Use forward kinematics transform to calculate relative error for end-effector pose. 



[//]: # (Image References)

[image1]: ./resources/schematic.jpg
[image2]: ./resources/DHTable.jpg
[image3]: ./resources/theta1.jpg
[image4]: ./resources/theta2.jpg
[image5]: ./resources/theta3.jpg

### **Step - 1**: Using the kr210.urdf.xacro file, analyze Kuka KR210 and derive DH parameters. 

[urdfpath]:./kuka_arm/urdf/kr210.urdf.xacro

- Referring the listed joints section in [URDF file][urdfpath] to draw a schmatic of KR210 as follows:

![KR210 schematic with axis and other relevant DH parameters][image1]

- From the schematic, prepared DH Table shown below using the algorithm mentioned [here](https://classroom.udacity.com/nanodegrees/nd209/parts/7b2fd2d7-e181-401e-977a-6158c77bf816/modules/8855de3f-2897-46c3-a805-628b5ecf045b/lessons/87c52cd9-09ba-4414-bc30-24ae18277d24/concepts/2f59c902-9c32-4b26-9e52-5e495ec14dba).

![DH Table following the Craig JJ convention][image2]

### **Step - 2**: Define the transform matrices for each joint starting from base to gripper link using DH parameters.

- Following text shows the definition of transform matrices for each joint and homogenous transform matrix from base to gripper link using DH parameters in Python SymPy library.

```     
		# DH param symbols
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

        # Definition of individual transformation matrices
        
        T0_1 = Matrix([[             cos(q1),            -sin(q1),            0,              a0],
                       [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
                       [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),  cos(alpha0),  cos(alpha0)*d1],
                       [                   0,                   0,            0,               1]]

        T1_2 = Matrix([[             cos(q2),            -sin(q2),            0,              a1],
                       [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
                       [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),  cos(alpha1),  cos(alpha1)*d2],
                       [                   0,                   0,            0,               1]])
        

        T2_3 = Matrix([[             cos(q3),            -sin(q3),            0,              a2],
                       [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
                       [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),  cos(alpha2),  cos(alpha2)*d3],
                       [                   0,                   0,            0,               1]])

        T3_4 = Matrix([[             cos(q4),            -sin(q4),            0,              a3],
                       [ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
                       [ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3),  cos(alpha3),  cos(alpha3)*d4],
                       [                   0,                   0,            0,               1]])


        T4_5 = Matrix([[             cos(q5),            -sin(q5),            0,              a4],
                       [ sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
                       [ sin(q5)*sin(alpha4), cos(q5)*sin(alpha4),  cos(alpha4),  cos(alpha4)*d5],
                       [                   0,                   0,            0,               1]])


        T5_6 = Matrix([[             cos(q6),            -sin(q6),            0,              a5],
                       [ sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
                       [ sin(q6)*sin(alpha5), cos(q6)*sin(alpha5),  cos(alpha5),  cos(alpha5)*d6],
                       [                   0,                   0,            0,               1]])

        T6_G = Matrix([[             cos(q7),            -sin(q7),            0,              a6],
                       [ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
                       [ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6)*d7],
                       [                   0,                   0,            0,               1]])
        
        # Tranformation of base to gripper link is the resultant of multiplication of individual transformations
        T0_G = (T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G )
        
        # Given position and orientation, transformation matrix can be defined as follows:
        
        px = req.poses[x].position.x
        py = req.poses[x].position.y
        pz = req.poses[x].position.z
        
        # Orientation in euler angles from quaternion
        (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [req.poses[x].orientation.x, req.poses[x].orientation.y,
        req.poses[x].orientation.z, req.poses[x].orientation.w]) 
        
        # Position of end effector
        P_EE = Matrix([[px],[py],[pz]])       
        R_roll = Matrix([[ 1,              0,        0],
        [ 0,        cos(roll), -sin(roll)],
        [ 0,        sin(roll),  cos(roll)]])
        
        # Definition of rotation matrices from pitch, roll and yaw.
        R_pitch = Matrix([[ cos(pitch),        0,  sin(pitch)],
        [       0,        1,        0],
        [-sin(pitch),        0,  cos(pitch)]])

        R_yaw = Matrix([[ cos(yaw), -sin(yaw),        0],
        [ sin(yaw),  cos(yaw),        0],
        [ 0,              0,        1]])

        Rrpy = simplify(R_yaw * R_pitch * R_roll)
        
        # Using numpy, concatenating matrices
        
        T_final = np.vstack(np.hstack(Rrpy,P_EE),np.array([[0,0,0,1]])
        
        
```

### **Step - 3**: Understanding the concept of spherical wrist, divide the problem into inverse position and inverse orientation.

Simplifying the inverse kinematics problem for manipulator with spherical wrist, computing the center of the spherical wrist is easy since it is the location where the last three revolute joint Z-axis intersect. Hence from the end-effector position, wrist center can be computer as follows:

```
P_WC = simplify(P_EE - 0.303 * Rrpy * Matrix([[1],[0],[0]]))
```
Now, with wrist center known, first three joint angles can be computed using Geometry while last three by solving the equations derived from the transform matrices.
        
### **Step - 4**: Compute first 3 joint angles for inverse position problem using geometry.

Following are the figures explaining the geometry behind the equations to solve joint angles for first three joint.

Diagram representing the geometry for calculation for Theta 1 angle
![Diagram representing the geometry for calculation for Theta 1 angle][image3]



Diagram representing the geometry for calculation for Theta 2 angle
![Diagram representing the geometry for calculation for Theta 2 angle][image4]


Diagram representing the geometry for calculation for Theta 3 angle
![Diagram representing the geometry for calculation for Theta 3 angle][image5]

### **Step - 5**: Computer last 3 joint angles using known joint angles and transform matrices.
Following are the expression to find the last 3 angles.

```
R3_6 = T0_3.T[0:3,0:3] * Rrpy*(Rcorr).T

(theta4,theta5,theta6) = tf.transformations.euler_from_matrix(np.array(R3_6).astype(np.float64), "ryzy")
```

Where, 

R3_6 is the required rotation matrix, a function with 3 unknown joint angles.

T0_3 is the transform matrix derived using solved first three joint angles in Step 4.

Rrpy is the resultant rotation in world frame obtained from end-effector pose.

Rcorr is the calculated correction rotation to convert the Rrpy in world frame to frame for joint 3.

Transformations library API euler_from_matrix is used to obtain 3 angles, since they are in Y, Z and Y axis for joint 3 reference frame, 'ryzy' is used.


### **Step - 6**: Use forward kinematics transform to calculate relative error for end-effector pose.

Using forward kinematics and substituting the values of joint angles in transform, the pose of the end effector can be obtained. The expression looks as follows:

```
T0_G = (T0_1* T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_G ).evalf(subs={q1: theta1,q2: theta2,q3:theta3,q4: theta4,q5: theta5,q6:theta6})
(fk_roll, fk_pitch, fk_yaw) =tf.transformations.euler_from_matrix(np.array(T0_G[:3,:3]*R_correction))
(fk_x,fk_y,fk_z) = (T0_G[0,3],T0_G[1,3],T0_G[2,3])
```

Here, the R_correction is the rotation correction from reference frame of joint 6 to world reference frame.

The difference in roll, pitch, yaw and x,y,z positions coming from simulator and Forward kinematics gives us the error.
Averaging the errors for 100 poses, following are the values for corresponding variables:

delta_px = -1.46038017070196E-016	

delta_py = -3.37491182523096E-017	

delta_pz = 8.63506796930677E-017	

delta_roll = -1.5012210995102E-017	

delta_pitch = 7.27680907102691E-018	

delta_yaw = 5.70308851965852E-018

The error being in power of -16 or less shows the correct solution for inverse kinematics.
