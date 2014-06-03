#!/usr/bin/env python
"""
Delta Robot Kinematics code, liberally borrowed from tutorial:
http://forums.trossenrobotics.com/tutorials/introduction-129/delta-robot-kinematics-3276/
License: "You can freely use this code in your applications."

Note that the above tutorial defines the lengths f and e differently to me.  The
tutorial uses the side length of the equilateral triangle which has either the
servo output (for f) or the parallel link anchor (for e) in the middle of the side.

I use the convention that f or e is the displacement from the servo output/parallel
link anchor to the centre of the triangle, which is easier to measure.
"""

import math as maths
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


class SimulatedDeltaBot(object):

    def __init__(self, servo_link_length, parallel_link_length, servo_displacement, effector_displacement):
        self.e = effector_displacement
        self.f = servo_displacement
        self.re = parallel_link_length
        self.rf = servo_link_length

    def forward(self, theta1, theta2, theta3):
        """ 
        Takes three servo angles in degrees.  Zero is horizontal.
        return (x,y,z) if point valid, None if not 
        """
        t = self.f-self.e

        theta1, theta2, theta3 = maths.radians(theta1), maths.radians(theta2), maths.radians(theta3)

        # Calculate position of leg1's joint.  x1 is implicitly zero - along the axis
        y1 = -(t + self.rf*maths.cos(theta1))
        z1 = -self.rf*maths.sin(theta1)

        # Calculate leg2's joint position
        y2 = (t + self.rf*maths.cos(theta2))*maths.sin(maths.pi/6)
        x2 = y2*maths.tan(maths.pi/3)
        z2 = -self.rf*maths.sin(theta2)

        # Calculate leg3's joint position
        y3 = (t + self.rf*maths.cos(theta3))*maths.sin(maths.pi/6)
        x3 = -y3*maths.tan(maths.pi/3)
        z3 = -self.rf*maths.sin(theta3)

        # From the three positions in space, determine if there is a valid
        # location for the effector
        dnm = (y2-y1)*x3-(y3-y1)*x2
    
        w1 = y1*y1 + z1*z1
        w2 = x2*x2 + y2*y2 + z2*z2
        w3 = x3*x3 + y3*y3 + z3*z3

        # x = (a1*z + b1)/dnm
        a1 = (z2-z1)*(y3-y1)-(z3-z1)*(y2-y1)
        b1 = -((w2-w1)*(y3-y1)-(w3-w1)*(y2-y1))/2.0

        # y = (a2*z + b2)/dnm;
        a2 = -(z2-z1)*x3+(z3-z1)*x2
        b2 = ((w2-w1)*x3 - (w3-w1)*x2)/2.0

        # a*z^2 + b*z + c = 0
        a = a1*a1 + a2*a2 + dnm*dnm
        b = 2*(a1*b1 + a2*(b2-y1*dnm) - z1*dnm*dnm)
        c = (b2-y1*dnm)*(b2-y1*dnm) + b1*b1 + dnm*dnm*(z1*z1 - self.re*self.re)
 
        # discriminant
        d = b*b - 4.0*a*c
        if d < 0:
            return None # non-existing point

        z0 = -0.5*(b+maths.sqrt(d))/a
        x0 = (a1*z0 + b1)/dnm
        y0 = (a2*z0 + b2)/dnm
        return (x0,y0,z0)


if __name__ == '__main__':

    bot = SimulatedDeltaBot(servo_link_length = 85.0, parallel_link_length = 210.0,
                            servo_displacement = 72, effector_displacement = 20)

    step=10
    minServo = -10
    maxServo = 50

    points = []
    for t1 in range(minServo, maxServo, step):
        for t2 in range(minServo, maxServo, step):
            for t3 in range(minServo, maxServo, step):
                points.append(bot.forward(t1,t2,t3))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    surf = ax.scatter(xs=[x for x,y,z in points] ,ys=[y for x,y,z in points],zs=[z for x,y,z in points])
    plt.show()
