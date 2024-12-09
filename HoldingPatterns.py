import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.io


def central_angle(center, point1, point2):
    # Calculate vectors
    v1 = (point1[0] - center[0], point1[1] - center[1])
    v2 = (point2[0] - center[0], point2[1] - center[1])
    
    # Compute the angle using atan2
    angle_radians = np.arctan2(v2[1]-v1[1], v2[0]-v1[0])# - math.atan2(v1[1], v1[0])
    
    # Normalize the angle to [0, 2*pi] if necessary
    angle_radians = (angle_radians + 2 * math.pi) % (2 * math.pi)
    
    # Convert to degrees (optional)
    angle_degrees = np.rad2deg(angle_radians)
    print('CA:', angle_radians, angle_degrees)
    
    return angle_radians, angle_degrees

OG_TOA = 5.3
ev_toa = 1.5*OG_TOA
hdg = 90
velocity = 57.412
turn_radius = (velocity*1.94384)**2/(11.26*math.tan(np.deg2rad(73)))
turn_radius*=0.3048

x1 = np.linspace(229, 229+2*turn_radius-.001, 100)
x2 = np.linspace(229+2*turn_radius-.001, 229, 100)
h = turn_radius+229
# print(h)
k = -391
# y =[]
y1 = [k+np.sqrt(turn_radius**2 - (i-h)**2) for i in x1]


y1[0] = y1[-1]
yStraight1 = []
yStraight1.append(y1[-1])
xStraight1 = []
xStraight1.append(x1[-1])
xStraight2 = []
xStraight2.append(x1[0])
toa = 0

while toa <= 5:
    yNew = yStraight1[-1] - 1
    dist = yNew-yStraight1[0]
    yStraight1.append(yNew)
    xStraight1.append(x1[-1])
    xStraight2.append(x1[0])
    toa = np.abs(dist/velocity)
    # print(toa)
    # yStraight1.
print(toa)
# central_angle1 = np.pi + np.arctan2(yStraight1[-1]-k, xStraight1[-1] - h)

exitLength = math.pi*turn_radius
exitTOA = exitLength/velocity
print(exitTOA)
# print(k- np.sqrt(turn_radius**2 - (0-h)**2))
y2 = [k- np.sqrt(turn_radius**2 - (i-h)**2) + dist for i in x1]
y2[0] = y2[-1]
# print(y2)

toa = (2*toa)+2*exitTOA
print(toa)

plt.plot(x1, y1)
plt.plot(xStraight1, yStraight1)
plt.plot(xStraight2, yStraight1)
plt.axis('equal')
plt.plot(x1, y2)
plt.show()

# holding_data = {
#     'turn_x': x1,
#     'turn_x2': x2,
#     'top_y': y1,
#     'bot_y': y2,
#     'right_x': xStraight1,
#     'left_x': xStraight2,
#     'straight_y': yStraight1
# }

# scipy.io.savemat('4XHold.mat', holding_data)