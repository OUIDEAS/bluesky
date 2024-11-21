import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.io

OG_TOA = 5.3
ev_toa = 4*OG_TOA
hdg = 90
velocity = 57.412
turn_radius = (velocity*1.94384)**2/(11.26*math.tan(np.deg2rad(73)))
turn_radius*=0.3048

x1 = np.linspace(0, 2*turn_radius, 100)
x2 = np.linspace(2*turn_radius, 0, 100)
h = turn_radius
k = 564
# y =[]
y1 = [k+np.sqrt(turn_radius**2 - (i-h)**2) for i in x1]



yStraight1 = []
yStraight1.append(y1[-1])
xStraight1 = []
xStraight1.append(x1[-1])
xStraight2 = []
xStraight2.append(x1[0])
toa = 0

while toa <= ev_toa-3.88:
    yNew = yStraight1[-1] - 1
    dist = yNew-yStraight1[0]
    yStraight1.append(yNew)
    xStraight1.append(x1[-1])
    xStraight2.append(x1[0])
    toa = np.abs(dist/velocity)
    # print(toa)
    # yStraight1.
central_angle = np.pi + np.arctan2(yStraight1[-1]-k, xStraight1[-1] - h)
exitLength = 2*np.pi*turn_radius * (central_angle/(2*np.pi))
exitTOA = exitLength/velocity
print(exitTOA)
# print(k- np.sqrt(turn_radius**2 - (0-h)**2))
y2 = [k- np.sqrt(turn_radius**2 - (i-h)**2) + dist for i in x1]
# print(y2)


plt.plot(x1, y1)
plt.plot(xStraight1, yStraight1)
plt.plot(xStraight2, yStraight1)
plt.axis('equal')
plt.plot(x1, y2)
plt.show()

holding_data = {
    'turn_x': x1,
    'turn_x2': x2,
    'top_y': y1,
    'bot_y': y2,
    'right_x': xStraight1,
    'left_x': xStraight2,
    'straight_y': yStraight1
}

scipy.io.savemat('4XHold.mat', holding_data)