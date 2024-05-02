import numpy as np
import matplotlib.pyplot as plt

corridor_linex = np.zeros(3200)
y = [i for i in range(-100, 3100)]
corridor_linex2 = [1500 for i in range(-100, 3100)]

evehicle_x1 = [750 for i in range(-100, 3100)]
evehicle_x2 = [1000 for i in range(-100, 3100)]
cpx = [750, 750, 1500, 1500]
cpy = [-50, 3050, 3050, -50]
plt.grid()
# plt.plot(corridor_linex, y, 'r')
plt.plot(corridor_linex2, y, 'r', linewidth = 2, label = 'UAM Flight Corridor Border')
# plt.plot(evehicle_x1, y, 'b--')
plt.plot(evehicle_x2, y, 'b--')

plt.plot(750, 3050, color = 'purple', marker = '*', markersize = 10, label = 'Vertiport')
plt.scatter(cpx, cpy, marker = 'o', color = 'orange', s=50)
# plt.fill_betweenx(y, evehicle_x1, corridor_linex, color='green', alpha = 0.5, label = 'Area for Control Point Placement')
plt.fill_betweenx(y, evehicle_x2, corridor_linex2, color='green', alpha = 0.5)
plt.fill_betweenx(y, evehicle_x2, evehicle_x1, color='blue', alpha = 0.5, label = 'Emergency Vehicle Clearance')
plt.legend(loc = 'upper center', fontsize = '7.5')
plt.axvline(x=min(corridor_linex), color='none')
plt.axvline(x=max(evehicle_x1), color='none')
plt.axhline(y=min(y), color='none')
plt.axhline(y=max(y), color='none')
plt.xlim([730, 1600])
# plt.axis('equal')
plt.show()
