import numpy as np
import matplotlib.pyplot as plt

# # corner_point = [1000, 0]
# # start_point = [750, -220]
# # angle = np.arctan2(corner_point[1] - start_point[1], corner_point[0]-start_point[0])
# # print(np.rad2deg(angle))
corridor_linex = np.zeros(900)
y = [i for i in range(0, 900)]
corridor_linex2 = [1500 for i in range(0, 900)]

evehicle_x1 = [500 for i in range(0, 900)]
evehicle_x2 = [1000 for i in range(0, 900)]


x = [i for i in range(750, 1500)]
y_sl = [i*np.tan(np.deg2rad(11.3))-200 for i in x]
# plt.plot(x, y_sl)


cpx = [750, 750, 1500, 1500]
cpy = [0, 900, 900, 0]
plt.grid()
plt.plot(corridor_linex, y, 'r')
plt.plot(corridor_linex2, y, 'r', linewidth = 2, label = 'UAM Flight Corridor Border')
plt.plot(evehicle_x1, y, 'b--')
plt.plot(evehicle_x2, y, 'b--')

plt.plot(750, 900, color = 'purple', marker = '*', markersize = 10, label = 'Vertiport')
# plt.scatter(cpx, cpy, marker = 'o', color = 'orange', s=50)
plt.fill_betweenx(y, evehicle_x1, corridor_linex, color='green', alpha = 0.5, label = 'Area for Control Point Placement')
plt.fill_betweenx(y, evehicle_x2, corridor_linex2, color='green', alpha = 0.5)
plt.fill_betweenx(y, evehicle_x2, evehicle_x1, color='blue', alpha = 0.5, label = 'Emergency Vehicle Clearance')
plt.legend(loc = 'lower left', fontsize = '7.5')
plt.axvline(x=min(corridor_linex), color='none')
plt.axvline(x=max(evehicle_x1), color='none')
plt.axhline(y=min(y), color='none')
plt.axhline(y=max(y), color='none')
plt.xlim([730, 1600])
plt.axis('equal')
plt.xlabel('X (ft)')
plt.ylabel('Y (ft)')
plt.show()


# ac = [i for i in range(0, 50)]
# ETA = [300 + 30*i for i in range(0, 50)]
# ETA_Real = []
# j = 0
# for i in ETA:
#     ETA_Real.append(i+10*j)
#     j+=1
# plt.scatter(ac, ETA, color = 'blue', label = 'Original ETA')
# plt.scatter(ac, ETA_Real, color = 'red', label = 'Actual ETA')
# plt.grid()
# plt.legend()
# plt.title('ETA of Aircraft in Fleet')
# plt.xlabel('Aircraft ID #')
# plt.ylabel('ETA (s)')
# plt.show()