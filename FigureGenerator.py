import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.patches import Circle
fig, ax = plt.subplots()
# # corner_point = [1000, 0]
# # start_point = [750, -220]
# # angle = np.arctan2(corner_point[1] - start_point[1], corner_point[0]-start_point[0])
# # print(np.rad2deg(angle))
# corridor_linex = np.zeros(900)
# y = [i for i in range(0, 900)]
# corridor_linex2 = [1500 for i in range(0, 900)]

# evehicle_x1 = [500 for i in range(0, 900)]
# evehicle_x2 = [1000 for i in range(0, 900)]


# x = [i for i in range(750, 1500)]
# y_sl = [i*np.tan(np.deg2rad(11.3))-200 for i in x]
# # plt.plot(x, y_sl)


# cpx = [750, 750, 1500, 1500]
# cpy = [0, 900, 900, 0]
# plt.grid()
# plt.plot(corridor_linex, y, 'r')
# plt.plot(corridor_linex2, y, 'r', linewidth = 2, label = 'UAM Flight Corridor Border')
# plt.plot(evehicle_x1, y, 'b--')
# plt.plot(evehicle_x2, y, 'b--')
# # plt.scatter(750, -100, marker = '^', color = 'green', s = 100, zorder = 20)
# plt.scatter(750, 900, color = 'purple', marker = '*', s = 100, label = 'Goal Point', zorder = 20)
# plt.scatter(750, 400, color = 'black', marker = '^', s = 100, zorder = 20, label = 'Emergency Aircraft') 
# # plt.scatter(cpx, cpy, marker = 'o', color = 'orange', s=50)
# plt.fill_betweenx(y, evehicle_x1, corridor_linex, color='green', alpha = 0.5, label = 'Area for Control Point Placement')
# plt.fill_betweenx(y, evehicle_x2, corridor_linex2, color='green', alpha = 0.5)
# plt.fill_betweenx(y, evehicle_x2, evehicle_x1, color='blue', alpha = 0.5, label = 'Emergency Vehicle Clearance')

# plt.axvline(x=min(corridor_linex), color='none')
# plt.axvline(x=max(evehicle_x1), color='none')
# plt.axhline(y=min(y), color='none')
# plt.axhline(y=max(y), color='none')
# plt.xlim([730, 1600])
# ax.add_patch(Circle((750, 400), 150, color='red', fill=True, label = 'Safety Radius'))
# plt.legend(loc = 'lower left', fontsize = '7.5')
# plt.axis('equal')
# plt.xlabel('X (ft)')
# plt.ylabel('Y (ft)')
# plt.show()


# ac = [i for i in range(0, 6)]
# ETA = [300 + 30*i for i in range(0, 6)]
# ETA_Real = []
# ETA_Other = []
# j = 0
# for i in ETA:
#     ETA_Real.append(i+5*j)
#     ETA_Other.append((i**2)/300)
#     j+=1
# ac_evens = [i for i in ac if i%2 == 0]
# plt.scatter(ac, ETA, color = 'blue', label = 'Original ETA')
# plt.scatter(ac, ETA_Real, color = 'red', label = 'ETA with Bezier Curve', marker = '^')
# plt.scatter(ac, ETA_Other, label = 'ETA with Traditional Avoidance', marker = '*', color = 'green')
# plt.grid()
# plt.legend()
# plt.xticks(ac_evens)
# # plt.title('ETA of Aircraft in Fleet')
# plt.xlabel('Aircraft ID #')
# plt.ylabel('ETA (s)')
# plt.show()

# times = ['10:00', '11:00', '12:30', '10:30', '12:00', '2:00']
# # plt.grid()
# plt.scatter('10:00', 'Original', marker = '*', color = 'black', s = 100, label = 'Aircraft A', zorder = 20)
# plt.scatter('10:30', 'Delayed', marker = '*', color = 'black', s = 100, zorder = 20)

# plt.scatter('11:00', 'Original', marker = '*', color = 'red', s = 100, label = 'Aricraft B', zorder = 20)
# plt.scatter('11:20', 'Delayed', marker = '*', color = 'red', s = 100, zorder = 20)

# plt.scatter('12:00', 'Original', marker = '*', color = 'green', s = 100, label = 'Aircraft C', zorder = 20)
# plt.scatter('12:10', 'Delayed', marker = '*', color = 'green', s = 100, zorder = 20)
# plt.legend()
# plt.plot(['10:00', '10:30'], ['Delayed', 'Delayed'],'k--', zorder = 20)
# plt.plot(['10:00', '10:00'], ['Original', 'Delayed'],'k', zorder = 20)

# plt.plot(['11:00', '11:20'], ['Delayed', 'Delayed'],'r--', zorder = 20)
# plt.plot(['11:00', '11:00'], ['Original', 'Delayed'],'r', zorder = 20)

# plt.plot(['12:00', '12:10'], ['Delayed', 'Delayed'],'g--', zorder = 20)
# plt.plot(['12:00', '12:00'], ['Original', 'Delayed'],'g', zorder = 20)
# # ax.text('1:00', 'Delayed', r'60 Minute Delay', va = 'baseline')
# ax.annotate('30 Minute Delay', xy=('10:30', 'Delayed'), xytext=(35, 10), textcoords='offset points', ha='center', color='black')
# ax.annotate('20 Minute Delay', xy=('11:20', 'Delayed'), xytext=(35, 10), textcoords='offset points', ha='center', color='red')
# ax.annotate('10 Minute Delay', xy=('12:10', 'Delayed'), xytext=(35, 10), textcoords='offset points', ha='center', color='green')
# # ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
# ax.invert_xaxis()
# # plt.scatter(times, y_points)
# plt.axis('equal')
# plt.xlabel('ETA')
# plt.ylabel('Status')
# plt.show()

# evt = 900/242
# b1a = evt*1.5
# print((2*b1a)/(evt))

# tr = 111.6**2/(11.26*np.tan(np.deg2rad(30)))
# print('TR', tr)
# h, k  = tr+750, -1282

# x_entry = [i for i in np.linspace(750, 2*tr+750, 300)]

# y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
# y_entry[-1] = -1282
# y_entry[0] = -1282
# print(y_entry[0])
# print(len(y_entry))
# d = 188*5

# y_w = [i for i in np.linspace(y_entry[-1], y_entry[-1]-d, 2)]
# x_sr = [x_entry[-1] for i in y_w]
# x_sl = [x_entry[0] for i in y_w]

# h2, k2  = tr+750, -1282-d
# y_exit = [k2-np.sqrt(tr**2 - (x-h2)**2) for x in x_entry]
# y_exit[-1] = -1282-d
# y_exit[0] = -1282-d

# plt.plot(x_entry, y_entry, color = 'blue', label = 'Maximum Bank Angle Turn')
# # plt.scatter(h, k)
# # plt.plot([750, h, h+tr], [-1282, k, -1282])
# plt.grid()
# plt.plot(x_sr, y_w, color = 'orange', label = 'Straight Segment')
# plt.plot(x_sl, y_w, color = 'orange')
# plt.plot(x_entry, y_exit, color = 'blue')
# # plt.scatter(h2, k2)
# plt.axis('equal')
# plt.legend()
# plt.show()
x_labels = ['AC0', 'AC1', 'AC2', 'AC3', 'AC4']

ETA = [[348.5412648256665, 350.4037558669337, 364.8412648256665], 
       [353.9630462925526, 353.9630462925526, 370.3830462925526], 
       [359.3848286351499, 359.3848286351499, 375.92482863514994], 
       [364.806611843445, 364.806611843445, 381.46661184344504], 
       [370.2283959266545, 370.2283959266545, 387.0083959266545]]
pc = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]

for i in range(0, 5):
    pc[i][0] = 100*((ETA[i][1]-ETA[i][0])/ETA[i][0])
    pc[i][1] = 100*((ETA[i][2]-ETA[i][0])/ETA[i][0])



# Separate the data for each column
eta_0 = [row[0] for row in ETA]
eta_1 = [row[1] for row in ETA]
eta_2 = [row[2] for row in ETA]

# Bar width
bar_width = 0.35

# X positions for the bars
x = np.arange(len(x_labels))


bars0 = ax.bar(x-bar_width/3, eta_0, width=bar_width, color='green', label='Nominal ETA')
bars1 = ax.bar(x+1*bar_width/3, eta_1, width=bar_width, color='blue', label='Bezier Alternate Maneuver ETA')
bars2 = ax.bar(x+bar_width, eta_2, width=bar_width, color='red', label='Holding Pattern Alternate Maneuver ETA')

# Adding the ETA[row][0] values as text on both bars
# for i in range(len(x_labels)):
#     ax.text(x[i] - bar_width/2, eta_1[i] + 0.1, f'{eta_0[i]:.1f}', ha='center', color='white', fontweight='bold')
#     ax.text(x[i] + bar_width/2, eta_2[i] + 0.1, f'{eta_0[i]:.1f}', ha='center', color='white', fontweight='bold')

# Set x-ticks and labels
ax.set_xticks(x)
ax.set_xticklabels(x_labels)


# Labeling
ax.set_xlabel('AC')
ax.set_ylabel('ETA Values')
ax.set_title('Bar Chart with Different Colored Bars')
ax.legend()

# Display the plot
plt.show()
