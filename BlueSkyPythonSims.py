import bluesky as bs
from bluesky import traffic
import numpy as np
import time
from bluesky.simulation import ScreenIO
import matplotlib.pyplot as plt
# from bluesky.network.client import Client
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit

class ScreenDummy(ScreenIO):
    """
    Dummy class for the screen. Inherits from ScreenIO to make sure all the
    necessary methods are there. This class is there to reimplement the echo
    method so that console messages are printed.
    """
    def echo(self, text='', flags=0):
        """Just print echo messages"""
        print("BlueSky console:", text)




#Start the Sim

bs.init(mode ='sim', detached=True)
bs.sim.ffmode = True
bs.scr = ScreenDummy()

n = 20

# bs.traf.mcre(n, actype ="A320")

bs.traf.mcre(n, 'B737')
# bs.traf.cre(acid = 'AC1', actype = 'B737', aclat  = 39.0, aclon = -1*82.0, achdg = 0, acalt= 1500, acspd= 300)

# bsclient.get_trajectories()

for acid in bs.traf.id:
    #Initial Path Waypoints
    bs.stack.stack(f'ORIG {acid} KCMH;'
                f'ADDWPT {acid} HEH FL150;'
                f'ADDWPT {acid} BUD FL100;'
                f'ADDWPT {acid} TVT FL200;'
                f'ADDWPT {acid} JPU FL150;'
                f'VNAV {acid} ON')

bs.stack.stack('DT 1;FF')

# we'll run the simulation for up to 4000 seconds
t_max = 1000

ntraf = bs.traf.ntraf
print(ntraf)
n_steps = int(t_max + 1)
t = np.linspace(0, t_max, n_steps)

# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))

# print(np.ndim(res))
# print(res[0][0][0])
# iteratively simulate the traffic
bs.traf
start_time = time.time()
for i in range(n_steps):
    # Perform one step of the simulation
    bs.sim.step()
    

    # save the results from the simulator in the results array,
    # here we keep the latitude, longitude, altitude and TAS
    res[i] = [bs.traf.lat,
                bs.traf.lon,
                bs.traf.alt,
                bs.traf.tas,
                bs.traf.hdg]
    # plt.plot(bs.traf.lon, bs.traf.lat, marker = '^')
    # plt.pause(0.0001)
end_time = time.time()
print("SIM TIME", end_time-start_time)

'''Plotting'''

# for idx, acid in enumerate(bs.traf.id):
#     fig = plt.figure(figsize=(10, 40))
#     ax1 = plt.subplot2grid((8, 1), (0, 0), rowspan=2)
#     ax2 = plt.subplot2grid((4, 1), (1, 0))
#     ax3 = plt.subplot2grid((4, 1), (2, 0))
#     ax4 = plt.subplot2grid((4, 1), (3, 0))
#     ax1.plot(res[:, 1, idx], res[:, 0, idx])
#     ax1.set_xlabel('lon')
#     ax1.set_ylabel('lat')
#     ax1.grid()

#     ax2.plot(t, res[:, 2, idx])
#     ax2.set_xlabel('t [s]')
#     ax2.set_ylabel('alt [m]')
#     ax2.grid()

#     ax3.plot(t, res[:, 3, idx])
#     ax3.set_xlabel('t [s]')
#     ax3.set_ylabel('TAS [m/s]')
#     ax3.grid()

#     ax4.plot(t, res[:, 4, idx])
#     ax4.set_xlabel('t [s]')
#     ax4.set_ylabel('heading')
#     ax4.grid()

#     fig.suptitle(f'Trajectory {acid}')


# markers = ['o', 's', '^', 'v', '>', '<', 'P', 'X', 'D', '*', '+', 'x', '|', '_', '1', '2', '3', '4', 'h', 'H']
used_colors = []
for idx, acid in enumerate(bs.traf.id):
    available_colors = [color for color in range(1, 101) if color not in used_colors]
    color_idx = np.random.choice(available_colors)
    color = plt.cm.tab20(color_idx)  # Use a colormap to generate a color
    used_colors.append(color_idx)  # Add the used color index to the list
    if idx%2 == 0:
        marker = 'o'
    elif idx%3 == 0:
        marker = '^'
    else:
        marker = 'x'
    # marker = np.random.choice()  # Randomly select a marker style
    color = np.random.rand(3,)  # Randomly select a color
    plt.plot(res[::10, 1, idx], res[::10, 0, idx], marker = marker, linestyle = '', label=f'AC{idx}', color=color)
# plt.plot(res[:, 1], res[:, 0])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.legend()

plt.show()
