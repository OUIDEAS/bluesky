import bluesky as bs
from bluesky.traffic.trafficgroups import TrafficGroups
from bluesky.navdatabase import loadnavdata
import numpy as np
import time
from bluesky.simulation import ScreenIO
import matplotlib.pyplot as plt
from bluesky.tools import geo, aero, areafilter, plotter
# from PyQt5.QtCore import Qt, QTimer
# from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTextEdit



class ScreenDummy(ScreenIO):
    """
    Dummy class for the screen. Inherits from ScreenIO to make sure all the
    necessary methods are there. This class is there to reimplement the echo
    method so that console messages are printed.
    """
    def echo(self, text='', flags=0):
        """Just print echo messages"""
        print("BlueSky console:", text)
    # def handle(self, cmdu):
    #     """
    #     Override the handle method to intercept and process commands.
    #     """
    #     if cmdu[:5] == 'ECHO ':
    #         self.echo(cmdu[5:])
    #     else:
    #         super().handle(cmdu)


def getWPLatLon(wpts):
    j = -1
    keys = []
    for i in wpts:
        for k in datasets_string:
            if j < 1:
                j+=1
            else:
                j=-0
            key = k[:2] + 'id'
            if np.isin(i, datasets[j][key]):
                keys.append(key)
                # print(i, key)
           
    wptlon = []
    wptlat = []
    for i in range(len(wpts)):
        if keys[i] == 'apid':
            # print(wpts[i])
            wptidx = bs.navdb.getaptidx(wpts[i])
            wptlon.append(aptdata['aplon'][wptidx])
            wptlat.append(aptdata['aplat'][wptidx])
        else:
            # print(wpts[i])
            wptidx = bs.navdb.getwpidx(wpts[i])
            wptlon.append(wptdata['wplon'][wptidx])
            wptlat.append(wptdata['wplat'][wptidx])

    return wptlat, wptlon
def shape4PlotAirport(ap1, ap2):
    #Check which dataset each wp is in
    j = -1
    for i in datasets_string:
        j+=1
        key = i[:2] + 'id'
        if np.isin(ap1, datasets[j][key]):
            ap1_type = key
            # print(key)
        if np.isin(ap2, datasets[j][key]):
            ap2_type = key
            # print(key)
    
    # Airport 1

    if ap1_type == 'apid':
        ap1idx = bs.navdb.getaptidx(ap1)
        ap1lat = aptdata['aplat'][ap1idx]
        ap1lon = aptdata['aplon'][ap1idx]
        print(aptdata['apid'][ap1idx])
    else:
        ap1idx = bs.navdb.getwpidx(ap1)
        ap1lat = wptdata['wplat'][ap1idx]
        ap1lon = wptdata['wplon'][ap1idx]

    

    # Airport 2
    if ap2_type == 'apid':
        ap2idx = bs.navdb.getaptidx(ap2)
        ap2lat = aptdata['aplat'][ap2idx]
        ap2lon = aptdata['aplon'][ap2idx]
    else:
        ap2idx = bs.navdb.getwpidx(ap2)
        # print(wptdata['wpid'][ap2idx])
        ap2lat = wptdata['wplat'][ap2idx]
        ap2lon = wptdata['wplon'][ap2idx]


    lons = [ap1lon, ap2lon, ap2lon, ap1lon, ap1lon]
    lats = [ap1lat, ap1lat, ap2lat, ap2lat, ap1lat]
    plt.plot(lons, lats)

def shape4PlotCoords(lon1, lat1, lon2, lat2, shape):
    if shape == 'square' or shape == 'SQUARE':
        lons = [lon1, lon2, lon2, lon1, lon1]
        lats = [lat1, lat1, lat2, lat2, lat1]
        plt.plot(lons, lats)

def plotWPTCoords(wpts):
     #Check which dataset each wp is in
    j = -1
    keys = []
    for i in wpts:
        for k in datasets_string:
            if j < 1:
                j+=1
            else:
                j=-0
            key = k[:2] + 'id'
            if np.isin(i, datasets[j][key]):
                keys.append(key)
                # print(i, key)
           
    wptlon = []
    wptlat = []
    for i in range(len(wpts)):
        if keys[i] == 'apid':
            # print(wpts[i])
            wptidx = bs.navdb.getaptidx(wpts[i])
            wptlon.append(aptdata['aplon'][wptidx])
            wptlat.append(aptdata['aplat'][wptidx])
        else:
            # print(wpts[i])
            wptidx = bs.navdb.getwpidx(wpts[i])
            wptlon.append(wptdata['wplon'][wptidx])
            wptlat.append(wptdata['wplat'][wptidx])
    # print(wptlon, wptlat)
    plt.scatter(wptlon, wptlat, marker = 's', label = 'Waypoint')
            

#Start the Sim
start_time = time.time()
bs.init(mode ='sim', detached=True)
bs.sim.ffmode = True
bs.scr = ScreenDummy()

#Creat Log
bs.stack.stack(f'CRELOG MYLOG 1')
bs.stack.stack(f'SCEN PYTHONSCEN; SAVEIC PYTHONSCEN')
#Set Sim Time

#Read from datasets
wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
datasets_string = ['wptdata', 'aptdata']
datasets = [wptdata, aptdata]

#Make Aircraft
n = 5

mytraf = bs.traf.mcre(n, 'B737', 150, 300 )
bs.stack.stack(f'MYLOG add traf.id, traf.lat, traf.lon, traf.alt, traf.tas, traf.hdg')
bs.stack.stack(f'MYLOG ON')
# '''Code to create a group of aircraft, kind of pointless'''
# # Get the list of aircraft identifiers
# aircraft_ids = [acid for acid in bs.traf.id]
# print(aircraft_ids)

# # Create a TrafficGroups object
# traf_group = TrafficGroups()

# # Retrieve the indices of the aircraft with the specified identifiers
# aircraft_indices = [bs.traf.id.index(ac_id) for ac_id in aircraft_ids if ac_id in bs.traf.id]

# # Add the aircraft to a group
# success, message = traf_group.group('Fleet 1', *aircraft_indices)
# if success:
#     print("Aircraft added to 'NewGroup' successfully.")
# else:
#     print("Failed to add aircraft to 'NewGroup':", message)

#Set Waupoints and Flight Level (Altitude)
wpts = ['NIKOE', 'LARGY', 'KCMH', 'COJEZ', 'JASLO']
fls = ['FL150', 'FL100', 'FL200', 'FL150' , 'FL100']
for acid in bs.traf.id:
    #Initial Path Waypoints 
    
    bs.stack.stack(
                #f'BANK {acid}, 60;'
                # f'ORIG {acid} KCMH;'
                # f'ADDWPTMODE {acid} TURNSPD;'
                # f'ADDWPTMODE {acid} FLYOVER;'
                # f'ADDWPT {acid} TURNSPD 5'
                f'ADDWPT {acid} {wpts[0]} {fls[0]};'
                f'{acid} AFTER {wpts[0]} ADDWPT {wpts[1]} {fls[1]};'
                f'{acid} AFTER {wpts[1]} ADDWPT {wpts[2]} {fls[2]};'
                f'{acid} AFTER {wpts[2]} ADDWPT {wpts[3]} {fls[3]};'
                f'{acid} AFTER {wpts[3]} ADDWPT {wpts[4]} {fls[4]};'
                f'{acid} RTA {wpts[4]} 300.0;'
                f'VNAV {acid} ON')
    #bs.stack.stack(f'LISTRTE {acid}') #List Waypoints
bs.stack.stack(f'DT 1; FF') #Set DT to 1 second
bs.stack.stack(f'ECHO TEST')
# we'll run the simulation for up to 4000 seconds
t_max = 500

ntraf = bs.traf.ntraf

n_steps = int(t_max + 1)
t = np.linspace(0, t_max, n_steps)


# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))
# print(res)
# print(n_steps)
# iteratively simulate the traffic
lats, lons = getWPLatLon(wpts)

count = np.zeros(len(bs.traf.id))
# print(bs.traf.id[0])
# print(count)

#The sim
for i in range(n_steps):
    # Perform one step of the simulation
    bs.sim.step()
    # if i >= 100 and count == 0:
    #     print("BREAK PRINT")

    #     # for acid in bs.traf.id:
    #     # print("REMOVING WAYPOINTS")
    #     # bs.stack.stack(f'DELWPT BUD')
    #     # print("WAYPOINT LIST WIHTOUT BUD")
    #     bs.stack.stack(f'LISTRTE {acid}')
    #     count = 1
    
    # save the results from the simulator in the results array,
    # here we keep the latitude, longitude, altitude and TAS
    res[i] = [bs.traf.lat,
                bs.traf.lon,
                bs.traf.alt,
                bs.traf.tas,
                bs.traf.hdg]
    for j in range(len(bs.traf.id)):

        if np.isclose(bs.traf.lat[j], lats[-1], atol = 0.05) and np.isclose(bs.traf.lon[j], lons[-1], atol = 0.05) and count[j] == 0:
            bs.stack.stack(f'ECHO {bs.traf.id[j]} {i}')
            count[j] = 1

end_time = time.time()
bs.scr.update()

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
    plt.plot(res[::25, 1, idx], res[::25, 0, idx], marker = marker, label=f'AC{idx}', color=color)

shape4PlotAirport('COJEZ', 'NIKOE')
plotWPTCoords(wpts)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.legend()

plt.show()
