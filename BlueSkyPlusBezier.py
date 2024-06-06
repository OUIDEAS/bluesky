import bluesky as bs
from bluesky.traffic.trafficgroups import TrafficGroups
from bluesky.navdatabase import loadnavdata
import numpy as np
import time
import math
from bluesky.simulation import ScreenIO
import matplotlib.pyplot as plt
from bluesky.tools import geo, aero, areafilter, plotter
import re, sys, io
from pyproj import Proj
import MichaelBezierOptimizeLatLon as MBO

class ScreenDummy(ScreenIO):
    """
    Dummy class for the screen. Inherits from ScreenIO to make sure all the
    necessary methods are there. This class is there to reimplement the echo
    method so that console messages are printed.
    """
    def echo(self, text='', flags=0):
        """Just print echo messages"""
        print("BlueSky console:", text)



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

def rwgs84(latd):
    """ Calculate the earths radius with WGS'84 geoid definition
        In:  lat [deg] (latitude)
        Out: R   [m]   (earth radius) """
    lat    = np.radians(latd)
    a      = 6378137.0       # [m] Major semi-axis WGS-84
    b      = 6356752.314245  # [m] Minor semi-axis WGS-84
    coslat = np.cos(lat)
    sinlat = np.sin(lat)

    an     = a * a * coslat
    bn     = b * b * sinlat
    ad     = a * coslat
    bd     = b * sinlat

    # Calculate radius in meters
    r = np.sqrt((an * an + bn * bn) / (ad * ad + bd * bd))

    return r

def qdrdist(latd1, lond1, latd2, lond2):
    """ Calculate bearing and distance, using WGS'84
        In:
            latd1,lond1 en latd2, lond2 [deg] :positions 1 & 2
        Out:
            qdr [deg] = heading from 1 to 2
            d [nm]    = distance from 1 to 2 in nm """

    # Haversine with average radius for direction

    # Check for hemisphere crossing,
    # when simple average would not work

    # res1 for same hemisphere
    res1 = rwgs84(0.5 * (latd1 + latd2))

    # res2 :different hemisphere
    a    = 6378137.0       # [m] Major semi-axis WGS-84
    r1   = rwgs84(latd1)
    r2   = rwgs84(latd2)
    res2 = 0.5 * (abs(latd1) * (r1 + a) + abs(latd2) * (r2 + a)) / \
        (np.maximum(0.000001,abs(latd1) + abs(latd2)))

    # Condition
    sw   = (latd1 * latd2 >= 0.)

    r    = sw * res1 + (1 - sw) * res2

    # Convert to radians
    lat1 = np.radians(latd1)
    lon1 = np.radians(lond1)
    lat2 = np.radians(latd2)
    lon2 = np.radians(lond2)

    #root = sin1 * sin1 + coslat1 * coslat2 * sin2 * sin2
    #d    =  2.0 * r * np.arctan2(np.sqrt(root) , np.sqrt(1.0 - root))
    # d =2.*r*np.arcsin(np.sqrt(sin1*sin1 + coslat1*coslat2*sin2*sin2))

    # Corrected to avoid "nan" at westward direction
    d = r*np.arccos(np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1) + \
                 np.sin(lat1)*np.sin(lat2))

    # Bearing from Ref. http://www.movable-type.co.uk/scripts/latlong.html

    # sin1 = np.sin(0.5 * (lat2 - lat1))
    # sin2 = np.sin(0.5 * (lon2 - lon1))

    coslat1 = np.cos(lat1)
    coslat2 = np.cos(lat2)

    qdr = np.degrees(np.arctan2(np.sin(lon2 - lon1) * coslat2,
                                coslat1 * np.sin(lat2) -
                                np.sin(lat1) * coslat2 * np.cos(lon2 - lon1)))
    d
    return qdr, d #m    


def Meters_To_WSG84(waypoints, home):
        # convert position back to LAT/LONg
        # from pyproj import Proj
        p = Proj(proj='utm',zone=17, ellps='WGS84', preserve_units=False)
        homeX, homeY = p(home[1], home[0])
        waypoints = np.array(waypoints)
        asize = waypoints.shape
        if (len(asize) > 1):
            waypoints_LongLat = []
            for pt in waypoints:
                x = (pt[0] + homeX)
                y = (pt[1] + homeY)
                lon, lat = p(x,y,inverse=True)
                altitude = 150
                if(len(pt)>2):
                    altitude = pt[2]
                waypoints_LongLat.append([lat, lon, altitude])
            return waypoints_LongLat

        else:
            x = (waypoints[0] + homeX)
            y = (waypoints[1] + homeY)
            lon, lat = p(x,y,inverse=True)
            altitude = 0
            if(len(waypoints)>2):
                altitude = waypoints[2]
            return [lat, lon, altitude]

nm  = 1852.  # m       1 nautical mile\
tr = 111.6**2/(11.26*math.tan(np.deg2rad(25)))
print(tr)
#Start the Sim
start_time = time.time()
bs.init(mode ='sim', detached=True)
# bs.sim.ffmode = True
bs.scr = ScreenDummy()

#Pattern for pulling from terminal
pattern = r'Dist = (\d+\.\d+) nm'
#Creat Log
# bs.stack.stack(f'CRELOG MANUALPOSOrderedFlightAttempt 1')
# bs.stack.stack(f'SCEN MANUALPOSOrderedFlight2; SAVEIC MANUALPOSOrderedFlight2')
#Set Sim Time

#Read from datasets
wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
datasets_string = ['wptdata', 'aptdata']
datasets = [wptdata, aptdata]

#Make Aircraft

mytraf = bs.traf.cre('AC0', 'EC35', 39.42, -83.2, 0, 80, 111.6)
mytraf = bs.traf.cre('AC1', 'EC35', 39.4172, -83.2, 0, 80, 111.6)
mytraf = bs.traf.cre('AC2', 'EC35', 39.4144 , -83.2, 0, 80, 111.6)
mytraf = bs.traf.cre('AC3', 'EC35', 39.4116, -83.2, 0, 80, 111.6)
mytraf = bs.traf.cre('AC4', 'MV25', 39.4088 , -83.2, 0, 80, 111.6)

for acid in bs.traf.id:
    bs.stack.stack(f'BANK {acid} 45')

# bs.stack.stack(f'MANUALPOSOrderedFlightAttempt add traf.id, traf.lat, traf.lon, traf.alt, traf.tas, traf.hdg')
# bs.stack.stack(f'MANUALPOSOrderedFlightAttempt ON')



bs.stack.stack(f'DT .001')
#bs.stack.stack(f'DT 0.025') #Set DT to .025 second

# we'll run the simulation for up to 4000 seconds
# h, toi_dist = qdrdist(39.42, -83.2, 39.422055, -83.2)
# print('DIST: ', toi_dist)
# bs.stack.stack(f'DIST 39.42, -83.2, 39.415, -83.2')

t_max = 1000000
mytraf2 = bs.traf.cre('EM0', 'EC35', 39.4075, -83.2, 0, 80, 124.422)
ntraf = bs.traf.ntraf

n_steps = int(t_max + 1)
t = np.linspace(0, t_max, n_steps)


# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))
em = np.zeros((n_steps, 1))


count = np.zeros(len(bs.traf.id))

ETA = np.zeros((len(bs.traf.id), n_steps))
distances = np.zeros((len(bs.traf.id), n_steps))

TOI = np.zeros((len(bs.traf.id)-1, n_steps))
counter = [0, 0, 0, 0, 0]
wpts_xac = {}
wpts_yac = {}
#The sim
for i in range(n_steps):
    # print(i)
    # Perform one step of the simulation
    bs.sim.step()
    
    
    # save the results from the simulator in the results array,
    # here we keep the latitude, longitude, altitude and TAS
    res[i] = [bs.traf.lat,
                bs.traf.lon,
                bs.traf.alt,
                bs.traf.tas,
                bs.traf.hdg]
    if i!=0:

        
        for j in range(len(bs.traf.id)-1):
            
            h, toi_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], bs.traf.lat[-1], bs.traf.lon[-1])
            TOI[j][i] = toi_dist/((bs.traf.tas[-1] - bs.traf.tas[j])/1.94384)
            # print(TOI[j][i], bs.traf.id[j])
            if TOI[j][i] <= 10 and counter[j] == 0:
                counter[j] = 1
                if j%2 == 0:
                    lr = 1
                    koz_x = bs.traf.lon[j]+.000686
                else:
                    lr = 1
                    koz_x = bs.traf.lon[j]-.000686
                print(lr)
                target_toa1 = 1.5*3.71
                target_toa2 = 1.5*3.71
                #2200 ft = .00604 lat 750ft = .002055 lon 690ft = .0019 lon 450ft = .001234 lat 250ft = .000686 lon 50ft = .0001375 lat
                nodes1 = [np.array([bs.traf.lon[j], bs.traf.lon[j]+(lr*.0018), bs.traf.lon[j]+(lr*.0019)]).flatten(), 
                          np.array([bs.traf.lat[j]+.00604, bs.traf.lat[j]+.00604, bs.traf.lat[j]+.00604+.001234]).flatten()]
                nodes2 = [np.array([bs.traf.lon[j]+(lr*.0019), bs.traf.lon[j]+(lr*.0015), bs.traf.lon[j]]).flatten(), 
                          np.array([bs.traf.lat[j]+.00604+.001234, bs.traf.lat[j]+.00604+.0015, bs.traf.lat[j]+.00604+(2*.001234)]).flatten()]
                
                koz_bot = bs.traf.lat[j]+.00604+.0001375
                koz_top = bs.traf.lat[j]+.00604+(2*.001234)-.0001375
                print(koz_x, koz_bot, koz_top)
                wpts, wpts_x, wpts_y = MBO.toCallOutside(bs.traf.tas[j]/1.94384, np.deg2rad(20), target_toa1, target_toa2,
                                                          bs.traf.hdg[j], nodes1, nodes2, koz_x, koz_bot, koz_top, lr)
                # bs.stack.stack(f'DELWPT {bs.traf.id[j]} 40.4 -83.2')
                wpts.append([40.4, -83.2])
                wpts_xac[f'AC{j}'] = wpts_x
                wpts_yac[f'AC{j}'] = wpts_y

            # if toi_dist <= 46 and count[j] == 1:
            #     count[j] = 2
                for k in range(5, len(wpts_xac[bs.traf.id[j]])):
                    bs.stack.stack(f'ADDWPT {bs.traf.id[j]} {wpts[k][0]} {wpts[k][1]} 150')

    

end_time = time.time()
bs.scr.update()
print("SIM TIME", end_time-start_time)

'''Plotting'''


# print(ETA)
# print(res[4])
# print(res[5])
# markers = ['o', 's', '^', 'v', '>', '<', 'P', 'X', 'D', '*', '+', 'x', '|', '_', '1', '2', '3', '4', 'h', 'H']
used_colors = []
# for idx, acid in enumerate(bs.traf.id):
#     available_colors = [color for color in range(1, 101) if color not in used_colors]
#     color_idx = np.random.choice(available_colors)
#     color = plt.cm.tab20(color_idx)  # Use a colormap to generate a color
#     used_colors.append(color_idx)  # Add the used color index to the list
#     if idx%2 == 0:
#         marker = 'o'
#     elif idx%3 == 0:
#         marker = '^'
#     else:
#         marker = 'x'
#     # marker = np.random.choice()  # Randomly select a marker style
#     color = np.random.rand(3,)  # Randomly select a color
#     if acid == 'AC4' or acid == 'EM0':
#         print(idx)
#         plt.plot(res[::25, 1, idx], res[::25, 0, idx], marker = marker, label=f'{acid}', color=color

for idx, acid in enumerate(bs.traf.id):
    available_colors = [color for color in range(1, 101) if color not in used_colors]
    color_idx = np.random.choice(available_colors)
    color = plt.cm.tab20(color_idx)
    used_colors.append(color_idx)

    marker = 'o' if idx % 2 == 0 else '^' if idx % 3 == 0 else 'x'
    color = np.random.rand(3,)
    # if acid in wpts_xac and acid in wpts_yac:
    #     valid_wpts_x = [x for x in wpts_xac[acid] if x != 0]
    #     valid_wpts_y = [y for y in wpts_yac[acid] if y != 0]
    # plt.scatter(wpts_), valid_wpts_y, label=f'{acid} waypoints')

    # if acid == 'AC4' or acid == 'EM0':
    plt.plot(res[::25, 1, idx], res[::25, 0, idx], marker=marker, label=f'{acid}', color=color)

# print(wpts_xac.keys())
# plt.scatter(wpts_xac[bs.traf.id[4]], wpts_yac[bs.traf.id[4]])
plt.axis('equal')
# shape4PlotAirport('COJEZ', 'NIKOE')
# plotWPTCoords(wpts)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.legend()

plt.show()

# fig = plt.figure(2)
# ax = fig.add_subplot()
# for i in range(0, 4):
#     ax.scatter(wpts_xac[bs.traf.id[i]], wpts_yac[bs.traf.id[i]])
# ax.axis('equal')
# plt.show()
