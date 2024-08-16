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
from matplotlib.patches import Circle
import Scenario2Outside as MBO
import scipy.io

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
        print(len(waypoints))
        # convert position back to LAT/LONg
        p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)  
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


def LongLat_To_WSG84_Meters(waypoints, home):
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    waypoints = np.array(waypoints)
    asize = waypoints.shape
    
    waypoints_meters = []
    if(len(asize) > 1): #m x 3 matrix       
        for pt in waypoints:
            if len(pt) > 2:
                newpt = __WSG84_To_Meters_Single(pt[0], home,projobject=p)
                if(asize[1]==3):#altitude
                    newpt=[newpt[0],newpt[1],pt[1], pt[2]]
            else:
                newpt = __WSG84_To_Meters_Single(pt, home, projobject=p)
                if(asize[1]==3):#altitude
                    newpt=[newpt[0],newpt[1],pt[1]]

            waypoints_meters.append(newpt)
    else:
        waypoints_meters=(__WSG84_To_Meters_Single([waypoints[0],waypoints[1]], home,projobject=p))
        altitude=0
        if(len(waypoints)==3):
            altitude = waypoints[2]
        waypoints_meters=[waypoints_meters[0],waypoints_meters[1], altitude]

    return waypoints_meters

def __WSG84_To_Meters_Single(waypoint, home, projobject):
    p = projobject
    homeX, homeY = p(home[1], home[0])  # Note the order: lon, lat
    lon = waypoint[1]
    lat = waypoint[0]
    x, y = p(lon, lat)  # Note the order: lon, lat
    x = (x - homeX)
    y = (y - homeY)
    return [x, y]

class gnss_converter():
    def __init__(self):
        self.R = 6378137.0
 
    def merc_y(self, lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * self.R
 
    def merc_x(self, lon):
        return math.radians(lon) * self.R

def get_utm_zone(longitude):
    return int((longitude + 180) / 6) + 1

def calc_Index(path, pos, idx, id):
    # print(idx)
    if idx == 0:
        dx = [pos[0] - x for x in path[0]]
        dy = [pos[1] - y for y in path[1]]
        idx = np.argmin(np.hypot(dx, dy))
        

    d_current = np.hypot(path[1][idx] - pos[1], path[0][idx] - pos[0])
    while True:
        if idx < len(path[0])-1:
            d_next = np.hypot(path[1][idx+1] - pos[1], path[0][idx+1] - pos[0])
            if d_current < d_next:
                break
        else:
            d_next = np.hypot(path[1][-1]-pos[1], path[0][-1]-pos[0])
            if d_current == d_next:
                break
        idx = idx+ 1 if (idx+1) < len(path[0]) else idx
        d_current = d_next
    
    while 10 > np.hypot(path[1][idx] - pos[1], path[0][idx] - pos[0]):
        if idx+1 >= len(path[0]):
            break
        idx+=1
    # print(idx, d_current)
    return idx, d_current

'''
Scenario 1 Options:
BezAM -> Uses Bezier alternate maneuver to avoid emergency aircraft
Hold -> Uses a holding pattern to avoid emergency aircraft

'''

scenario = 'Hold'
subscen = 'Single'

nm  = 1852.  # m       1 nautical mile\
# tr = 111.6**2/(11.26*math.tan(np.deg2rad(25)))
# print(tr)
#Start the Sim
start_time = time.time()
bs.init(mode ='sim', detached=True)
# bs.sim.ffmode = True
bs.scr = ScreenDummy()

#Pattern for pulling from terminal
pattern = r'Dist = (\d+\.\d+) nm'
#Create Log
# bs.stack.stack(f'CRELOG HoldFullTestScen1 1')
# bs.stack.stack(f'SCEN HoldFullTestScen1; SAVEIC HoldFullTestScen1')
bs.stack.stack(f'CRELOG Scenario2Bez 1')
bs.stack.stack(f'SCEN Scenario2Bez; SAVEIC Scenario2Bez')
#Set Sim Time

#Read from datasets
wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
datasets_string = ['wptdata', 'aptdata']
datasets = [wptdata, aptdata]

#Make Aircraft

mytraf = bs.traf.cre('AC0', 'M250', 39.5824, -82.2, 0, 80, 57.412)
mytraf = bs.traf.cre('AC1', 'M250', 39.5796, -82.2, 0, 80, 57.412)
mytraf = bs.traf.cre('AC2', 'M250', 39.5768, -82.2, 0, 80, 57.412)
mytraf = bs.traf.cre('AC3', 'M250', 39.574, -82.2, 0, 80, 57.412)
mytraf = bs.traf.cre('AC4', 'M250', 39.5712, -82.2, 0, 80, 57.412)

for acid in bs.traf.id:
    bs.stack.stack(f'BANK {acid} 30')
    bs.stack.stack(f'CONFLICTDETECTION {acid} OFF')
bs.stack.stack(f'ASAS OFF')

# bs.stack.stack(f'MANUALPOSOrderedFlightAttempt add traf.id, traf.lat, traf.lon, traf.alt, traf.tas, traf.hdg')
# bs.stack.stack(f'MANUALPOSOrderedFlightAttempt ON')



bs.stack.stack(f'DT .01')

#bs.stack.stack(f'DT 0.025') #Set DT to .025 second

# we'll run the simulation for up to 4000 seconds
# h, toi_dist = qdrdist(39.42, -83.2, 39.422055, -83.2)
# print('DIST: ', toi_dist)
# bs.stack.stack(f'DIST 39.42, -83.2, 39.415, -83.2')
utm_zone = get_utm_zone(-82.2)
p = Proj(proj='utm',zone=utm_zone,ellps='WGS84')
t_max = 7000
# t_max = 4000
# mytraf2 = bs.traf.cre('EM0', 'M250', 39.4075, -82.2, 0, 80, 64)
mytraf2 = bs.traf.cre('EM0', 'M250', 39.6, -82.225, 90, 80, 64)

ntraf = bs.traf.ntraf

n_steps = int(t_max + 1)
t = np.linspace(0, t_max, n_steps)


# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))
em = np.zeros((n_steps, 1))

homexy = [0, 0, 0, 0, 0]
homell = [0, 0, 0, 0, 0]
req_head = [[0,0], [0,0], [0,0], [0,0], [0,0]]
bas = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
guide = [0, 0, 0, 0, 0]
index = [0, 0, 0, 0, 0]

count = np.zeros(len(bs.traf.id))
violation_count = [0, 0, 0, 0, 0]
dist_to_EV = [[], [], [], [], []]


ETA = [[0, 0],[0, 0],[0, 0],[0, 0],[0, 0]]
distances = np.zeros((len(bs.traf.id), n_steps))
TOI = np.zeros((len(bs.traf.id)-1, n_steps))
counter = [0, 0, 0, 0, 0]

nodes = np.zeros((len(bs.traf.id)-1, 2), dtype = object)
waypts = [0, 0, 0, 0, 0]

# if scenario == 'BezAM':
#     start_end = [[0,0], [0,0], [0,0], [0,0], [0,0]]
#     times = [[0,0], [0,0], [0,0], [0,0], [0,0]]
# else:
start_end = [[[0,0],[0,0]], [[0,0],[0,0]], [[0,0],[0,0]], [[0,0],[0,0]], [[0,0],[0,0]]]
entry_point = [[0,0], [0,0], [0,0], [0,0],[0,0]]
print(entry_point[0][0])
times = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]]
# for i in range(0, 5):
#     print(times[i][2] )

#Handy Conversions
nmFeet = 6076.12
msFeet = 3.28084
nmMeters = 1852

lr = 1

sp_diff = 64-57.412
# print(nodes)
# nodes = np.array([[[0], [0]], [[0],[0]]], dtype = object)
# print(nodes[0][0])
# nodes[0][0][0]= [1, 2, 23, 4, 53]
# print(nodes[0][0])
#The sim

#home will be end point:
home_end = [39.599988478353495, -82.20158972364715]

ev_Goal = [39.6, -82.1994]

h, ev_dist = qdrdist(bs.traf.lat[-1], bs.traf.lon[-1], ev_Goal[0], ev_Goal[1])
print('EV GOAL STUFF', ev_dist, ev_dist/64)
ev_TOA = ev_dist/64

h, ac0_dist = qdrdist(bs.traf.lat[0], bs.traf.lon[0], home_end[0], home_end[1])
print('AC GOAL STUFF', ac0_dist, ac0_dist/57.412)

h, goal_dist = qdrdist(ev_Goal[0], ev_Goal[1], home_end[0], home_end[1])
print('GOAL DISTS', goal_dist)

'''
Initial ETAS
'''
for i in range(len(bs.traf.id)-1):
    goalHead, goalDist = qdrdist(bs.traf.lat[i], bs.traf.lon[i], home_end[0], home_end[1])
    print(goalDist)
    ETA[i][0] = goalDist/57.412
    print('INITIAL ETA:', ETA[i][0], bs.traf.id[i])
    
g0 = 0
pick = 0
evL = 0
k = 0
evT = 0
if scenario == 'BezAM':
    for i in range(n_steps):
        # print(i)
        
        # Perform one step of the simulation
        bs.sim.step()
        if np.isclose(39.6, bs.traf.lat[-1], atol = 0.00001) and np.isclose(-82.1994, bs.traf.lon[-1]):
            bs.traf.lat[-1] = 39.6
            bs.traf.lon[-1] = -82.1994
            bs.traf.tas[-1] = 0
            bs.traf.alt[-1] = 0
            evL = 1

        if g0 == 4:
            bs.traf.tas[pick] = 0
            bs.traf.alt[pick] = 0
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

                h, ev_dist = qdrdist(bs.traf.lat[-1], bs.traf.lon[-1], ev_Goal[0], ev_Goal[1])
                ev_TOA = ev_dist/64

                

                h, ac0_d2g = qdrdist(bs.traf.lat[0], bs.traf.lon[0], home_end[0], home_end[1])
                TOI[j][i] = (toi_dist)/(sp_diff)
                # print(ev_TOA)
                # print(TOI[j][i], bs.traf.id[j])
                # print('AIRCRAFT:', bs.traf.id[j], 'DIST TO EV:', toi_dist, 'TOI:', TOI[j][i])

                if j!=0 and bs.traf.lat[j] > bs.traf.lat[j-1] and guide[j] == 0:
                    # print('here!', bs.traf.id[j])
                    toi_dist = 999999
                    counter[j] = 1

                # if TOI[j][i] <= 10 and counter[j] == 0:
                if (j == 0 and ev_TOA <= 10 and counter[j] == 0):
                    # print('TOIDIST', toi_dist)
                    print(ac0_d2g)
                    counter[j] = 1
                    # if j%2 == 0 or j == 0:
                    #     lr = 1
                    koz_x = 55
                    # else:
                    #     lr = -1
                    #     koz_x = -7
                    # print(lr)
                    target_toa1 = 1.5*ev_TOA
                    target_toa2 = 1.5*ev_TOA

                    homexy[j] = [0, 0]
                    homell[j] = [bs.traf.lat[j], bs.traf.lon[j]] #save as lonlat for WSG84
                    # sjdlfh = Meters_To_WSG84([-129, 563], homell[j])
                    # print('AQGHDKASHDKJSA', sjdlfh)
                    # nodes1 = [np.array([homexy[j][0], 213, 221]).flatten(), np.array([ac0_d2g-137.5, ((ac0_d2g-137.5)+275)/60+(ac0_d2g-137.5), ac0_d2g])]
                    # nodes2 = [np.array([221, 236,homexy[j][0]]).flatten(), np.array([ac0_d2g, ac0_d2g+137.5, ac0_d2g])]

                    # koz_bot = 648+182.88
                    # koz_top = 60.96+648
                    # corridor = 228
                    nodes1 = [np.array([homexy[j][0], 213, 200]).flatten(), np.array([ac0_d2g-137.5-100, (ac0_d2g)/60+400, ac0_d2g+50])]
                    nodes2 = [np.array([200, 220, homexy[j][0]-129]).flatten(), np.array([ac0_d2g+50, ac0_d2g+400, ac0_d2g])]

                    koz_bot = ac0_d2g+45.5
                    koz_top = ac0_d2g-45.5
                    corridor = 228
                    wpts, wpts_x, wpts_y, b1, b2, nodes1, nodes2 = MBO.toCallOutside(bs.traf.tas[j], np.deg2rad(4.5), target_toa1, target_toa2,
                                                            bs.traf.hdg[j], nodes1, nodes2, koz_x, koz_bot, koz_top, lr, [bs.traf.lat[j], bs.traf.lon[j]], corridor, ev_TOA)
                    nodes[j][0] = nodes1
                    nodes[j][1] = nodes2
                elif j!=0 and np.isclose(start_end[0][0][0], bs.traf.lat[j], atol = 0.000001) and counter[j] == 0 and counter[0] != 0 :
                #     # print('TOIDIST', toi_dist)
                #     print(ac0_d2g)
                    counter[j] = 1
                    # if j%2 == 0 or j == 0:
                    #     lr = 1
                    koz_x = 55
                    # else:
                    #     lr = -1
                    #     koz_x = -7
                    # print(lr)
                    target_toa1 = 1.5*ev_TOA
                    target_toa2 = 1.5*ev_TOA

                    homexy[j] = [0, 0]
                    homell[j] = [bs.traf.lat[j], bs.traf.lon[j]] #save as lonlat for WSG84
                    # sjdlfh = Meters_To_WSG84([-129, 563], homell[j])
                    # print('AQGHDKASHDKJSA', sjdlfh)
                    # nodes1 = [np.array([homexy[j][0], 213, 221]).flatten(), np.array([ac0_d2g-137.5, ((ac0_d2g-137.5)+275)/60+(ac0_d2g-137.5), ac0_d2g])]
                    # nodes2 = [np.array([221, 236,homexy[j][0]]).flatten(), np.array([ac0_d2g, ac0_d2g+137.5, ac0_d2g])]

                    # koz_bot = 648+182.88
                    # koz_top = 60.96+648
                    # corridor = 228
                    nodes1 = [np.array([homexy[j][0], 213, 200]).flatten(), np.array([ac0_d2g-137.5-100, (ac0_d2g)/60+400, ac0_d2g+50])]
                    nodes2 = [np.array([200, 220, homexy[j][0]-129]).flatten(), np.array([ac0_d2g+50, ac0_d2g+400, ac0_d2g])]

                    koz_bot = ac0_d2g+45.5
                    koz_top = ac0_d2g-45.5
                    corridor = 228
                    wpts, wpts_x, wpts_y, b1, b2, nodes1, nodes2 = MBO.toCallOutside(bs.traf.tas[j], np.deg2rad(4.5), target_toa1, target_toa2,
                                                            bs.traf.hdg[j], nodes1, nodes2, koz_x, koz_bot, koz_top, lr, [bs.traf.lat[j], bs.traf.lon[j]], corridor, ev_TOA)
                    nodes[j][0] = nodes1
                    nodes[j][1] = nodes2

                if toi_dist <= 45.72:
                    violation_count[j]+=1
                    dist_to_EV[j].append(toi_dist)
                    # violation_count[j][1] = toi_dist
                    # print(toi_dist)
                
                if counter[j] ==1:
                    # print('TOIDIST', toi_dist)
                    counter[j] = 2
                    start_end[j][0][0], start_end[j][0][1]  = bs.traf.lat[j], bs.traf.lon[j]
                    times[j][0] = i/100
                    pos = [bs.traf.lat[j], bs.traf.lon[j]]
                    new = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
                    # print('POSITION IS:', new, 'POSITION IN LATLON IS:', pos, 'HOME IS:', homell[j])
                    entry, toa, wpts_x, wpts_y = MBO.EntryExitOutside(nodes1=nodes[0][0], nodes2=nodes[0][1], pos = new, velocity = bs.traf.tas[j], lr = lr, id = j)
                    end_point = [home_end[0], home_end[1]]
                    end_ll = Meters_To_WSG84(end_point, homell[j])
                    req_head[j][0] = entry[2]
                    bas[j][0] = entry[-1][0]
                    entry_point[j] = Meters_To_WSG84([entry[0][-1], entry[1][-1]], homell[j])
                    # print(entry_point[j])
                    waypts[j] = [wpts_x, wpts_y]
                    # if lr == -1:
                    #     toa = 21.57
                    # print('END POINT LL:', end_ll, 'END POINT XY:', end_point)
                    EVh, EV_dist = qdrdist(bs.traf.lat[-1], bs.traf.lon[-1], ev_Goal[0], ev_Goal[1])
                    print('EV TOA TO END:', EV_dist/64, EV_dist)
                    if EV_dist/64 < toa or evL == 1:
                        print('PATH IS VALID AND SHOULD ALLOW FOR PASSAGE', 'TIME STEP:', i)
                        # if lr == 1:
                        bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[-1][0]+1}; HDG {bs.traf.id[j]} {entry[2]}')
                        print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} AT {entry[-1][0]} DEGREES UNTIL A HEADING OF {entry[2]}')
                    elif j!=0:
                        print('PATH IS VALID AND SHOULD ALLOW FOR PASSAGE', 'TIME STEP:', i)
                        # if lr == 1:
                        bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[-1][0]+1}; HDG {bs.traf.id[j]} {entry[2]}')
                        print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} AT {entry[-1][0]} DEGREES UNTIL A HEADING OF {entry[2]}')
                    # # else:
                    #     bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[-1][0]}; HDG {bs.traf.id[j]} {90-entry[2]}')
                    #     print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} AT {entry[-1][0]} DEGREES UNTIL A HEADING OF {entry[2]}')
                    else:
                        print('GET FUCKED LOSER')
                
                # print(bs.traf.hdg[j])
                
                # if np.isclose(req_head[j][0], bs.traf.hdg[j], atol = 0.1) and counter[j] == 2:
                if guide[j] == 0 and np.isclose(entry_point[j][0], bs.traf.lat[j], atol = 0.0001) and np.isclose(entry_point[j][1], bs.traf.lon[j], atol = 0.0001):
                    guide[j] = 1
                    bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                    counter[j] = 3
                    print(f'INITIATING PURE PURSUIT GUIDANCE FOR {bs.traf.id[j]}')
                    print(ev_TOA)
                # elif lr == -1 and np.isclose(90-req_head[j][0], bs.traf.hdg[j]-360, atol = 0.1) and counter[j] == 2:
                #     guide[j] = 1
                #     bs.stack.stack(f'BANK {bs.traf.id[j]} 30')
                #     counter[j] = 3
                

                if guide[j] == 1 and counter[j] == 3:
                    pos = [bs.traf.lat[j], bs.traf.lon[j]]
                    m_pos = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
                    # print(waypts[j][1])
                    # print(m_pos[1])
                    # print(waypts[j])
                    index[j], d = calc_Index(waypts[j], m_pos, index[j], bs.traf.id[j])
                    
                    # if index[j] >= 2:
                    #     index[j] = -1
                    head = np.arctan2(waypts[j][1][index[j]] - m_pos[1], waypts[j][0][index[j]]-m_pos[0])

                    bs.stack.stack(f'HDG {bs.traf.id[j]} {90-np.rad2deg(head)}')
                    
                    # print(bs.traf.hdg[j], bs.traf.hdg[j]-360, 90-req_head[j][1], req_head[j][1], np.rad2deg(head))
                    # if index[j] == len(waypts[j][0]):
                    #     guide[j] = 2
                    
                if guide[j] == 1 and np.isclose(home_end[0],bs.traf.lat[j], atol = 1e-20) and np.isclose(home_end[1], bs.traf.lon[j], atol= 1e-20):
                    # bs.stack.stack(f'ALT {bs.traf.id[j]} 0 10000; SPD {bs.traf.id[j]} 0; HDG {bs.traf.id[j]} 0')
                    # bs.traf.tas[j] = 0
                    # bs.traf.alt[j] = 0
                    print(f'INITIATING LANDING FOR {bs.traf.id[j]}')
                    guide[j] = 4
                    g0 = guide[j]
                    pick = j
                    start_end[j][1][0], start_end[j][1][1]  = bs.traf.lat[j], bs.traf.lon[j]
                    times[j][1] = i/100
                # if guide[j] == 2:
                #     bs.traf.lat

                # if guide[j] == 3 and np.isclose(0, bs.traf.hdg[j], atol = 0.1):
                #     # goalHead, goalDist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], home_end[0], home_end[1])
                #     # print(goalDist)
                #     print(f'{bs.traf.id[j]} IS BACK ON THE NOMINAL PATH')
                #     start_end[j][1] = [bs.traf.lat[j], bs.traf.lon[j]]
                #     times[j][1] = i/100
                #     # ETA[j][1] = goalDist/57.412
                #     # print(f'ETA FOR {bs.traf.id[j]} BEFORE ALTERNATE MANEUVER IS {ETA[j][0]}. ETA AFTER ALTERNATE MANEUVER IS {ETA[j][1]}')
                #     guide[j] = 0
                

                    
        # time.sleep(0.025)

elif scenario == 'Hold':
    for i in range(n_steps):
        # print(i)
        # Perform one step of the simulation
        bs.sim.step()
        
        if np.isclose(39.6, bs.traf.lat[-1], atol = 0.00001) and np.isclose(-82.1994, bs.traf.lon[-1]):
            bs.traf.lat[-1] = 39.6
            bs.traf.lon[-1] = -82.1994
            bs.traf.tas[-1] = 0
            bs.traf.alt[-1] = 0
            evL = 1
        if g0 == 4:
            bs.traf.tas[pick] = 0
            bs.traf.alt[pick] = 0
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
                h, ev_dist = qdrdist(bs.traf.lat[-1], bs.traf.lon[-1], ev_Goal[0], ev_Goal[1])
                ev_TOA = ev_dist/64
                
                TOI[j][i] = (toi_dist)/(sp_diff)
                # print(TOI[j][i], bs.traf.id[j])
                # print('AIRCRAFT:', bs.traf.id[j], 'DIST TO EV:', toi_dist, 'TOI:', TOI[j][i])
                if j!=0 and bs.traf.lat[j] >= bs.traf.lat[j-1] and guide[j] == 0:
                    # print('here!', bs.traf.id[j])
                    toi_dist = 999999
                    counter[j] = 1

                # if TOI[j][i] <= 10 and counter[j] == 0:
                #     # print('TOIDIST', toi_dist)
                #     counter[j] = 1
                #     if j%2 == 0 or j == 0:
                #         lr = 1
                #     else:
                #         lr = -1
                #     print(lr)

                if toi_dist <= 45.72:
                    violation_count[j]+=1
                    dist_to_EV[j].append(toi_dist)
                    # violation_count[j][1] = toi_dist
                    # print(toi_dist)
                
                if (j == 0 and ev_TOA <= 10 and counter[j] == 0):
                    counter[j] =2
                    guide[j] = 1
                    # print(bs.traf.lat)
                    start_end[j][0] = [bs.traf.lat[j], bs.traf.lon[j]]
                    times[j][0] = i/100
                    bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                    if lr == 1:
                        bs.stack.stack(f'HDG {bs.traf.id[j]} 90; HDG {bs.traf.id[j]} 180')
                    else:
                        bs.stack.stack(f'HDG {bs.traf.id[j]} -90; HDG {bs.traf.id[j]} -180')
                    print(f'INITIATING HOLDING PATTERN FOR {bs.traf.id[j]} AT TIMESTAMP {i}')
                    k = i
                    evT = ev_TOA
                    print(ev_TOA)
                elif j!=0 and counter[j] == 1:
                    k = i
                    counter[j] =2
                    guide[j] = 1
                    # print(bs.traf.lat)
                    start_end[j][0] = [bs.traf.lat[j], bs.traf.lon[j]]
                    times[j][0] = i/100
                    bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                    if lr == 1:
                        bs.stack.stack(f'HDG {bs.traf.id[j]} 90; HDG {bs.traf.id[j]} 180')
                    else:
                        bs.stack.stack(f'HDG {bs.traf.id[j]} -90; HDG {bs.traf.id[j]} -180')
                    print(f'INITIATING HOLDING PATTERN FOR {bs.traf.id[j]} AT TIMESTAMP {i}')
                    print(ev_TOA)
                
                # if counter[j] == 2:
                #     print(i/100 - times[j][0], bs.traf.hdg[j], bs.traf.id[j])
                    # print(bs.traf.hdg[j])
                
                if np.isclose(180, bs.traf.hdg[j], atol = 0.01) and counter[j] == 2:
                    # print(bs.traf.lat)
                    times[j][1] = i/100
                    counter[j] = 3
                if subscen == 'Single':
                    if j!=0 and np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and counter[j] == 3 and bs.traf.lat[j]<bs.traf.lat[j-1]:
                        print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH AT TIMESTAMP {i}')
                        print(ev_TOA)
                        k = i
                        counter[j] = 4
                        if lr ==1:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} 270; HDG {bs.traf.id[j]} 0')
                        else:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} -270; HDG {bs.traf.id[j]} 0')
                        
                    elif j ==0 and np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and i/100 - times[j][1] >= 5 and counter[j] == 3:
                        print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH AT TIMESTAMP {i}')
                        print(ev_TOA)
                        k = i
                        counter[j] = 4
                        if lr ==1:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} 270; HDG {bs.traf.id[j]} 0')
                        else:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} -270; HDG {bs.traf.id[j]} 0')
                else:
                    if np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and i/100 - times[j][1] >= 5 and counter[j] == 3:
                        print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH')
                        print(ev_TOA)
                        counter[j] = 4
                        if lr ==1:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} 270; HDG {bs.traf.id[j]} 0')
                        else:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} -270; HDG {bs.traf.id[j]} 0')
                # if bs.traf.id[j] == 'AC4':
                #     print(counter[j], bs.traf.id[j])
                # if counter[j] == 4:
                    # print(start_end[j][0][0], bs.traf.lat[j])

                if np.isclose(start_end[j][0][0], bs.traf.lat[j], atol = 0.0001) and counter[j] == 4 and bs.traf.hdg[j] == 0:
                    start_end[j][1] = [bs.traf.lat[j], bs.traf.lon[j]]
                    times[j][2] = i/100
                    print(f'HOLDING PATTERN FOR {bs.traf.id[j]} ENDED AT TIMESTAMP {i}')
                    counter[j] = 5
                    guide[j] = 0
                    # guide[j] = 4
                    g0 = guide[j]
                    pick = j
                    start_end[j][1][0], start_end[j][1][1]  = bs.traf.lat[j], bs.traf.lon[j]
                    times[j][1] = i/100
                    pos_now = __WSG84_To_Meters_Single(start_end[j][1], start_end[j][0], p)
                    goal = __WSG84_To_Meters_Single(home_end, start_end[j][0], p)
                    head = np.rad2deg(np.arctan2(goal[1]-pos_now[1], goal[0]-pos_now[1]))
                    bs.stack.stack(f'HDG {bs.traf.id[j]} {0-head+90}')

                if guide[j] == 1 and np.isclose(home_end[0],bs.traf.lat[j], atol = 1e-20) and np.isclose(home_end[1], bs.traf.lon[j], atol= 1e-20):
                    # bs.stack.stack(f'ALT {bs.traf.id[j]} 0 10000; SPD {bs.traf.id[j]} 0; HDG {bs.traf.id[j]} 0')
                    # bs.traf.tas[j] = 0
                    # bs.traf.alt[j] = 0
                    print(f'INITIATING LANDING FOR {bs.traf.id[j]} at {i}')
                    guide[j] = 4
                    g0 = guide[j]
                    pick = j
                    # start_end[j][1][0], start_end[j][1][1]  = bs.traf.lat[j], bs.traf.lon[j]
                    # times[j][1] = i/100
                

        # time.sleep(0.0025)
print(ev_TOA)

end_time = time.time()
bs.scr.update()
print("SIM TIME", end_time-start_time)


'''
Computing ETA
'''
point_dist = [0, 0, 0, 0, 0]
nominal_time = [0, 0, 0, 0, 0]
real_time = [0, 0, 0, 0, 0]
delay = [0, 0, 0, 0, 0]



if scenario == 'BezAM':
    for i in range(0, 5):
        # print(i, counter[i], start_end[i][0][0])
        if counter[i] == 0:
            start_end[i][0] = [0, 0]
            start_end[i][1] = [0, 0]
        # print(i, counter[i])
        heh, point_dist[i] = qdrdist(start_end[i][0][0], start_end[i][0][1], start_end[i][1][0], start_end[i][1][1])
        nominal_time[i] = point_dist[i]/57.412
        real_time[i] = times[i][1]-times[i][0]
        delay[i] = real_time[i]-nominal_time[i]
        ETA[i][1] = ETA[i][0]+delay[i]
        print(f'DATA FOR {bs.traf.id[i]}')
        print(f'ALTERNATE MANEUVER START/END LATLON{start_end[i]}, DISTANCE BETWEEN POINTS {point_dist[i]}, STRAIGHT LINE TRAVEL TIME {nominal_time[i]}, ALTERNATE MANEUVER TRAVEL TIME {real_time[i]}, ABSORBED DELAY {delay[i]}')
else:
    for i in range(0,5):
        delay[i] = times[i][2] - times[i][0] 
        print(f'DATA FOR {bs.traf.id[i]}, {delay[i]}')
        ETA[i][1] = ETA[i][0]+delay[i]
        print(start_end[j])
        print(times[j])
print('ETA FOR EACH AIRCRAFT:', ETA)
# for i in range(len(violation_count)):
#     violation_count[i]/=100
# mins = []
# for i in dist_to_EV:
#     mins.append(np.argmin(i))

# print('TOTAL TIME EACH AIRCRAFT VIOLATED THE RADIUS:', violation_count)
# print('MINIMUM DISTANCE EACH AIRCRAFT GOT TO THE EV:', mins )
'''Plotting'''
# gpts = []
# for i in range(0, len(waypts[0][0])-1):
#     gpts.append([waypts[0][0][i], waypts[0][1][i]])
# print(gpts)
# goal = Meters_To_WSG84(gpts, homell[0])
# goalx = []
# goaly = []
# for i in range(0, len(goal)-1):
#     goalx.append(goal[i][1])
#     goaly.append(goal[i][0])
# print(goal[0])
# # print(ETA)
# print(res[4])
# print(res[5])
# markers = ['o', 's', '^', 'v', '>', '<', 'P', 'X', 'D', '*', '+', 'x', '|', '_', '1', '2', '3', '4', 'h', 'H']

# gpts = []
# j = 1
# for i in range(0, len(waypts[j][0])):
#     gpts.append([waypts[j][0][i], waypts[j][1][i]])
# # print(gpts)
# goal = Meters_To_WSG84(gpts, homell[j])
# goalx = []
# goaly = []
# for i in range(0, len(goal)):
#     goalx.append(goal[i][1])
#     goaly.append(goal[i][0])

# enpts = []
# for i in range(0, len(entry[0])):
#     enpts.append([entry[0][i], entry[1][i]])
# enp = Meters_To_WSG84(enpts, homell[j])
# enpx = []
# enpy = []
# for i in range(0, len(enp)):
#     enpx.append(enp[i][1])
#     enpy.append(enp[i][0])




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
# bez1_ll = Meters_To_WSG84(waypts[4], homell[4])
goalAC = plt.Circle((-82.20158972364715,39.599988478353495), 0.0005, color = 'g', fill = True, label = 'Landing Pad')
goalEV = plt.Circle((-82.1994, 39.6), 0.0005, color = 'r', fill = True, label = 'Emergency Vehicle Landing Pad')
fig = plt.figure(1)
ax = fig.add_subplot(111)
# for idx, acid in enumerate(bs.traf.id):
#     available_colors = [color for color in range(1, 101) if color not in used_colors]
#     color_idx = np.random.choice(available_colors)
#     color = plt.cm.tab20(color_idx)
#     used_colors.append(color_idx)
#     # plt.scatter(bez1_ll[1], bez1_ll[0])
#     marker = 'o' if idx % 2 == 0 else '^' if idx % 3 == 0 else 'x'
#     color = np.random.rand(3,)

#     # if acid == 'AC4' or acid == 'EM0':
#     ax.plot(res[::10, 1, idx], res[::10, 0, idx], marker=marker, label=f'{acid}', color=color, zorder = 100)
# goal = 

plt.plot(res[::10, 1, -1], res[::10, 0, -1], marker='x', label='Emergency Vehicle', color='red', zorder = 10)
plt.plot(res[::10, 1, 0], res[::10, 0, 0], marker='o', label='Aircraft 0', color='blue', zorder = 10)
plt.plot(res[::10, 1, 1], res[::10, 0, 1], marker='o', label='Aircraft 1', color='green', zorder = 10)
plt.plot(res[::10, 1, 2], res[::10, 0, 2], marker='o', label='Aircraft 2', color='orange', zorder = 10)
plt.plot(res[::10, 1, 3], res[::10, 0, 3], marker='o', label='Aircraft 3', color='purple', zorder = 10)
plt.plot(res[::10, 1, 4], res[::10, 0, 4], marker='o', label='Aircraft 4', color='brown', zorder = 10)

# ax.scatter(goalx, goaly, color = 'green', marker = '*', s = 30, zorder = 50)
ax.add_patch(goalEV)
ax.add_patch(goalAC)
ax.text(bs.traf.lon[-1]-0.005, bs.traf.lat[-1]+0.001, f'Emergency Vehicle TOA To Goal {ev_TOA:.3g}', fontsize = 7)
# # ax.text(bs.traf.lon[0]-0.0075, bs.traf.lat[0]-0.001, f'AC0 Begins Return To Nominal Path At t = {k/100}', fontsize = 7)
# ax.text(bs.traf.lon[1]-0.0075, bs.traf.lat[1]-0.001, f'AC1 Begins Return To Nominal Path At t = {k/100}, After AC0 Has Passed It', fontsize = 7)
# ax.text(bs.traf.lon[2]-0.0075, bs.traf.lat[2]+0.001, f'AC2 Begins Holding Pattern At t = 35.21, After It Has Passed AC1', fontsize = 7)

# # ax.scatter()\size = 7)
# # ax.scatter()\
# # ax.plot([bs.traf.lon[0]-0.001, bs.traf.lon[1]+0.001], [bs.traf.lat[0], bs.traf.lat[1]], color = 'blue', linestyle = 'dashed', linewidth = 2)
# ax.plot([bs.traf.lon[1]-0.005, bs.traf.lon[1]], [39.595, 39.595], color = 'black', linestyle = 'dashed', linewidth = 2)
# ax.plot([bs.traf.lon[1]-0.005, bs.traf.lon[1]], [bs.traf.lat[1], bs.traf.lat[1]], color = 'black', linestyle = 'dashed', linewidth = 2)
# ax.text(bs.traf.lon[1]-0.005, (bs.traf.lat[1]+39.595)/2, f'5 Second Travel Time', fontsize = 7)
# ax.text(bs.traf.lon[1]+0.001, bs.traf.lat[1]-0.001, f'AC1 Is Passing AC0', fontsize = 7)
plt.axis('equal')
# shape4PlotAirport('COJEZ', 'NIKOE')
# plotWPTCoords(wpts)
plt.title(f'Flight Snapshot at t = {t_max/100}')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.legend()

plt.show()

# if scenario == 'Hold':
#     data = {
#         'AC0': np.array([res[::10, 1, 0], res[::10, 0, 0], res[::10, 3, 0], res[::10, 4, 0]]),
#         'AC1': np.array([res[::10, 1, 1], res[::10, 0, 1], res[::10, 3, 1], res[::10, 4, 1]]),
#         'AC2': np.array([res[::10, 1, 2], res[::10, 0, 2], res[::10, 3, 2], res[::10, 4, 2]]),
#         'AC3': np.array([res[::10, 1, 3], res[::10, 0, 3], res[::10, 3, 3], res[::10, 4, 3]]),
#         'AC4': np.array([res[::10, 1, 4], res[::10, 0, 4], res[::10, 3, 4], res[::10, 4, 4]]),
#         'EM0': np.array([res[::10, 1, 5], res[::10, 0, 5], res[::10, 3, 5], res[::10, 4, 5]])
#     }
#     scipy.io.savemat('Scen2DataHold.mat', data)
# else:
#     data = {
#     'AC0': np.array([res[::10, 1, 0], res[::10, 0, 0], res[::10, 3, 0], res[::10, 4, 0]]),
#     'AC1': np.array([res[::10, 1, 1], res[::10, 0, 1], res[::10, 3, 1], res[::10, 4, 1]]),
#     'AC2': np.array([res[::10, 1, 2], res[::10, 0, 2], res[::10, 3, 2], res[::10, 4, 2]]),
#     'AC3': np.array([res[::10, 1, 3], res[::10, 0, 3], res[::10, 3, 3], res[::10, 4, 3]]),
#     'AC4': np.array([res[::10, 1, 4], res[::10, 0, 4], res[::10, 3, 4], res[::10, 4, 4]]),
#     'EM0': np.array([res[::10, 1, 5], res[::10, 0, 5], res[::10, 3, 5], res[::10, 4, 5]]),
#     'Waypts': np.array([goalx, goaly]),
#     'Entry': np.array([enpx, enpy])
#     }
#     scipy.io.savemat('Scen2Data.mat', data)
