import bluesky as bs
from bluesky.traffic.trafficgroups import TrafficGroups
from bluesky.navdatabase import loadnavdata
import numpy as np
import time
import math
from math import pi as pi
from bluesky.simulation import ScreenIO
import matplotlib.pyplot as plt
from bluesky.tools import geo, aero, areafilter, plotter
import re, sys, io
from pyproj import Proj
import AircraftSpacingOptimize as MBO
import scipy.io
import pandas as pd
import json
import argparse
import scipy
import os
import SequencingOutside as SQO

class ScreenDummy(ScreenIO):
    """
    Dummy class for the screen. Inherits from ScreenIO to make sure all the
    necessary methods are there. This class is there to reimplement the echo
    method so that console messages are printed.
    """
    def echo(self, text='', flags=0):
        """Just print echo messages"""
        print("BlueSky console:", text)

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
    return qdr, d #m    

def Meters_To_WSG84(waypoints, home):
        # convert position back to LAT/LONg
        p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)  
        homeX, homeY = p(home[1], home[0])
        waypoints = np.array(waypoints)
        asize = waypoints.shape
        lats = []
        lons = []
        if (len(asize) > 1):
            waypoints_LongLat = []
            for pt in range(0, len(waypoints[0])):
                x = (waypoints[0][pt] + homeX)
                y = (waypoints[1][pt] + homeY)
                # print(x,y)
                lon, lat = p(x,y,inverse=True)
                altitude = 150
                if(len(waypoints)>2):
                    altitude = waypoints[2][pt]
                waypoints_LongLat.append([lat, lon, altitude])
                lats.append(lat)
                lons.append(lon)
            return waypoints_LongLat, lats, lons

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
            idx = idx
            break
        idx+=1
    # print(idx, d_current)
    return idx, d_current

def calc_Index_ll(path, pos, idx, id):
    # print(idx)
    if idx == 0:
        dy = [pos[0] - x for x in path[0]]
        dx = [pos[1] - y for y in path[1]]
        idx = np.argmin(np.hypot(dx, dy))
        

    d_current = np.hypot(path[0][idx] - pos[0], path[1][idx] - pos[1])
    while True:
        if idx < len(path[0])-1:
            d_next = np.hypot(path[0][idx+1] - pos[0], path[1][idx+1] - pos[1])
            if d_current < d_next:
                break
        else:
            d_next = np.hypot(path[0][-1]-pos[0], path[1][-1]-pos[1])
            if d_current == d_next:
                break
        idx = idx+ 1 if (idx+1) < len(path[0]) else idx
        d_current = d_next
    
    while 10 > np.hypot(path[0][idx] - pos[0], path[1][idx] - pos[1]):
        if idx+1 >= len(path[0]):
            idx = idx
            break
        idx+=1
    # print(idx, d_current)
    return idx, d_current

def convert_numpy(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.generic):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(i) for i in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_numpy(i) for i in obj)
    else:
        return obj

def to_json(obj, name, path):
    expanded_path = os.path.expanduser(path)

    os.makedirs(expanded_path, exist_ok=True)

    data_serializable = convert_numpy(obj)

    output_file = os.path.join(expanded_path, f"{name}.json")

    with open(output_file, "w") as f:
        json.dump(data_serializable, f, indent=4)

    print(f"{name} data saved to {output_file}")

def tracon_points(r, n):

    x = [math.cos(2*pi/n*x)*r for x in range(0,n+1)]
    y = [math.sin(2*pi/n*x)*r for x in range(0,n+1)]

    return x, y

def get_ll(x, y, home_ll, picked):
    i = np.random.randint(0, 200)
    j = np.random.randint(500, 3000)
    c = 0  
    
    while c == 0:
        if i in picked:
            i = np.random.randint(0, 200)
        else:
            picked.append(i)
            c=1

    xp = x[i]
    yp = y[i]

    if xp>0 and yp>0:
        xp+=j
        yp+=j
    elif xp>0 and yp<0:
        xp+=j
        yp-=j
    elif xp<0 and yp>0:
        xp-=j
        yp+= j
    elif xp<0 and yp<0: 
        xp-=j
        yp-=j

    hdg = np.rad2deg(np.arctan2(yp, xp))

    correct = 90-hdg
    diff = correct-hdg
    hdg+=diff
    if i <= 99:
        hdg+=180
    else:
        hdg-=180
   

    lat, lon, alt = Meters_To_WSG84([xp, yp], home_ll)
    return lat, lon, alt, hdg, picked

def find_quadrant(picked):
    if picked<=50:
        quad = 'I'
    elif 50<picked<=100:
        quad = 'II'
    elif 100<picked<=151:
        quad = 'III'
    else:
        quad = 'IV'
    return quad, quad

def get_tran(gate, lat, lon, t_tran, fixes, v):
    # fixes -> R, 1,2,2,1
     # fixes_l = [favus_path, teeze_path, elupy_path, xavyr_path]
    for i in range(len(lat)):
        if gate[i] == 'I' or gate[i] == 'II':
            h, d_r1 = qdrdist(lat[i], lon[i], fixes[0][0][0], fixes[0][0][1])
            h2, d_r2 = qdrdist(lat[i], lon[i], fixes[1][0][0], fixes[1][0][1])
            t_r1 = d_r1/v[i]
            t_r2 = d_r2/v[i]
            t_tran[i] = [t_r1, t_r2]

        elif gate[i][0] == 'IV' or gate[i] == 'III':
            h2, d_r1 = qdrdist(lat[i], lon[i], fixes[3][0][0], fixes[3][0][1])
            h, d_r2 = qdrdist(lat[i], lon[i], fixes[2][0][0], fixes[2][0][1])
            t_r1 = d_r1/v[i]
            t_r2 = d_r2/v[i]
            t_tran[i] = [t_r1, t_r2]

    # print(t_tran)
    return t_tran

def makeWPTS(wpts):
    lats = [i[0] for i in wpts]
    lons = [i[1] for i in wpts]
    return [lats, lons]

def findWPs(acid, run, fixes_l):
    index = [i for i, entry in enumerate(run) if entry[1] == acid]
    print(acid, index[0])
    # fixes_l = [favus_path, teeze_path, elupy_path, xavyr_path]
    # TR/G1R1, TL/G2R2, BL/G3R2, BR/G4R1
    seq_dat = run[index[0]]
    if seq_dat[3] == 'II':
        if seq_dat[4] == 'R2':
            print('TEEZE')
            wp_list = makeWPTS(fixes_l[1])
        else:
            print('FAVUS')
            wp_list = makeWPTS(fixes_l[0])
    if seq_dat[3] == 'I':
        if seq_dat[4] == 'R1':
            print('FAVUS')
            wp_list = makeWPTS(fixes_l[0])
        else:
            print('TEEZE')
            wp_list = makeWPTS(fixes_l[1])
    if seq_dat[3] == 'III':
        if seq_dat[4] == 'R2':
            print('ELUPY')
            wp_list = makeWPTS(fixes_l[2])
        else:
            print('XAVYR')
            wp_list = makeWPTS(fixes_l[3])
    if seq_dat[3] == 'IV':
        if seq_dat[4] == 'R1':
            print('XAVYR')
            wp_list = makeWPTS(fixes_l[3])
        else:
            print('ELUPY')
            wp_list = makeWPTS(fixes_l[2])
    # print(wp_list)
    return index, wp_list

def WPGuidance(wp_list, idx, pos, d, gate, i, id):
    # print(wp_list[idx])
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    kcmh = [39.997472, -82.891194]
    m_pos = __WSG84_To_Meters_Single(pos, kcmh, p)
    wp_pos = __WSG84_To_Meters_Single([wp_list[0][idx], wp_list[1][idx]], kcmh, p)
    head = np.rad2deg(np.arctan2(wp_pos[1]-m_pos[1], wp_pos[0]-m_pos[0]))
    h, dn = qdrdist(pos[0], pos[1], wp_list[0][idx], wp_list[1][idx])
    # if i%100 == 0:
    #     print(f'{id} at timestep {i} with a new distance of {dn} and old of {d}')
    # print(bs.traf.id[j], h, head, idx)
    # correct = 90-head
    # diff = correct-head
    # head+=diff
    # if gate == 'I' or gate == 'II':
    #     head+=180
    # else:
    #     head-=180
    # print(head)
    if dn > d and dn <= 50:
        idx+=1
        if idx >= len(wp_list):
            return head, dn, idx-1, 'done'
    return head, dn, idx, 'not_done'
    
    



start_time = time.time()
bs.init(mode ='sim', detached=True)
bs.scr = ScreenDummy()

#Pattern for pulling from terminal
pattern = r'Dist = (\d+\.\d+) nm'
#Create Log
# bs.stack.stack(f'CRELOG HoldFullTestScen1 1')
# bs.stack.stack(f'SCEN HoldFullTestScen1; SAVEIC HoldFullTestScen1')
# bs.stack.stack(f'CRELOG HoldSingleBypass 1')
# bs.stack.stack(f'SCEN HoldSingleBypass; SAVEIC HoldSingleBypass')
#Set Sim Time

#Read from datasets
wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
datasets_string = ['wptdata', 'aptdata']
datasets = [wptdata, aptdata]

# bs.traf.mcre(10, 'M250')
#R2 is 28 L
nm = 1852
points = 200
tracon = 10*nm
center = 30*nm # OPTIMIZATION HORIZON
freeze = 20.7*nm
other = 36*nm

'''
KCMH Airport and Meter Fixes for different STAR Routes
All fixes lead to either Runway 28L/R, aka R1 and R2 when sorted
'''

kcmh = [39.997472, -82.891194]
favus = [40.058917, -82.639417]
teeze = [40.0865, -82.848111]
elupy = [39.831083, -83.074778]
xavyr = [39.906167, -82.647333]


# TEEZE WAYPOINTS
melzz = [40.313417, -83.247194]
dubln = [40.202944, -83.132306]
trlgy = [40.167944, -83.061139]
polrs = [40.111194, -82.953444]
taces = [40.090111, -82.916333]

teeze_path = [melzz, dubln, trlgy, polrs, taces, teeze]

# FAVUS WAYPOINTS
bugzz = [40.565, -82.454056]
cbuss = [40.325306, -82.5405]
molls = [40.132139, -82.609194]

favus_path =[cbuss, molls, favus]

# ELUPY WAYPOINTS
jaktz = [39.591028, -83.419583]
rscot = [39.722389,  -83.286306]
obetz = [39.787667, -83.158389]
edwib = [39.877472, -82.984861]
gagbe = [39.907167, -82.927278]
jesce = [39.903556, -82.858889]

elupy_path = [rscot, obetz, elupy, edwib, gagbe, jesce]

# XAVYR WAYPOINTS
scrlt = [39.502917, -82.350833]
brtus = [39.730944, -82.473083]

xavyr_path = [brtus, xavyr]

fixes = [favus, dubln, elupy, xavyr]

fixes_l = [favus_path, teeze_path, elupy_path, xavyr_path]

trac_x,trac_y = tracon_points(tracon, points)
cent_x, cent_y = tracon_points(center, points)
freeze_x, freeze_y = tracon_points(freeze, points)
other_x, other_y = tracon_points(other, points)
# print(freeze_y.index(np.max(freeze_y)), freeze_y.index(np.min(freeze_y)))

max_centx = np.max(cent_x)
min_centx = np.min(cent_x)
max_freezey = np.max(freeze_y)
min_freezey = np.min(freeze_y)
# for i in range(0, len(x)):
#     if x[i] <=0:
#         print('negtive',i)
#     elif i > 0 and x[i-1] <=0 and x[i] >=0:
#         print('positive after negative',i)
picked = []
ntraf = 20
for i in range(0, ntraf):
    lat, lon, alt, hdg, picked = get_ll(other_x, other_y, kcmh, picked)
    # print(hdg)
    bs.traf.cre(f'AC{i}', 'M250', lat, lon, hdg, 80, 57.412)
    print(bs.traf.id[i],bs.traf.hdg[i])
print('PICKED VALUES:', picked)
for acid in bs.traf.id:
        bs.stack.stack(f'BANK {acid} 73')
        bs.stack.stack(f'CONFLICTDETECTION {acid} OFF')
bs.stack.stack(f'ASAS OFF')
trac_wpts, tac_lat, trac_lon = Meters_To_WSG84([trac_x,trac_y], kcmh)
cent_wpts, cent_lat, cent_lon = Meters_To_WSG84([cent_x, cent_y], kcmh)
freeze_wpts, freeze_lat, freeze_lon = Meters_To_WSG84([freeze_x, freeze_y], kcmh)
# print(cent_lat)
max_centlon = np.max(cent_lon)
min_centlon = np.min(cent_lon)
max_freezelat = np.max(freeze_lat)
min_freezelat = np.min(freeze_lat)
# print(min_centlon, max_centlon, min_freezelat, max_freezelat)
# print(circ_wpts)
utm_zone = get_utm_zone(-82.2)
p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)





t_max = 30001
t_max = 6010
# t_max = 2
dt = 0.1
bs.stack.stack(f'DT {dt}')
# we'll run the simulation for up to 4000 seconds
# t_max = 1
# ntraf = bs.traf.ntraf
# print(ntraf)
t_tran = np.zeros((ntraf, 2))
# print(t_tran)
# t_tran[0] = [5, 10]
# print(t_tran, t_tran[0])
gate = np.zeros((ntraf, 1), dtype = str)

n_steps = int(t_max + 1)
# print(n_steps)
t = np.linspace(0, t_max, n_steps)
# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))
current_xy = np.zeros((ntraf, n_steps, 2))
# print(current_xy[0][0])
# print(res[0])
sorter = np.zeros((ntraf, 1))
c = 1
sort_group1 = []
sort_group2 = []
step = [0, 0]
closed = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
guide = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
waypts = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
d = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
d2 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
idx = [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
etas = []
fcfs = []

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
    if i == 3000:
        print(f'Optimization Interval Closed! Number Of Aircraft In Interval Is {len(locals()[f'sort_group{c}'])}')
        step[c] = 1
        c+=1
        

    for j in range(len(bs.traf.id)):
        pos = __WSG84_To_Meters_Single([bs.traf.lat[j], bs.traf.lon[j]], kcmh, p)
        current_xy[j][i] = pos
        # if i %100 == 0:
        #     print(f'min_centx: {min_centx}, p[0]: {pos[0]}, max_centx: {max_centx}, min_freezey: {min_freezey}, p[1]: {pos[1]}, max_freezey: {max_freezey}')
        # if min_centlon < bs.traf.lon[j] < max_centlon and min_freezelat < bs.traf.lat[j] < max_freezelat and sorter[j] ==  0:
        h, dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
        # print(dist, center)
        # if min_centx <= pos[0] <= max_centx and min_freezey <= pos[1] <= max_freezey and sorter[j] ==  0:
        if dist <= center and sorter[j] == 0:
            # print(f'{bs.traf.id[j]} crossed Optimization Horizon at timestep {i}!')
            # print(f'Postion is {bs.traf.lat[j], bs.traf.lon[j]}, horizon bounds are\nLON: {min_centlon},{max_centlon} LAT: {min_freezelat},{max_freezelat}')
            print(f'{bs.traf.id[j]} crossed Optimization Horizon at timestep {i}!')
            print(f'Postion is {pos[0], pos[1]}, horizon bounds are\nX: {min_centx, max_centx} Y: {min_freezey, max_freezey}')
            print(f'Postion is {bs.traf.lat[j], bs.traf.lon[j]}, horizon bounds are\nLON: {min_centlon,max_centlon} LAT: {min_freezelat,max_freezelat}\n')
            sorter[j] = 1
            h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
            eta = eta_dist/bs.traf.tas[j]
            locals()[f'sort_group{c}'].append([bs.traf.id[j], bs.traf.lat[j], bs.traf.lon[j], bs.traf.hdg[j], eta_dist, eta, i])
            fcfs.append(bs.traf.id[j])
            closed[j] = c
        if closed[j] == c-1:
            gate[j], g = find_quadrant(picked[j])
            h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
            eta = eta_dist/bs.traf.tas[j]
            t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas)
            etas.append([eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], eta])
            closed[j]+=1
        if step[c-1] == 2 and guide[j] == -1:
            index, waypts[j] = findWPs(bs.traf.id[j], run, fixes_l)
            h, d[j] = qdrdist(waypts[j][0][0], waypts[j][1][0], bs.traf.lat[j], bs.traf.lon[j])
            idx[j] = 0
            guide[j] = 1
        if guide[j] == 1:
            pos = [bs.traf.lat[j], bs.traf.lon[j]]
            # print(bs.traf.id[j], d[j])
            # m_pos = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
            # print(waypts[j][1])
            # print(m_pos[1])
            # print(waypts[j])
            head, d[j], idx[j], stat = WPGuidance(waypts[j], idx[j], [bs.traf.lat[j], bs.traf.lon[j]], d[j], gate[j], i, bs.traf.id[j])
            # print(d[j])
            # head = np.arctan2(waypts[j][1][index[j]] - pos[1], waypts[j][0][index[j]]-pos[0])
            h, dtw = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
            d2[j] = [bs.traf.id[j], dtw]
            bs.stack.stack(f'HDG {bs.traf.id[j]} {90-head}')
        # if  step[c-1] ==2 and i%1000 == 0:
        #     gate[j], g = find_quadrant(picked[j])
        #     t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas)
        #     etas.append([eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], eta])
    if step[c-1] == 1:
        t_it = 5
        dtmax = 8
        sta_p, run, cp3, staff = SQO.sort(etas, t_it, t_tran, dtmax)
        step[c-1] = 2
        etas = []
    # if step [c-1] == 2 and i%1000 == 0 and etas:
    #     sta_p, run, cp3, staff = SQO.sort(etas, t_it, t_tran, dtmax)
    #     etas = []
    #     step[c-1]=2
    #     fcfs, order = SQO.step0(etas)
    # if step[c-1] == 2:
    #     # t_tran = [[15, 30], [30, 15], [15, 30], [30, 15]]
        
    #     print(f'AIRCRAFT TRANSITION TIMES: {t_tran}')
    #     t_it = 5
    #     cp, rta_p, non_rta = SQO.step12(fcfs, t_it, t_tran)
    #     step[c-1]=3
    # if step[c-1] == 3:
    #     sta_p, non_sta = SQO.step3(cp, rta_p, non_rta, order, t_it)
    #     step[c-1]=4
    # if step[c-1] == 4:
    #     dtmax = 8
    #     dels, ddfC, ddfT = SQO.step4(sta_p, dtmax)
    #     step[c-1] = 5
    # if step[c-1] == 5:
    #     SQO.step5()
            
print(run)
run_id = [i[1] for i in run]
print(run_id)
# print(sorted(d2, key=lambda x: x[1]))
            


end_time = time.time()
print("SIM TIME", end_time-start_time)
# print(sort_group1)
'''Plotting'''

for idx, acid in enumerate(bs.traf.id):

    if idx%2 == 0:
        marker = 'o'
    elif idx%3 == 0:
        marker = '^'
    else:
        marker = 'x'

    plt.plot(res[::10, 1, idx], res[::10, 0, idx], marker = marker, linestyle = '', label=f'AC{idx}')
    plt.plot([kcmh[1], res[0, 1, idx]], [kcmh[0], res[0, 0, idx]], linestyle = '--', color = 'black', linewidth = 2)

for i in range(0, points, 2):
    if i == 0:
        plt.scatter(trac_wpts[i][1], trac_wpts[i][0], color = 'red', label = 'TRACON AIRSPACE', marker = 's')
        plt.scatter(cent_wpts[i][1], cent_wpts[i][0], color = 'blue', label = 'OPTIMIZATION HORIZON', marker = 's')
        plt.scatter(freeze_wpts[i][1], freeze_wpts[i][0], color = 'green', label = 'FREEZE HORIZON', marker = 's')
    else:
        plt.scatter(trac_wpts[i][1], trac_wpts[i][0], color = 'red', marker = 's')
        plt.scatter(cent_wpts[i][1], cent_wpts[i][0], color = 'blue', marker = 's')
        plt.scatter(freeze_wpts[i][1], freeze_wpts[i][0], color = 'green', marker = 's')
plt.scatter(kcmh[1], kcmh[0], color = 'black', s = 100, marker = 's', label = 'KCMH')
plt.scatter(xavyr[1], xavyr[0], color = 'cyan', s = 100, marker = 'o', label = 'XAVYR')
plt.scatter(elupy[1], elupy[0], color = 'magenta', s = 100, marker = 'o', label = 'ELUPY')
plt.scatter(favus[1], favus[0], color = 'brown', s = 100, marker = 'o', label = 'FAVUS')
plt.scatter(dubln[1], dubln[0], color = 'orange', s = 100, marker = 'o', label = 'DUBLN')
plt.scatter([favus_path[0][1], favus_path[1][1]], [favus_path[0][0], favus_path[1][0]], marker = 's', s = 100, color = 'brown')
plt.scatter(xavyr_path[0][1], xavyr_path[0][0], marker = 's', s = 100, color = 'cyan')
plt.scatter([teeze_path[0][1], teeze_path[2][1], teeze_path[3][1], teeze_path[4][1], teeze_path[5][1]], 
            [teeze_path[0][0], teeze_path[2][0], teeze_path[3][0], teeze_path[4][0], teeze_path[5][0]], marker = 's', s = 100, color = 'orange')
plt.scatter([elupy_path[0][1], elupy_path[1][1], elupy_path[3][1], elupy_path[4][1], elupy_path[5][1]],
            [elupy_path[0][0], elupy_path[1][0], elupy_path[3][0], elupy_path[4][0], elupy_path[5][0]], marker = 's', s=  100, color = 'magenta')
plt.axis('equal')
# plt.gca().set_aspect('equal', adjustable='datalim')

plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.xlim(-83.3, -82.5)
plt.legend(loc = 'upper right')
plt.show()