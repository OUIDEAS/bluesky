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
import BezIntersectOptimizeP1P2 as BIO
import sympy as sp
from scipy.optimize import fsolve

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

def reformat_list(data):
    reformatted = []
    for sublist in data:
        if sublist != False:
            x_vals = [point[0] for point in sublist]
            y_vals = [point[1] for point in sublist]
            reformatted.append([x_vals, y_vals])
        else:
            reformatted.append(False)
    return reformatted

def tracon_points(r, n):

    x = [math.cos(2*pi/n*x)*r for x in range(0,n+1)]
    y = [math.sin(2*pi/n*x)*r for x in range(0,n+1)]

    return x, y

def get_ll(x, y, home_ll, picked):
    i = np.random.randint(0, 360)
    j = np.random.randint(500, 3000)
    c = 0 
    # i = 202
    # j = 500
    while c == 0:
        if any(abs(i - p) <= 3 for p in picked):
            i = np.random.randint(0, 360)
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
    # diff = correct-hdg
    # hdg+=diff
    # if i <= 180:
    hdg=correct+180
    # else:
        # hdg-=180
    if hdg >=360:
        hdg-=360

    lat, lon, alt = Meters_To_WSG84([xp, yp], home_ll)
    return lat, lon, alt, hdg, picked

def find_quadrant(picked):
    if picked<=90:
        quad = 'I'
    elif 90<picked<=180:
        quad = 'II'
    elif 180<picked<=270:
        quad = 'III'
    else:
        quad = 'IV'
    # print('in GATE FUNCTION')
    return quad, quad

def get_eta(gate, lat, lon, eta, fixes, v):
    # fixes -> R, 1,2,2,1
     # fixes_l = [favus_path, teeze_path, elupy_path, xavyr_path]
    # print(len(lat))
    for i in range(len(lat)):
        # print(gate[i][0], i)
        if gate[i][0] == 'I' or gate[i][0] == 'II':
            h, d_r1 = qdrdist(lat[i], lon[i], fixes[0][0], fixes[0][1])
            h2, d_r2 = qdrdist(lat[i], lon[i], fixes[1][0], fixes[1][1])
            t_r1 = d_r1/v[i]
            t_r2 = d_r2/v[i]
            # print(t_r1, t_r2)
            eta[i] = [t_r1, t_r2]

        elif gate[i][0] == 'IV' or gate[i][0] == 'III':
            h2, d_r1 = qdrdist(lat[i], lon[i], fixes[3][0], fixes[3][1])
            h, d_r2 = qdrdist(lat[i], lon[i], fixes[2][0], fixes[2][1])
            t_r1 = d_r1/v[i]
            t_r2 = d_r2/v[i]
            eta[i] = [t_r1, t_r2]

    # print(t_tran)
    # print(eta)
    return eta

def get_tran(gate, lat, lon, t_tran, fixes, v, star_routes):
    # fixes -> R, 1,2,2,1
    # fixes_l = [favus_path, teeze_path, elupy_path, xavyr_path]
    for i in range(len(lat)):
        if gate[i] == 'I' or gate[i] == 'II':
            # print(gate, gate[i], star_routes)
            d_r1 = star_routes[0][2]+star_routes[2][3]
            d_r2 = star_routes[1][1]+star_routes[1][2]+star_routes[1][3]+star_routes[1][4]
            t_r1 = d_r1/v[i]
            t_r2 = d_r2/v[i]
            t_tran[i] = [t_r1, t_r2]

        elif gate[i] == 'III' or gate[i] == 'IV':
            d_r2 = star_routes[2][2]+star_routes[2][3]+star_routes[2][4]
            d_r1 = star_routes[3][1]+star_routes[3][2]
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
            wps = fixes_l[1]
        else:
            print('FAVUS')
            wp_list = makeWPTS(fixes_l[0])
            wps = fixes_l[0]
    if seq_dat[3] == 'I':
        if seq_dat[4] == 'R1':
            print('FAVUS')
            wp_list = makeWPTS(fixes_l[0])
            wps = fixes_l[0]
        else:
            print('TEEZE')
            wp_list = makeWPTS(fixes_l[1])
            wps = fixes_l[1]
    if seq_dat[3] == 'III':
        if seq_dat[4] == 'R2':
            print('ELUPY')
            wp_list = makeWPTS(fixes_l[2])
            wps = fixes_l[2]
        else:
            print('XAVYR')
            wp_list = makeWPTS(fixes_l[3])
            wps = fixes_l[3]
    if seq_dat[3] == 'IV':
        if seq_dat[4] == 'R1':
            print('XAVYR')
            wp_list = makeWPTS(fixes_l[3])
            wps = fixes_l[3]
        else:
            print('ELUPY')
            wp_list = makeWPTS(fixes_l[2])
            wps = fixes_l[2]
    # print(wp_list)
    return index, wp_list, wps

def WPGuidance(wp_list, idx, pos, d, gate, i, id, h_current, lr):
    # print(wp_list[idx])
    turn = np.deg2rad(73)
    stat = 'not_done'
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    kcmh = [39.997472, -82.891194]
    m_pos = __WSG84_To_Meters_Single(pos, kcmh, p)
    # print(wp_list[0])
    if idx >= len(wp_list[0]):
        idx-=1
        stat = 'done'
    wp_pos = __WSG84_To_Meters_Single([wp_list[0][idx], wp_list[1][idx]], kcmh, p)
    head = np.rad2deg(np.arctan2(wp_pos[1]-m_pos[1], wp_pos[0]-m_pos[0]))
    h, dn = qdrdist(pos[0], pos[1], wp_list[0][idx], wp_list[1][idx])
    # head = np.pi/2-(np.arctan2(y_int-center[1], x_int-center[0]))
    # rt_cmd = np.clip((head-np.deg2rad(90-hdg))/dt, -turn, turn)
    # head = head+lr*rt_cmd*dt
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
    # # print(head)
    # print(f'{id}, {dn}, {d}, {idx}')
    if dn > d and dn <= 50:
        idx+=1
        # if head >= 180 + h_current or head<= h_current-180:
        #     idx-=1
        if idx > len(wp_list[0]):
            return head, dn, idx-1, stat
    return head, dn, idx, stat
    
def interp_wps(waypts, n, v):
    interps = np.zeros((len(waypts)-1, 2, n))
    dists = []
    travel = []
    # interps[0][0] = (8, 10)
    # print(interps[0][0])
    # interps[0][0][0] = 9
    # interps[0][1][0] = 10
    # print(interps[0])
    # print(interps)
    total = 0
    for i in range(0, len(waypts)-1):
        dx = waypts[i+1][0] - waypts[i][0]
        dy = waypts[i+1][1] - waypts[i][1]
        m = (waypts[i+1][1] - waypts[i][1])/(waypts[i+1][0] - waypts[i][0])
        d = np.hypot(dx, dy)
        dists.append(d)
        travel.append(d/v)
        total+= d/v
        x = np.linspace(waypts[i][0], waypts[i+1][0], n, endpoint=False)
        interps[i][0][0], interps[i][1][0] = x[0], waypts[i][1]
        for j in range(1, len(x)):
            y = waypts[i][1] + m * (x[j]-waypts[i][0])
            interps[i][0][j] = x[j]
            interps[i][1][j] = y
    print(f'Total STAR Route Travel Time: {total}')
    return interps, dists, travel, total

def carrot_following(nodes, ac_ll, kcmh, tr, t_guess, hdg, dt, lr, gate):
    turn = np.deg2rad(73)
    # print(lr)
    # if hdg<0:
    #     hdg+=360
    # print(h, hdg)
    t_vals = np.linspace(0,1,500)
    
    bx_vals = nodes[0][0] * (1 - t_vals)**2 + 2 * (1 - t_vals) * t_vals * nodes[0][1] + t_vals**2 * nodes[0][2]
    by_vals = nodes[1][0] * (1 - t_vals)**2 + 2 * (1 - t_vals) * t_vals * nodes[1][1] + t_vals**2 * nodes[1][2]

    center = LongLat_To_WSG84_Meters(ac_ll, kcmh)
    y_c = center[1]
    x_c = center[0]
    t_guess +=0.0525

    Bx = lambda t: nodes[0][0] * (1 - t)**2 + 2 * (1 - t) * t * nodes[0][1] + t**2 * nodes[0][2]
    By = lambda t: nodes[1][0] * (1 - t)**2 + 2 * (1 - t) * t * nodes[1][1] + t**2 * nodes[1][2]
    circle_eq = lambda t: (Bx(t)-x_c)**2+(By(t)-y_c)**2 - (tr)**2
    
  
    sol = fsolve(circle_eq, t_guess)

    sol[0] = np.clip(sol[0], 0, 1)
    t_new = sol[0]
    x_int = Bx(t_new)
    y_int = By(t_new)


    head = np.pi/2-(np.arctan2(y_int-center[1], x_int-center[0]))
    # rt_cmd = np.clip((head-np.deg2rad(90-hdg))/dt, -turn, turn)
    # head = head+lr*rt_cmd*dt
    
    # fig, ax = plt.subplots(figsize=(8, 8))
    # plt.axis('equal')
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.scatter(center[0], center[1], color = 'black')
    # plt.plot(bx_vals, by_vals, color = 'blue')
    # plt.title(f'{t_new}, {hdg}, {np.rad2deg(head)}')
    # circle = plt.Circle((x_c, y_c), tr, color='r', fill=False, label="Circle")
    # ax.add_artist(circle)
    # plt.scatter(x_int, y_int)
    # plt.arrow(x_c, y_c, 2500*np.cos(np.deg2rad(90-hdg)), 2500*np.sin(np.deg2rad(90-hdg)), width=50,
    #             length_includes_head=False,
    #             head_width=1000,
    #             head_length=1000,
    #             head_starts_at_zero=False,
    #             facecolor='black',
    #             edgecolor='black',
    #             label = f'Heading of {hdg:0.2f}')
    # plt.show()
    return np.rad2deg(head), t_new

def get_wpIndex(pos, wp_list, gate):
    dists = 0
    tol = 0.005
    if gate =='I':
        for i in range(len(wp_list[0])-1):
            if wp_list[0][i+1]-tol<=pos[0]<=wp_list[0][i]+tol and pos[1] >= wp_list[1][i+1]:
                return i+1
        return len(wp_list[0])-1
    elif gate == 'II':
        for i in range(len(wp_list[0])-1):
            if wp_list[0][i+1]-tol<=pos[0]<=wp_list[0][i]+tol and pos[1] <= wp_list[1][i+1]:
                return i+1
        return len(wp_list[0])-1
    elif gate == 'III':
        for i in range(len(wp_list[0])-1):
            if i ==4:
                if wp_list[0][i+1]-tol<=pos[0]<=wp_list[0][i]+tol and pos[1] <= wp_list[1][i+1]:
                    return i+1
            if wp_list[0][i]-tol<=pos[0]<=wp_list[0][i+1]+tol and pos[1] <= wp_list[1][i+1]:
                return i+1
        return len(wp_list[0])-1
    elif gate =='IV':
        for i in range(len(wp_list[0])-1):
            if wp_list[0][i]-tol<=pos[0]<=wp_list[0][i+1]+tol and pos[1] >= wp_list[1][i+1]:
                return i+1
        return len(wp_list[0])-1
        
    


start_time = time.time()
bs.init(mode ='sim', detached=True)
bs.scr = ScreenDummy()

#Pattern for pulling from terminal
pattern = r'Dist = (\d+\.\d+) nm'
#Create Log
bs.stack.stack(f'CRELOG SequencingFlightTest3 1')
bs.stack.stack(f'SCEN SequencingFlightTest3; SAVEIC SequencingFlightTest3')
# bs.stack.stack(f'CRELOG HoldSingleBypass 1')
# bs.stack.stack(f'SCEN HoldSingleBypass; SAVEIC HoldSingleBypass')
# Set Sim Time

#Read from datasets
wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
datasets_string = ['wptdata', 'aptdata']
datasets = [wptdata, aptdata]

# bs.traf.mcre(10, 'M250')
#R2 is 28 L
nm = 1852
points = 360
tracon = 10*nm
center = 30*nm # OPTIMIZATION HORIZON
freeze = 20.7*nm
other = 36*nm

scen = 'Bez'

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
teeze_wps = LongLat_To_WSG84_Meters(teeze_path, kcmh)
# teeze_pts, teeze_dists, teeze_toa, teeze_total_t = interp_wps(teeze_wps, 100, 57.412)
teeze_pts, teeze_dists, teeze_toa, teeze_total_t = interp_wps(teeze_wps, 100, 54.84)



# FAVUS WAYPOINTS
bugzz = [40.565, -82.454056]
cbuss = [40.325306, -82.5405]
molls = [40.132139, -82.609194]
ordiw = [39.989306, -82.666056]
bazel = [39.981917, -82.704556]

favus_path =[cbuss, molls, favus, ordiw, bazel]
favus_wps = LongLat_To_WSG84_Meters(favus_path, kcmh)
# favus_pts, favus_dists, favus_toa, favus_total_t = interp_wps(favus_wps, 100, 57.412)
favus_pts, favus_dists, favus_toa, favus_total_t = interp_wps(favus_wps, 100, 54.84)

# ELUPY WAYPOINTS
jaktz = [39.591028, -83.419583]
rscot = [39.722389,  -83.286306]
obetz = [39.787667, -83.158389]
edwib = [39.877472, -82.984861]
gagbe = [39.907167, -82.927278]
jesce = [39.903556, -82.858889]

elupy_path = [rscot, obetz, elupy, edwib, gagbe, jesce]
elupy_wps = LongLat_To_WSG84_Meters(elupy_path, kcmh)
# elupy_pts, elupy_dists, elupy_toa, elupy_total_t = interp_wps(elupy_wps, 100, 57.412)
elupy_pts, elupy_dists, elupy_toa, elupy_total_t = interp_wps(elupy_wps, 100, 54.84)

# XAVYR WAYPOINTS
scrlt = [39.502917, -82.350833]
brtus = [39.730944, -82.473083]
guber = [39.963222, -82.670889]
bkeye = [39.982056, -82.706694]

xavyr_path = [brtus, xavyr, guber, bkeye]
xavyr_wps = LongLat_To_WSG84_Meters(xavyr_path, kcmh)
# xavyr_pts, xavyr_dists, xavyr_toa, xavyr_total_t = interp_wps(xavyr_wps, 100, 57.412)
xavyr_pts, xavyr_dists, xavyr_toa, xavyr_total_t = interp_wps(xavyr_wps, 100, 54.84)

fixes = [favus, dubln, elupy, xavyr]
star_routes = [favus_dists, teeze_dists, elupy_dists, xavyr_dists] 

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
    # bs.traf.cre(f'AC{i}', 'M250', lat, lon, hdg, 10000, 57.412)
    bs.traf.cre(f'AC{i}', 'M250', lat, lon, hdg, 5000, 54.84)
    print(bs.traf.id[i],bs.traf.hdg[i])
print('PICKED VALUES:', picked)
for acid in bs.traf.id:
        bs.stack.stack(f'BANK {acid} 73')
        bs.stack.stack(f'CONFLICTDETECTION {acid} OFF')
        bs.stack.stack(f'VS {acid} 2.2')
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



b_ll = []
e_ll = []

t_max = 15001
# t_max = 10000
# t_max = 1

# t_max = 6010
# t_max = 2
dt = 0.1
bs.stack.stack(f'DT {dt}')
# we'll run the simulation for up to 4000 seconds
# t_max = 1
# print(ntraf)
eta = np.zeros((ntraf, 2)) # AC -> Each Meter Gate
t_tran = np.zeros((ntraf, 2))
# print(t_tran)
# t_tran[0] = [5, 10]
# print(t_tran, t_tran[0])
gate = np.zeros((ntraf, 1), dtype = '<U3')

n_steps = int(t_max + 1)
# print(n_steps)
t = np.linspace(0, t_max, n_steps)
# allocate some empty arrays for the results
res = np.zeros((n_steps, 5, ntraf))
current_xy = np.zeros((ntraf, n_steps, 2))
# print(current_xy[0][0])
# print(res[0])
sorter = np.zeros((ntraf, 1))
t_it = 10
c = 1
sort_group1 = []
sort_group2 = []
step = [0, 0]
closed = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
guide = [0] * ntraf#np.zeros(ntraf)#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
waypts = [-1] * ntraf# [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
waypt_set = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
bez = [-1] * ntraf#np.zeros(ntraf)#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
entry = [[-1, -1, -1] for _ in range(ntraf)]#np.zeros(ntraf)#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
interps = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
end_point = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
nodes = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
times = [-1] * ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
diffs = [[-1,-1] for _ in range(ntraf)]#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
entry_point = [[-1,-1] for _ in range(ntraf)]#np.zeros((ntraf,2))#[]
t_g = [0]*ntraf
lr = [0]*ntraf
fl = [0]*ntraf
re_sched = []
# exit_point = np.zeros((20,2))#[]
start_end = [[0,0,0] for _ in range(ntraf)]
# start_end[0][1] = 3
# print(start_end[0], start_end[1], start_end)

d = np.zeros(ntraf)#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
d2 = [0,0]*ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
index = [0]*ntraf#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
run_id = np.zeros(ntraf)#[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
etas = []
fcfs = []
count = [0]*ntraf
final_order = []
passed = [False]*ntraf
ac_info = [[0, 0, 0, 0, 0, 0, True, 0, 0, 0] for _ in range(ntraf)]
plan_eta = [0 for _ in range(ntraf)]
diff = [0]*ntraf
buff = [0]*ntraf
all_dists = np.zeros((n_steps, ntraf, ntraf), dtype = int)

start_time = time.time()
if scen == 'WP':
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
            for k in range(ntraf):
                if k != j:
                    h, d = qdrdist(bs.traf.lat[j], bs.traf.lon[j], bs.traf.lat[k], bs.traf.lon[k])
                    all_dists[i][j][k] = d
            # if i %100 == 0:
            #     print(f'min_centx: {min_centx}, p[0]: {pos[0]}, max_centx: {max_centx}, min_freezey: {min_freezey}, p[1]: {pos[1]}, max_freezey: {max_freezey}')
            # if min_centlon < bs.traf.lon[j] < max_centlon and min_freezelat < bs.traf.lat[j] < max_freezelat and sorter[j] ==  0:
            h, dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
            # print(dist, center)
            print(dist, center, bs.traf.id[j])
            # if min_centx <= pos[0] <= max_centx and min_freezey <= pos[1] <= max_freezey and sorter[j] ==  0:
            if dist <= center and sorter[j] == 0:
                # print(f'{bs.traf.id[j]} crossed Optimization Horizon at timestep {i}!')
                # print(f'Postion is {bs.traf.lat[j], bs.traf.lon[j]}, horizon bounds are\nLON: {min_centlon},{max_centlon} LAT: {min_freezelat},{max_freezelat}')
                print(f'{bs.traf.id[j]} crossed Optimization Horizon at timestep {i}!')
                print(f'Postion is {pos[0], pos[1]}, horizon bounds are\nX: {min_centx, max_centx} Y: {min_freezey, max_freezey}')
                print(f'Postion is {bs.traf.lat[j], bs.traf.lon[j]}, horizon bounds are\nLON: {min_centlon,max_centlon} LAT: {min_freezelat,max_freezelat}\n')
                sorter[j] = 1
                if j == 0:
                    gate[j], g = find_quadrant(picked[j])
                else:
                    gate[j], g = find_quadrant(picked[j+(4*j)])
                # print(gate[j], g, bs.traf.id[j])
                eta = get_eta(gate, bs.traf.lat, bs.traf.lon, eta, fixes, bs.traf.tas) #ETA now has the shape [[1,2], []...]
                # h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # eta = eta_dist/bs.traf.tas[j]
                locals()[f'sort_group{c}'].append([bs.traf.id[j], bs.traf.lat[j], bs.traf.lon[j], bs.traf.hdg[j], eta[j], i])
                fcfs.append(bs.traf.id[j])
                closed[j] = c
            if closed[j] == c-1:
                # gate[j], g = find_quadrant(picked[j])
                # h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # eta = eta_dist/bs.traf.tas[j]
                # t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas)
                # etas.append([eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], eta])
                # closed[j]+=1
                gate[j], g = find_quadrant(picked[j])
                # print(gate[j][0], g, bs.traf.id[j])
                # print('GATE GATE GATE', gate)
                # h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # eta = eta_dist/bs.traf.tas[j]
                eta = get_eta(gate, bs.traf.lat, bs.traf.lon, eta, fixes, bs.traf.tas) #ETA now has the shape [[1,2], []...]
                p_eta = np.min(eta[j])
                t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas, star_routes)
                # print(eta[j])
                # print(t_tran[j])
                etas.append([p_eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], p_eta, eta[j]])
                closed[j]+=1
            if step[c-1] == 2 and guide[j] == -1:
                index[j], waypts[j] = findWPs(bs.traf.id[j], run, fixes_l)
                h, d[j] = qdrdist(waypts[j][0][0], waypts[j][1][0], bs.traf.lat[j], bs.traf.lon[j])
                index[j] = 0
                guide[j] = 1
            if guide[j] == 1:
                pos = [bs.traf.lat[j], bs.traf.lon[j]]
                # print(bs.traf.id[j], d[j])
                # m_pos = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
                # print(waypts[j][1])
                # print(m_pos[1])
                # print(waypts[j])
                head, d[j], index[j], stat = WPGuidance(waypts[j], index[j], [bs.traf.lat[j], bs.traf.lon[j]], d[j], gate[j], i, bs.traf.id[j], bs.traf.hdg[j])
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
            # t_it = 10
            dtmax = 13
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
if scen == 'Bez':
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
            
        lat = np.array(bs.traf.lat)
        lon = np.array(bs.traf.lon)
        for j in range(len(bs.traf.id)):
            # bs.stack.stack(f'SPD {bs.traf.id[j]} 111.6')
            pos = __WSG84_To_Meters_Single([bs.traf.lat[j], bs.traf.lon[j]], kcmh, p)
            current_xy[j][i] = pos
            # for k in range(ntraf):
            #     if k != j:
            #         h, d = qdrdist(bs.traf.lat[j], bs.traf.lon[j], bs.traf.lat[k], bs.traf.lon[k])
            #         all_dists[i][j][k] = d
            mask = np.arange(ntraf) != j  # Create a boolean mask to exclude j
            h, dists = qdrdist(lat[j], lon[j], lat[mask], lon[mask])  # Vectorized call
            all_dists[i, j, mask] = dists  # Assign results directly
            all_dists[i][j][j] = -1
            # if j== 1 and guide[j]==4: 
            #     print(bs.traf.id[j], bs.traf.hdg[j], guide[j])
            # if guide[j] == 3:
                # print(f'{pos}, {entry_point[j]}, {bs.traf.id[j]}')
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
                gate[j], g = find_quadrant(picked[j])
                # print(gate[j], g, bs.traf.id[j])
                eta = get_eta(gate, bs.traf.lat, bs.traf.lon, eta, fixes, bs.traf.tas) #ETA now has the shape [[1,2], []...]
                # h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # eta = eta_dist/bs.traf.tas[j]
                locals()[f'sort_group{c}'].append([bs.traf.id[j], bs.traf.lat[j], bs.traf.lon[j], bs.traf.hdg[j], eta[j], i])
                fcfs.append(bs.traf.id[j])
                closed[j] = c
            if closed[j] == c-1:
                # gate[j], g = find_quadrant(picked[j])
                # h, eta_dist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # eta = eta_dist/bs.traf.tas[j]
                # t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas)
                # etas.append([eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], eta])
                # closed[j]+=1
                gate[j], g = find_quadrant(picked[j])

                eta = get_eta(gate, bs.traf.lat, bs.traf.lon, eta, fixes, bs.traf.tas) #ETA now has the shape [[1,2], []...]
                p_eta = np.min(eta[j])
                t_tran = get_tran(gate, bs.traf.lat, bs.traf.lon, t_tran, fixes_l, bs.traf.tas, star_routes)
                # print(eta[j])
                # print(t_tran[j])
                etas.append([p_eta, bs.traf.id[j], fcfs.index(bs.traf.id[j]), g, t_tran[j], p_eta, eta[j]])
                closed[j]+=1
            if step[c-1] == 2 and guide[j] == 0:
                index[j], waypts[j], waypt_set[j] = findWPs(bs.traf.id[j], run, fixes_l)
                h, d[j] = qdrdist(waypts[j][0][0], waypts[j][1][0], bs.traf.lat[j], bs.traf.lon[j])
                index[j] = 0
                guide[j] = 1
            if guide[j] == 1:
                pos_ll = [bs.traf.lat[j], bs.traf.lon[j]]
                # print(bs.traf.id[j], pos_ll, bs.traf.hdg[j])
                interps[j], end_point[j], nodes[j], bez[j], entry[j], fl[j], times[j], diffs[j], lr[j], path = BIO.GenerateTrajectoryOutside([bs.traf.lat[j], bs.traf.lon[j]],
                                                                                                                                        run_id[j][0], waypt_set[j], 
                                                                                                                                        bs.traf.hdg[j], bs.traf.tas[j], kcmh, 
                                                                                                                                        gate[j], False, ac_info[j], t_it)
                b_ll.append(path[0])
                # if path[1] != False:
                e_ll.append(path[1])
                plan_eta[j] = diffs[j][1]
                diff[j] = diffs[j][0]
                if entry[j][0]:
                    # entry_point[j]= Meters_To_WSG84([entry[j][4][0], entry[j][4][1]], kcmh)
                    entry_point[j]= [entry[j][4][0], entry[j][4][1]]
                    if entry[j][2]<0:
                        entry[j][2]+=360
                guide[j] = 2
                start_end[j][0] = i
                if bs.traf.id[j] == bs.traf.id[-1]:
                    step[c-1] = 3
            if buff[j] == 1:
                print(f'REPLANNING PATH FOR {bs.traf.id[j]}')
                pos_ll = [bs.traf.lat[j], bs.traf.lon[j]]
                # print(bs.traf.id[j], pos_ll, bs.traf.hdg[j])
                interps[j], end_point[j], nodes[j], bez[j], entry[j], fl[j], times[j], diffs[j], lr[j], path = BIO.GenerateTrajectoryOutside([bs.traf.lat[j], bs.traf.lon[j]],
                                                                                                                                        ac_info[j][1], waypt_set[j], 
                                                                                                                                        bs.traf.hdg[j], bs.traf.tas[j], kcmh, 
                                                                                                                                        gate[j], False, ac_info[j], t_it)
                b_ll[j] = path[0]
                # if path[1] != False:
                e_ll[j] = path[1]
                plan_eta[j] = diffs[j][1]
                diff[j] = diffs[j][0]
                if entry[j][0]:
                    # entry_point[j]= Meters_To_WSG84([entry[j][4][0], entry[j][4][1]], kcmh)
                    entry_point[j]= [entry[j][4][0], entry[j][4][1]]
                    if entry[j][2]<0:
                        entry[j][2]+=360
                guide[j] = 2
                buff[j] = 2
                # start_end[j][0] = i
            if guide[j]==2 and entry[j][0] and buff[j] == 2:
                # if picked[j]==0 or picked[j] == 360: #np.isclose(picked[j], 0, atol = 2):
                #     bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[j][3]}; HDG {bs.traf.id[j]} {entry[j][2]+180}')
                #     print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} WITH A BANK OF {entry[j][3]} UNTIL A HEADING OF {entry[j][2]+180}')
                # else:
                # bs.stack.stack(f'ALT {bs.traf.id[j]}, {fl[j]}')
                bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[j][3]}; HDG {bs.traf.id[j]} {entry[j][2]}')
                print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} WITH A BANK OF {entry[j][3]} UNTIL A HEADING OF {entry[j][2]}')
                guide[j]=3
            # if guide[j] == 3 and entry[j][0]:
            #     print(f'HEADING {bs.traf.hdg[j]}, POS: {pos}, ENTRY POINT: {entry_point[j]}, DIST: {np.hypot(entry_point[j][1]-pos[1], entry_point[j][0]-pos[0])}')

            if np.isclose(entry_point[j][0], pos[0], atol = 50) and np.isclose(entry_point[j][1], pos[1], atol = 50) and guide[j]==3:
                guide[j] = 4
                # bs.stack.stack(f'ALT {bs.traf.id[j]}, {fl[j]}')
                bs.stack.stack(f'BANK {bs.traf.id[j]} 73')

                print(f'HEADNIG AT INTERCEPTION: {bs.traf.hdg[j]}, TARGET: {entry[j][2]}, DIFF: {np.abs(bs.traf.hdg[j]-entry[j][2])}')
                print(f'INITIATING NONLINEAR GUIDANCE FOR {bs.traf.id[j]} AT TIME STEP: {i}')
            elif guide[j] == 2 and not entry[j][0]:
                guide[j] = 4
                bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                print(f'INITIATING NONLINEAR GUIDANCE FOR {bs.traf.id[j]} AT TIME STEP: {i}')
            if guide[j]==4:

                pos_ll = [bs.traf.lat[j], bs.traf.lon[j]]
                tr = 0.3048*((bs.traf.tas[j]*1.94384)**2/(11.26*math.tan(np.deg2rad(73))))

                head, t_g[j] = carrot_following(nodes[j], pos_ll, kcmh, tr, t_g[j], bs.traf.hdg[j], dt, lr[j], gate[j])

                bs.stack.stack(f'HDG {bs.traf.id[j]} {head}')
                if t_g[j] >=1:
                    if gate[j] in ['I', 'IV'] and pos[0] <= bez[j][0][-1]:
                        passed[j] = True
                    elif gate[j] in ['II', 'III'] and pos[0] >= bez[j][0][-1]:
                        passed[j] = True

            
            if t_g[j] >= 1 and guide[j] == 4 and passed[j]:
                guide[j] = 5
                print(f'INITIATING WP GUIDANCE FOR {bs.traf.id[j]} AT TIMESTEP {i}')
                tag = i
               
                h, d[j] = qdrdist(waypts[j][0][index[j]], waypts[j][1][index[j]], bs.traf.lat[j], bs.traf.lon[j])
                index[j] = get_wpIndex([bs.traf.lat[j], bs.traf.lon[j]], waypts[j], gate[j])

                
                print(f'INITIAL WP INDEX IS {index[j]}')

            if guide[j] == 5:
                
                head, d[j], index[j], stat = WPGuidance(waypts[j], index[j], [bs.traf.lat[j], bs.traf.lon[j]], d[j], gate[j], i, bs.traf.id[j], bs.traf.hdg[j], lr[j])
                # if gate[j] == 'I':
                #     if index[j] == 0:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 8000')
                #         # bs.stack.stack(f'VS {bs.traf.id[j]} 1128.2')
                #     if index[j] >=1:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 5000')
                # if gate[j] == 'II': #[melzz, dubln, trlgy, polrs, taces, teeze]
                #     if index[j] == 0:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 10000')
                #     if index[j] == 2:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 8000')
                #         # bs.stack.stack(f'VS {bs.traf.id[j]} 1128.2')
                #     if index[j] >= 3:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 6000')
                # if gate[j] == 'III':#[rscot, obetz, elupy, edwib, gagbe, jesce]
                #     if index[j] == 0:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 10000')
                #     if index[j] == 2:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 8000')
                #         # bs.stack.stack(f'VS {bs.traf.id[j]} 1128.2')
                #     if index[j] == 3:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 7000')
                #         # bs.stack.stack(f'VS {bs.traf.id[j]} 1130.5')
                #     if index[j] >= 4:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 6000')
                # if gate[j] == 'IV': #[brtus, xavyr, guber, bkeye]
                #     if index[j] == 0:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 10000')
                #     if index[j] >= 1:
                #         bs.stack.stack(f'ALT {bs.traf.id[j]}, 5000')

                h, dtw = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                d2[j] = [bs.traf.id[j], dtw]

                bs.stack.stack(f'HDG {bs.traf.id[j]} {90-head}')
                # if index[j] >= len(waypts[j][0]):
                    # print(f'{bs.traf.id[j]}, {start_end[j]}, {stat}, {index[j]}, {len(waypts[j][0])}')
                if stat == 'done' and start_end[j][2] == 0:
                    start_end[j][1] = i
                    start_end[j][2] = 1
                    final_order.append([bs.traf.id[j], (start_end[j][1] - start_end[j][0])*dt, run_id[j][0], (start_end[j][1] - start_end[j][0])*dt - run_id[j][0], times[j], gate[j]])
                # h, dtw = qdrdist(bs.traf.lat[j], bs.traf.lon[j], kcmh[0], kcmh[1])
                # d2[j] = [bs.traf.id[j], dtw]
                # bs.stack.stack(f'HDG {bs.traf.id[j]} {90-head}')
            # if guide[j] >= 3 and i %10 == 0:
            #     plt.scatter(bs.traf.lon[j], bs.traf.lat[j], color = 'black')
        if step[c-1] == 1:
            # t_it = 5
            dtmax = 8
            sta_p, run, cp3, staff = SQO.sort(etas, t_it, t_tran, dtmax)
            step[c-1] = 2
            run_id = sorted(run, key = lambda x: int(x[1][2:]))
            # etas = []
        if step[c-1] == 3 or any(not sublist[6] for sublist in ac_info):#or step[c-1] == 4 or step[c-1] == 5:
            # p = 0
            for key, vals in staff.items():
                for i in range(len(vals)):
                    ind = int(vals[i][1][2:])
                    if i!= 0 and i!= len(vals)-1:
                        front = int(vals[i-1][1][2:])
                        back = int(vals[i+1][1][2:])
                        ac_info[ind] = [vals[i][1], plan_eta[ind], i, key, vals[i-1][1], plan_eta[front], (plan_eta[ind]-plan_eta[front])>=t_it-0.5, vals[i+1][1], plan_eta[back], (plan_eta[back]-plan_eta[ind])>=t_it]
                        if not ac_info[ind][6]:
                            buff[ind] = 1
                            re_sched.append([vals[i][1], plan_eta[ind], plan_eta[ind] + (t_it - (plan_eta[ind]-plan_eta[front])) + 2])
                        else:
                            buff[ind] = 2
                            
                    elif i == 0:
                        back = int(vals[i+1][1][2:])
                        ac_info[ind] = [vals[i][1], plan_eta[ind], i, key, vals[i][1], plan_eta[ind], True, vals[i+1][1], plan_eta[back], (plan_eta[back]-plan_eta[ind])>=t_it]
                        buff[ind] = 2
                    elif i == len(vals)-1:
                        front = int(vals[i-1][1][2:])
                        ac_info[ind] = [vals[i][1], plan_eta[ind], i, key, vals[i-1][1], plan_eta[front], (plan_eta[ind]-plan_eta[front])>=t_it-0.5, vals[i][1], plan_eta[ind], True]
                        if not ac_info[ind][6]:
                            buff[ind] = 1
                            re_sched.append([vals[i][1], plan_eta[ind], plan_eta[ind] + (t_it - (plan_eta[ind]-plan_eta[front])) + 2])
                        else:
                            buff[ind] = 2
            for i in ac_info:
                print(i[0], i[6], i[1]-i[5])
            print(buff)
            if step[c-1] == 3:
                step[c-1] = 4
            elif step[c-1] == 4:
                step[c-1] = 5
            else:
                step[c-1] = 6
        #     for i in range(len(run)):
        #         ind = run[i][1][2:]
        #         if i != 0:

        #             prev = 
        #             ac_info[ind] = [i, run[i][4], run[i-1], diffs[ind][1], ]
            
        # if i % 100 == 0:
        # if guide[0] == 4:
        
# for i in range(0,ntraf):
#     if entry[i][0]:
#         plt.plot(entry[i][1][0], entry[i][1][1], label = f'{bs.traf.id[i]} ENTRY PATH', linestyle = '--', linewidth = 2)
# for key, vals in sta_p.items():
#     for i in range(len(vals)):
#         ind = int(vals[i][1][2:])
#         if i!= 0 and i!= len(vals)-1:
#             front = int(vals[i-1][1][2:])
#             back = int(vals[i+1][1][2:])
#             ac_info[ind] = [vals[i][1], plan_eta[j], i, key, vals[i-1][1], plan_eta[front],(plan_eta[j]-plan_eta[front])>=5, vals[i+1][1], plan_eta[back], (plan_eta[back]-plan_eta[j])>=5]
#             if (plan_eta[j]-plan_eta[front])<5:
#                 buff[ind] = 1
#             else:
#                 buff[ind] = 2
#         elif i == 0:
#             back = int(vals[i+1][1][2:])
#             ac_info[ind] = [vals[i][1], plan_eta[j], i, key, vals[i][1], plan_eta[j], True, vals[i+1][1], plan_eta[back], (plan_eta[back]-plan_eta[j])>=5]
#             # if not ac_info[ind][6]:
#             #     buff[ind] = 1
#             # else:
#             buff[ind] = 2
#         elif i == len(vals)-1:
#             front = int(vals[i-1][1][2:])
#             ac_info[ind] = [vals[i][1], plan_eta[j], i, key, vals[i-1][1], plan_eta[front], (plan_eta[j]-plan_eta[front])>=5, vals[i][1], plan_eta[j], True]
#             if (plan_eta[j]-plan_eta[front])<5:
#                 buff[ind] = 1
#             else:
#                 buff[ind] = 2

print(ac_info)
print(run)
# print(sta_p)
print(final_order)
final_id = [i[0] for i in final_order]
print(final_id)
run_id = [i[1] for i in run]
print(run_id)
run_sort = sorted(run, key = lambda x: int(x[1][2:]))
print(run_sort[0])
print(f'AVERAGE COMPUTE TIME: {np.average(times)}, MAX COMPUTE TIME: {np.max(times)}')
print(f'AVERAGE ETA/TRAJECTORY TOA DIFFERENCE: {np.average(np.abs(diff))}, MAX DIFFERENCE: {np.max(np.abs(diff))}')
print(f'AMOUNT OF AIRCRAFT RESCHEDULES: {len(re_sched)}')
pct_error = []
for i in range(len(diff)):
    d = int(np.abs(diff[i]))
    t= int(run_sort[i][0])
    pct_error.append(100*(d/t))
    # pct_error.append(100*(np.abs(diffs[i])/run_sort[i][3]))
    if np.abs(diff[i]) > 3:
        print(bs.traf.id[i], bs.traf.lat[i], bs.traf.lon[i], bs.traf.hdg[i], run_sort[i][0], run_sort[i][3])
print(f'AVERAGE PERECENT ERROR: {np.average(pct_error)} MAX ERROR: {np.max(np.abs(pct_error))}')
# print(sorted(d2, key=lambda x: x[1]))
            
for i in range(len(bs.traf.id)):
    print(bs.traf.id[i], bs.traf.lat[i],',', bs.traf.lon[i], bs.traf.hdg[i])

end_time = time.time()
# print(all_dists)

masked_dists = np.where(all_dists == -1, np.inf, all_dists)  # Replace 0s with np.inf
valid_indices = np.where(np.logical_and(masked_dists != np.inf, np.arange(all_dists.shape[0])[:, None, None] < tag))

# Extract corresponding distances
valid_dists = masked_dists[valid_indices]

# Get indices of top 10 smallest values
sorted_indices = np.argsort(valid_dists)[:150]

# Extract the actual (i, j, k) indices from valid_indices
top_i = valid_indices[0][sorted_indices]
top_j = valid_indices[1][sorted_indices]
top_k = valid_indices[2][sorted_indices]
top_values = valid_dists[sorted_indices]

# Print results
for rank in range(len(top_values)):
    print(f"i: {top_i[rank]}, j: {top_j[rank]}, k: {top_k[rank]}, Distance: {top_values[rank]}")


print("SIM TIME", end_time-start_time)
print(sort_group1)

plt.figure()
'''Plotting'''
color_map = {}
for idx, acid in enumerate(bs.traf.id):

    if idx%2 == 0:
        marker = 'o'
    elif idx%3 == 0:
        marker = '^'
    else:
        marker = 'x'
    plt.text(res[0, 1, idx], res[0, 0, idx], f'AC{idx}', fontweight = 'bold', fontsize = 10)
    # plt.plot(res[::10, 1, idx], res[::10, 0, idx], marker = marker, linestyle = '')#, label=f'AX{idx}')
    # plt.plot([kcmh[1], res[0, 1, idx]], [kcmh[0], res[0, 0, idx]], linestyle = '--', color = 'black', linewidth = 2)
    line, = plt.plot(res[::10, 1, idx], res[::10, 0, idx], marker=marker, linestyle='')

    # Store the color assigned to this object
    color_map[idx] = line.get_color()

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
plt.scatter([favus_path[0][1], favus_path[1][1], favus_path[3][1], favus_path[4][1]],
            [favus_path[0][0], favus_path[1][0], favus_path[3][0], favus_path[4][0]], marker = 's', s = 100, color = 'brown')
plt.scatter([xavyr_path[0][1], xavyr_path[2][1], xavyr_path[3][1]],
            [xavyr_path[0][0], xavyr_path[2][0], xavyr_path[3][0]], marker = 's', s = 100, color = 'cyan')
plt.scatter([teeze_path[0][1], teeze_path[2][1], teeze_path[3][1], teeze_path[4][1], teeze_path[5][1]], 
            [teeze_path[0][0], teeze_path[2][0], teeze_path[3][0], teeze_path[4][0], teeze_path[5][0]], marker = 's', s = 100, color = 'orange')
plt.scatter([elupy_path[0][1], elupy_path[1][1], elupy_path[3][1], elupy_path[4][1], elupy_path[5][1]],
            [elupy_path[0][0], elupy_path[1][0], elupy_path[3][0], elupy_path[4][0], elupy_path[5][0]], marker = 's', s=  100, color = 'magenta')
plt.axis('equal')
# print(b_ll)
b2 = reformat_list(b_ll)
e2 = reformat_list(e_ll)
# plt.gca().set_aspect('equal', adjustable='datalim')
# c = 0
for i in b2:
    plt.plot(i[1], i[0], linewidth = 2, zorder = 100, color = 'black')
#     c+=1
# c = 0
for i in e2:
    if i != False:
        plt.plot(i[1], i[0], linestyle = '--', linewidth = 2, zorder = 100, color = 'black')
    # c+=1
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid()
plt.xlim(-83.3, -82.5)
plt.legend(loc = 'best')
# plt.show()
color_list = ['black', 'red', 'blue', 'green', 'yellow', 'orange', 'cyan', 'magenta', 'purple', 'brown']
c = 0
plt.figure()
for idx, acid in enumerate(bs.traf.id):
    # final_index = final_order.index(acid)
    # act_index = run.index(acid)
    try:
        final_index = next(i for i, item in enumerate(final_order) if item[0] == acid)
        act_index = next(i for i, item in enumerate(run) if item[1] == acid)  # Assuming run[i][1] is the ID
        # print(run[act_index])
        # print(final_order[final_index])
        if final_order[final_index][-1] in ['I', 'IV']:
            plt.scatter(1, run[act_index][0], marker = 's', color = color_map[idx], label = f'{acid} Delay {final_order[final_index][3]:0.2f}')
            if ac_info[idx][2]!=0:
                if ac_info[idx][1] - ac_info[idx][5] >= t_it-0.5:
                    plt.scatter(1.1, ac_info[idx][1], color = 'green')
                else:
                    plt.scatter(1.1, ac_info[idx][1], color = 'red')
            # plt.scatter(1.2, final_order[final_index][1], marker = 's', color = color_list[c])
            # plt.plot([1, 1.2], [run[act_index][0], final_order[final_index][1]], linestyle = '--', linewidth = 2, color = 'black')
            plt.plot([1, 1.1], [run[act_index][0], ac_info[idx][1]], linestyle = '--', linewidth = 2, color = 'black')
            plt.text(0.975, run[act_index][0], f'{acid}', weight = 'bold')
            # plt.text(1.25, final_order[final_index][1], f'{acid}', weight = 'bold')
        else:
            plt.scatter(1.3, run[act_index][0], marker = 's', color = color_map[idx], label = f'{acid} Delay {final_order[final_index][3]:0.2f}')
            if ac_info[idx][2]!=0:
                if ac_info[idx][1] - ac_info[idx][5] >= t_it-0.5:
                    plt.scatter(1.4, ac_info[idx][1], color = 'green')
                else:
                    plt.scatter(1.4, ac_info[idx][1], color ='red')

            # plt.scatter(1.6,final_order[final_index][1], marker = 's', color = color_list[c])
            plt.plot([1.3, 1.4], [run[act_index][0], ac_info[idx][1]], linestyle = '--', linewidth = 2, color = 'black')
            plt.text(1.25, run[act_index][0], f'{acid}', weight = 'bold')
            # plt.text(1.65, final_order[final_index][1], f'{acid}', weight = 'bold')
        # plt.text(2, final_order[act_index][1], f'Delay: {final_order[act_index][3]:0.2f}')
        c+=1
        if c == len(color_list):
            c = 0
    except:
        print(f'{acid} has not reached its final waypoints.')

for key, vals in staff.items():
    for i in range(len(vals)):
        ind = vals[i][1]
        if i!= 0:
            final_index = next(i for i, item in enumerate(final_order) if item[0] == ind)
            front = vals[i-1][1]
            front_index = next(i for i, item in  enumerate(final_order) if item[0] == front)
            act_index = next(i for i, item in enumerate(ac_info) if item[0] == ind)
            # back = int(vals[i+1][1][2:])
            # ac_info[ind] = [vals[i][1], plan_eta[ind], i, key, vals[i-1][1], plan_eta[front], (plan_eta[ind]-plan_eta[front])>=t_it-0.5, vals[i+1][1], plan_eta[back], (plan_eta[back]-plan_eta[ind])>=t_it]
            if key == 'R1':
                if final_order[final_index][1] - final_order[front_index][1]>=t_it-0.5:
                    plt.scatter(1.2, final_order[final_index][1], marker = 's', color = 'green')
                else:
                    plt.scatter(1.2, final_order[final_index][1], marker = 's', color = 'red')
                plt.plot([1.1, 1.2], [ac_info[act_index][1], final_order[final_index][1]], linestyle = '--', color = 'black', linewidth = 2)
            else:
                if final_order[final_index][1] - final_order[front_index][1]>=t_it-0.5:
                    plt.scatter(1.5, final_order[final_index][1], marker = 's', color = 'green')
                else:
                    plt.scatter(1.5, final_order[final_index][1], marker = 's', color = 'red')
                plt.plot([1.4, 1.5], [ac_info[act_index][1], final_order[final_index][1]], linestyle = '--', color = 'black', linewidth = 2)
        else:
            final_index = next(i for i, item in enumerate(final_order) if item[0] == ind)
            act_index = next(i for i, item in enumerate(ac_info) if item[0] == ind)
            if key == 'R1':
                plt.plot([1.1, 1.2], [ac_info[act_index][1], final_order[final_index][1]], linestyle = '--', color = 'black', linewidth = 2)
                plt.scatter(1.2, final_order[final_index][1], marker = 's', color = 'green')
            else: 
                plt.plot([1.3, 1.5], [ac_info[act_index][1], final_order[final_index][1]], linestyle = '--', color = 'black', linewidth = 2)
                plt.scatter(1.5, final_order[final_index][1], marker = 's', color = 'green')
                


plt.grid()
plt.xlim(0.9, 1.6)
plt.legend(loc = 'best')
plt.ylim(run[0][0]-25, run[-1][0]+25)
plt.title('R1                                                                                                                                        R2')
plt.text(0.98, run[0][0]-15, 'STA')
plt.text(1.18, run[0][0]-15, 'Actual')
# plt.xlabel(f'Scheduled Arrival Time Actual Arrival Time')
plt.ylabel('Arrival Time (s)')


plt.show()