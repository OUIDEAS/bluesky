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
import AircraftSpacingOptimize as MBO
import scipy.io
import pandas as pd
import json
import argparse
import scipy
import os

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

def convert_numpy(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(v) for v in obj]
    return obj

def to_json(obj, name, path):
    expanded_path = os.path.expanduser(path)

    os.makedirs(expanded_path, exist_ok=True)

    data_serializable = convert_numpy(obj)

    output_file = os.path.join(expanded_path, f"{name}.json")

    with open(output_file, "w") as f:
        json.dump(data_serializable, f, indent=4)

    print(f"{name} data saved to {output_file}")

def run_sim(scen, subscenario, spacing, t, expnum, exptype):
    '''
    Scenario 1 Options:
    BezAM -> Uses Bezier alternate maneuver to avoid emergency aircraft
    Hold -> Uses a holding pattern to avoid emergency aircraft

    '''

    scenario = scen
    subscen = subscenario

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
    # bs.stack.stack(f'CRELOG HoldSingleBypass 1')
    # bs.stack.stack(f'SCEN HoldSingleBypass; SAVEIC HoldSingleBypass')
    #Set Sim Time

    #Read from datasets
    wptdata, aptdata, awydata, firdata, codata, rwythresholds = loadnavdata.load_navdata()
    datasets_string = ['wptdata', 'aptdata']
    datasets = [wptdata, aptdata]

    #Make Aircraft
    home = [39.42, -82.2]
    des_dist = [0, -int(spacing)]
    ax1_ll = Meters_To_WSG84(des_dist, home)
    print('AX1 POSITION:', ax1_ll)

    mytraf = bs.traf.cre('AX0', 'M250', 39.42, -82.2, 0, 80, 57.412)
    mytraf = bs.traf.cre('AX1', 'M250', ax1_ll[0], -82.2, 0, 80, 57.412)
    # print(dist)
    # mytraf = bs.traf.cre('AC2', 'M250', 39.4144 , -82.2, 0, 80, 57.412)
    # mytraf = bs.traf.cre('AC3', 'M250', 39.4116, -82.2, 0, 80, 57.412)
    # mytraf = bs.traf.cre('AC4', 'M250', 39.4088 , -82.2, 0, 80, 57.412)

    bez_tot = []
    dub_tot = []

    for acid in bs.traf.id:
        bs.stack.stack(f'BANK {acid} 30')
        bs.stack.stack(f'CONFLICTDETECTION {acid} OFF')
    bs.stack.stack(f'ASAS OFF')

    # bs.stack.stack(f'MANUALPOSOrderedFlightAttempt add traf.id, traf.lat, traf.lon, traf.alt, traf.tas, traf.hdg')
    # bs.stack.stack(f'MANUALPOSOrderedFlightAttempt ON')
    dt = 0.01
    bs.stack.stack(f'DT {dt}')
    

    #bs.stack.stack(f'DT 0.025') #Set DT to .025 second

    # we'll run the simulation for up to 4000 seconds
    # h, toi_dist = qdrdist(39.42, -83.2, 39.422055, -83.2)
    # print('DIST: ', toi_dist)
    # bs.stack.stack(f'DIST 39.42, -83.2, 39.415, -83.2')
    utm_zone = get_utm_zone(-82.2)
    p = Proj(proj='utm',zone=utm_zone,ellps='WGS84')
    t_max = int(t)
    # t_max = 22500

    home = [ax1_ll[0], -82.2]
    des_dist = [0, -650]
    ex0_ll = Meters_To_WSG84(des_dist, home)
    print('EX0 START POSITION:', ex0_ll)

    if subscen == 'NotSingle':
        mytraf2 = bs.traf.cre('EX0', 'M250', ex0_ll[0], -82.2, 0, 80, 64)
    if subscen == 'Single':
        mytraf2 = bs.traf.cre('EX0', 'M250', 39.4186, -82.2, 0, 80, 64)
    cat = ['Fleet Aircraft', 'Fleet Aircraft', 'Emergency Vehicle']
    ac_data = []
    for i in range(len(bs.traf.id)-1):
        ac = {
            'ACID': bs.traf.id[i],
            'ExpNum': expnum,
            'Category': cat[i],
            'dt': dt,
            'ExpType': exptype,
            'ax_spacing': spacing
        }
        ac_data.append(ac)
    ac_df = pd.DataFrame(ac_data)


    ntraf = bs.traf.ntraf

    n_steps = t_max + 1
    t = np.linspace(0, t_max, n_steps)


    # allocate some empty arrays for the results
    res = np.zeros((n_steps, 5, ntraf))
    em = np.zeros((n_steps, 1))

    homexy = [0, 0]#, 0, 0, 0]
    homell = [0, 0]#, 0, 0, 0]
    req_head = [[0,0], [0,0]]#, [0,0], [0,0], [0,0]]
    bas = [[0, 0], [0, 0]]#, [0, 0], [0, 0], [0, 0]]
    guide = [0, 0]#, 0, 0, 0]
    index = [0, 0]#, 0, 0, 0]

    count = np.zeros(len(bs.traf.id))
    violation_count = [0, 0]#, 0, 0, 0]
    dist_to_EV = [[], []]#, [], [], []]


    ETA = [[0, 0],[0, 0]]#,[0, 0],[0, 0],[0, 0]]
    distances = np.zeros((len(bs.traf.id), n_steps))
    TOI = np.zeros((len(bs.traf.id)-1, n_steps))
    counter = [0, 0]#, 0, 0, 0]

    fleet_dist = np.zeros(n_steps)

    nodes = np.zeros((len(bs.traf.id)-1, 2), dtype = object)
    waypts = [0, 0]#, 0, 0, 0]
    if scenario == 'BezAM':
        start_end = [[0,0], [0,0]]#, [0,0], [0,0], [0,0]]
        times = [[0,0], [0,0]]#, [0,0], [0,0], [0,0]]
    else:
        start_end = [[[0,0],[0,0]], [[0,0],[0,0]]]#, [[0,0],[0,0]], [[0,0],[0,0]], [[0,0],[0,0]]]
        times = [[0,0,0], [0,0,0]]#, [0,0,0], [0,0,0], [0,0,0]]
    entry_point = [[0,0], [0,0]]#, [0,0], [0,0],[0,0]]

    # for i in range(0, 5):
    #     print(times[i][2] )
    events = []
    event_l = {
        'event': 'Simulation Start',
        'timeStamp': 0,
        'ACID': 'AX0, AX1, EX0',
        'ExpNum': expnum,
        'Category': 'All Aircraft',
        'ExpType': exptype,
        'ax_spacing': spacing
    }
    events.append(event_l)
    

    ev_stuff = []
    ev_entry = {
        'TOI': 9999,
        'TOI_Dist': 9999,
        'timeStamp': 0,
        'ax_spacing': spacing,
        'ACID': 'AX0, AX1, EX0',
        'ExpNum': expnum,
        'Category': 'All Aircraft',
        'ExpType': exptype
    }
    ev_stuff.append(ev_entry)
    # ev_specific = pd.DataFrame(ev_stuff)
    #Handy Conversions
    nmFeet = 6076.12
    msFeet = 3.28084
    nmMeters = 1852

    lr = 0

    sp_diff = 64-57.412
    fleet_rad = 8
    # print(nodes)
    # nodes = np.array([[[0], [0]], [[0],[0]]], dtype = object)
    # print(nodes[0][0])
    # nodes[0][0][0]= [1, 2, 23, 4, 53]
    # print(nodes[0][0])
    #The sim

    #home will be end point:
    home_end = [39.6, -82.2]
    k = 0
    k2 = [0, 0]#, 0, 0, 0]
    '''
    Initial ETAS
    '''
    for i in range(len(bs.traf.id)-1):
        goalHead, goalDist = qdrdist(bs.traf.lat[i], bs.traf.lon[i], home_end[0], home_end[1])
        print('Goal Distance:',goalDist)
        print('ID',bs.traf.id[i])
        ETA[i][0] = goalDist/57.412
        print('INITIAL ETA:', ETA[i][0], bs.traf.id[i])

    if scenario == 'BezAM':
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

                    TOI[j][i] = (toi_dist)/(sp_diff)
                    # print(TOI[j][i], bs.traf.id[j])
                    # print('AIRCRAFT:', bs.traf.id[j], 'DIST TO EV:', toi_dist, 'TOI:', TOI[j][i])
                    ev_entry = {
                        'TOI': TOI[j][i],
                        'TOI_Dist': toi_dist,
                        'timeStamp': i,
                        'ax_spacing': spacing,
                        'ACID': bs.traf.id[j],
                        'ExpNum': expnum,
                        'Category': 'Fleet Aircraft',
                        'ExpType': exptype
                    }
                    ev_stuff.append(ev_entry)
                    if j==1:
                        h, fd = qdrdist(bs.traf.lat[0], bs.traf.lon[0], bs.traf.lat[1], bs.traf.lon[1])
                        fleet_dist[i] = fd
                        # print(fd)
                        #calc distance
                        #append distance to list of dists
                        #if dist < radius (based on wingspan)
                        #add note event of aircraft being too close?
                        if fd <= fleet_rad:
                            note_event = {
                            'event': f'Fleet Aircraft Are Too Close: {fd}',
                            'timeStamp': i,
                            'ACID': 'AX0, AX1',
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                            events.append(note_event)
                            print(f'Fleet Aircraft Are Too Close: {fd} at timestamp {i}')

                    if j!=0 and bs.traf.lat[j] > bs.traf.lat[j-1] and guide[j] == 0:
                        # print('here!', bs.traf.id[j])
                        toi_dist = 999999
                        counter[j] = 1
                    if TOI[j][i] <= 10 and counter[j] == 0:
                        print(f'{bs.traf.id} GENERATING BEZIER PATH AT TIME STEP {i}')
                        note_event = {
                            'event': 'Bezier Generation',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # note_events = pd.DataFrame(events)
                        # print('TOIDIST', toi_dist)
                        k = i
                        counter[j] = 1
                        if j%2 == 0 or j == 0:
                            lr = 1
                            koz_x = 76
                        else:
                            lr = -1
                            koz_x = -76
                        print(lr)
                        target_toa1 = 1.5*((toi_dist+275)/bs.traf.tas[-1])
                        target_toa2 = 1.5*((toi_dist+275)/bs.traf.tas[-1])
                        print(target_toa1, target_toa2)
                        print(toi_dist, bs.traf.tas[-1])
                        # target_toa1 = 2*(275/64)
                        # target_toa2 = 2*(275/64)
                        homexy[j] = [0, 0]
                        homell[j] = [bs.traf.lat[j], bs.traf.lon[j]] #save as lonlat for WSG84
                        print('BEZ CURVE HOME',homell[j])
                        nodes1 = [np.array([homexy[j][0], lr*(221), lr*221]).flatten(), np.array([574.12, (574.12+275)/60+574.12, 711.62])]
                        nodes2 = [np.array([lr*221, lr*(225),homexy[j][0]]).flatten(), np.array([711.62, 750, 849.12])]

                        koz_bot = 589.12
                        koz_top = 834.12
                        corridor = lr*228
                        wpts, wpts_x, wpts_y, b1, b2, nodes1, nodes2, bezData = MBO.toCallOutside(bs.traf.tas[j], np.deg2rad(4.5), target_toa1, target_toa2,
                                                                bs.traf.hdg[j], nodes1, nodes2, koz_x, koz_bot, koz_top, lr,
                                                                [bs.traf.lat[j], bs.traf.lon[j]], corridor, i, bs.traf.id[j], expnum, exptype, spacing, homell[j])
                        bez_tot.append(bezData)
                        nodes[j][0] = nodes1
                        nodes[j][1] = nodes2
                        waypts[j] = [wpts_x, wpts_y]
                    if toi_dist <= 45.72:
                        violation_count[j]+=1
                        dist_to_EV[j].append(toi_dist)
                        # violation_count[j][1] = toi_dist
                        # print(toi_dist)
                    
                    if toi_dist <= 45.72 and counter[j] ==1:
                        print(f'{bs.traf.id} GENERATING INTERCEPTION PATHS AT TIME STEP: {i}')
                        note_event = {
                            'event': 'Dubins Generation',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # note_events = pd.DataFrame(events)
                        # print('TOIDIST', toi_dist)
                        k = i
                        counter[j] = 2
                        start_end[j][0] = [bs.traf.lat[j], bs.traf.lon[j]]
                        times[j][0] = i/100
                        pos = [bs.traf.lat[j], bs.traf.lon[j]]
                        new = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
                        # print('POSITION IS:', new, 'POSITION IN LATLON IS:', pos, 'HOME IS:', homell[j])
                        entry, exit, toa, wpts_x, wpts_y, dubData = MBO.EntryExitOutside(nodes1=nodes[j][0], nodes2=nodes[j][1], 
                                                                                         pos = new, velocity = bs.traf.tas[j], lr = lr, id = bs.traf.id[j],
                                                                                         timeStamp = i, expnum = expnum, exptype = exptype, spacing = spacing)
                        dub_tot.append(dubData)
                        end_point = [exit[0][0], exit[1][0]]
                        end_ll = Meters_To_WSG84(end_point, homell[j])
                        entry_point[j] = Meters_To_WSG84([entry[0][-1], entry[1][-1]], homell[j])
                        req_head[j][0] = entry[2]
                        req_head[j][1] = exit[2]
                        bas[j][0] = entry[-1][0]
                        bas[j][1] = exit[-1]
                        waypts[j] = [wpts_x, wpts_y]
                        if lr == -1:
                            toa = 21.57
                        # print('END POINT LL:', end_ll, 'END POINT XY:', end_point)
                        EVh, EV_dist = qdrdist(bs.traf.lat[-1], bs.traf.lon[-1], end_ll[0], end_ll[1])
                        print('EV TOA TO END:', EV_dist/64, EV_dist)
                        print(entry[-1][0])
                        if EV_dist/64 < toa:
                            note_event = {
                            'event': 'Path Flight Initiated',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                            events.append(note_event)
                            # note_events = pd.DataFrame(events)
                            print('PATH IS VALID AND SHOULD ALLOW FOR PASSAGE', 'TIME STEP:', i)
                            if lr == 1:
                                bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[-1][0]}; HDG {bs.traf.id[j]} {entry[2]}')
                                print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} AT {entry[-1][0]} DEGREES UNTIL A HEADING OF {entry[2]}')
                            else:
                                bs.stack.stack(f'BANK {bs.traf.id[j]} {entry[-1][0]}; HDG {bs.traf.id[j]} {90-entry[2]}')
                                print(f'INITIATING SINGLE BANK TURN FOR {bs.traf.id[j]} AT {entry[-1][0]} DEGREES UNTIL A HEADING OF {entry[2]}')
                        else:
                            print('GET FUCKED LOSER')
                    
                    # print(bs.traf.hdg[j])
                    
                    if lr == 1 and np.isclose(entry_point[j][0], bs.traf.lat[j], atol = 0.0001) and np.isclose(entry_point[j][1], bs.traf.lon[j], atol = 0.00001) and counter[j] == 2:
                        guide[j] = 1
                        k = i
                        note_event = {
                            'event': 'Initiating Pure Pursuit',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # note_events = pd.DataFrame(events)
                        bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                        print(f'INITIATING PURE PURSUIT GUIDANCE FOR {bs.traf.id[j]} AT TIME STEP: {i}')
                        counter[j] = 3
                    elif lr == -1 and np.isclose(entry_point[j][0], bs.traf.lat[j], atol = 0.001) and np.isclose(entry_point[j][1], bs.traf.lon[j], atol = 0.001) and counter[j] == 2:
                        guide[j] = 1
                        k = i
                        bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                        print(f'INITIATING PURE PURSUIT GUIDANCE FOR {bs.traf.id[j]} AT TIME STEP: {i}')
                        counter[j] = 3
                    

                    if guide[j] == 1 and counter[j] == 3:
                        pos = [bs.traf.lat[j], bs.traf.lon[j]]
                        m_pos = __WSG84_To_Meters_Single(waypoint =pos, home = homell[j], projobject=p)
                        # print(waypts[j][1])
                        # print(m_pos[1])
                        # print(waypts[j])
                        index[j], d = calc_Index(waypts[j], m_pos, index[j], bs.traf.id[j])
                        
                        head = np.arctan2(waypts[j][1][index[j]] - m_pos[1], waypts[j][0][index[j]]-m_pos[0])

                        bs.stack.stack(f'HDG {bs.traf.id[j]} {90-np.rad2deg(head)}')
                        # print(req_head[j][1])
                        # print(bs.traf.hdg[j], bs.traf.hdg[j]-360, 90-req_head[j][1], req_head[j][1], np.rad2deg(head))
                        if lr == 1:
                            if np.isclose(90-req_head[j][1], bs.traf.hdg[j]-360, atol = 1):
                                guide[j] = 2
                        else:
                            if np.isclose(90-req_head[j][1], bs.traf.hdg[j], atol = 5):
                                guide[j] = 2
                        
                    if guide[j] == 2:
                        note_event = {
                            'event': 'Nominal Path Interception Initiated',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # note_events = pd.DataFrame(events)
                        k = i
                        bs.stack.stack(f'BANK {bs.traf.id[j]} {exit[-1]}; HDG {bs.traf.id[j]} 0')
                        print(f'INITIATING NOMINAL PATH INTERCEPTION FOR {bs.traf.id[j]} AT {exit[-1]} AT TIME STEP: {i}')#
                        guide[j] = 3

                    if guide[j] == 3 and np.isclose(0, bs.traf.hdg[j], atol = 0.1):
                        k = i
                        # goalHead, goalDist = qdrdist(bs.traf.lat[j], bs.traf.lon[j], home_end[0], home_end[1])
                        # print(goalDist)
                        note_event = {
                            'event': 'Nominal Path Intercepted',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # note_events = pd.DataFrame(events)
                        print(f'{bs.traf.id[j]} IS BACK ON THE NOMINAL PATH AT TIME STEP {i}')
                        start_end[j][1] = [bs.traf.lat[j], bs.traf.lon[j]]
                        times[j][1] = i/100
                        # ETA[j][1] = goalDist/57.412
                        # print(f'ETA FOR {bs.traf.id[j]} BEFORE ALTERNATE MANEUVER IS {ETA[j][0]}. ETA AFTER ALTERNATE MANEUVER IS {ETA[j][1]}')
                        guide[j] = 0

                        
            # time.sleep(0.025)

    elif scenario == 'Hold':
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

                    TOI[j][i] = (toi_dist)/(sp_diff)
                    ev_entry = {
                        'TOI': TOI[j][i],
                        'TOI_Dist': toi_dist,
                        'timeStamp': i,
                        'ax_spacing': spacing,
                        'ACID': bs.traf.id[j],
                        'ExpNum': expnum,
                        'Category': 'Fleet Aircraft',
                        'ExpType': exptype
                    }
                    ev_stuff.append(ev_entry)
                    # print(TOI[j][i], bs.traf.id[j])
                    # print('AIRCRAFT:', bs.traf.id[j], 'DIST TO EV:', toi_dist, 'TOI:', TOI[j][i])
                    if j!=0 and bs.traf.lat[j] > bs.traf.lat[j-1] and guide[j] == 0:
                        # print('here!', bs.traf.id[j])
                        toi_dist = 999999
                        counter[j] = 1

                    if TOI[j][i] <= 10 and counter[j] == 0:
                        # print('TOIDIST', toi_dist)
                        counter[j] = 1
                        if j%2 == 0 or j == 0:
                            lr = 1
                        else:
                            lr = -1
                        print(lr)
                        note_event = {
                            'event': 'TOI is 10 or less',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                    if toi_dist <= 45.72:
                        violation_count[j]+=1
                        dist_to_EV[j].append(toi_dist)
                        # violation_count[j][1] = toi_dist
                        # print(toi_dist)
                    
                    if (toi_dist <= 45.72 or toi_dist == 999999) and counter[j] ==1: 
                        k = i
                        counter[j] =2
                        guide[j] = 1
                        note_event = {
                            'event': 'Holding Pattern Initiated',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                        }
                        events.append(note_event)
                        # print(bs.traf.lat)
                        start_end[j][0] = [bs.traf.lat[j], bs.traf.lon[j]]
                        times[j][0] = i/100
                        bs.stack.stack(f'BANK {bs.traf.id[j]} 73')
                        if lr == 1:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} 90; HDG {bs.traf.id[j]} 180')
                        else:
                            bs.stack.stack(f'HDG {bs.traf.id[j]} -90; HDG {bs.traf.id[j]} -180')
                        print(f'INITIATING HOLDING PATTERN FOR {bs.traf.id[j]} AT TIMESTAMP {i}')
                    
                    # if counter[j] == 2:
                    #     print(i/100 - times[j][0], bs.traf.hdg[j], bs.traf.id[j])
                        # print(bs.traf.hdg[j])
                    
                    if np.isclose(180, bs.traf.hdg[j], atol = 0.01) and counter[j] == 2:
                        # print(bs.traf.lat)
                        times[j][1] = i/100
                        counter[j] = 3

                    if subscen == 'Single':
                        if j!=0 and np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and counter[j] == 3 and bs.traf.lat[j]<bs.traf.lat[j-1]:
                            k = i
                            k2[j] = i/100 - times[j][1]
                            note_event = {
                            'event': 'Returning to Nominal Path',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                            events.append(note_event)
                            print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH AT TIMESTAMP {i}')
                            counter[j] = 4
                            if lr ==1:
                                bs.stack.stack(f'HDG {bs.traf.id[j]} 270; HDG {bs.traf.id[j]} 0')
                            else:
                                bs.stack.stack(f'HDG {bs.traf.id[j]} -270; HDG {bs.traf.id[j]} 0')
                        elif j ==0 and np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and i/100 - times[j][1] >= 5 and counter[j] == 3:
                            note_event = {
                            'event': 'Returning to Nominal Path',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                            events.append(note_event)
                            k = i
                            print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH AT TIMESTAMP {i}')
                            counter[j] = 4
                            if lr ==1:
                                bs.stack.stack(f'HDG {bs.traf.id[j]} 270; HDG {bs.traf.id[j]} 0')
                            else:
                                bs.stack.stack(f'HDG {bs.traf.id[j]} -270; HDG {bs.traf.id[j]} 0')
                    else:
                        if np.isclose(180.0, bs.traf.hdg[j], atol = 0.01) and i/100 - times[j][1] >= 5 and counter[j] == 3:
                            note_event = {
                            'event': 'Returning to Nominal Path',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                            events.append(note_event)
                            print(f'{bs.traf.id[j]} RETURNING TO NOMINAL PATH')
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
                        note_event = {
                            'event': 'Holding Pattern Ended',
                            'timeStamp': i,
                            'ACID': bs.traf.id[j],
                            'ExpNum': expnum,
                            'Category': 'Fleet Aircraft',
                            'ExpType': exptype,
                            'ax_spacing': spacing
                            }
                        events.append(note_event)
                        start_end[j][1] = [bs.traf.lat[j], bs.traf.lon[j]]
                        times[j][2] = i/100
                        print(f'HOLDING PATTERN FOR {bs.traf.id[j]} ENDED AT TIMESTAMP {i}')
                        counter[j] = 5
                        guide[j] = 0

                    

            # time.sleep(0.0025)


    end_time = time.time()
    bs.scr.update()
    print("SIM TIME", end_time-start_time)


    '''
    Computing ETA
    '''
    point_dist = [0, 0]# 0, 0, 0]
    nominal_time = [0, 0]#, 0, 0, 0]
    real_time = [0, 0]#, 0, 0, 0]
    delay = [0, 0]#, 0, 0, 0]


    # print('TIMES',times, times[0][2])
    # print('START_END',start_end)
    if scenario == 'BezAM':
        for i in range(0, 2):
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
        for i in range(0,2):
            delay[i] = times[i][2] - times[i][0] 
            real_time[i] = delay[i]
            nominal_time[i] = 0
            print(f'DATA FOR {bs.traf.id[i]}, {delay[i]}')
            ETA[i][1] = ETA[i][0]+delay[i]
            print(start_end[i])
            print(times[i])
    print('ETA FOR EACH AIRCRAFT:', ETA)
    print(start_end,times)

    delayData = {
        'AX0':{
            'delay': delay[0],
            'nom_time':nominal_time[0],
            'am_time': real_time[0],
            'point_dist': point_dist[0],
            'start_lat': start_end[0][0][0],
            'start_lon': start_end[0][0][1],
            'end_lat': start_end[0][1][0],
            'end_lon': start_end[0][1][1],
            'ACID': 'AX0',
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        },
        'AX1':{
            'delay': delay[1],
            'nom_time':nominal_time[1],
            'am_time': real_time[1],
            'point_dist': point_dist[1],
            'start_lat': start_end[1][0][0],
            'start_lon': start_end[1][0][1],
            'end_lat': start_end[1][1][0],
            'end_lon': start_end[1][1][1],
            'ACID': 'AX1',
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        }
    }
    # for i in range(len(violation_count)):
    #     violation_count[i]/=100
    # mins = []
    # for i in dist_to_EV:
    #     mins.append(np.argmin(i))

    # print('TOTAL TIME EACH AIRCRAFT VIOLATED THE RADIUS:', violation_count)
    # print('MINIMUM DISTANCE EACH AIRCRAFT GOT TO THE EV:', mins )
    '''Plotting'''
    # fig = plt.figure(1)
    # ax = fig.add_subplot(111)
    # # gpts = []
    # # for i in range(0, len(waypts[4][0])-1):
    # #     gpts.append([waypts[4][0][i], waypts[4][1][i]])
    # # # print(gpts)
    # # goal = Meters_To_WSG84(gpts, homell[4])
    # # goalx = []
    # # goaly = []
    # # for i in range(0, len(goal)-1):
    # #     goalx.append(goal[i][1])
    # #     goaly.append(goal[i][0])
    # # plt.scatter(goalx, goaly, color = 'green', marker = '*', s = 30, zorder = 50)

    # # gpts = []
    # # for i in range(0, len(waypts[3][0])-1):
    # #     gpts.append([waypts[3][0][i], waypts[3][1][i]])
    # # # print(gpts)
    # # goal = Meters_To_WSG84(gpts, homell[3])
    # # goalx = []
    # # goaly = []
    # # for i in range(0, len(goal)-1):
    # #     goalx.append(goal[i][1])
    # #     goaly.append(goal[i][0])
    # # plt.scatter(goalx, goaly, color = 'green', marker = '*', s = 30, zorder = 50)

    # # gpts = []
    # # for i in range(0, len(waypts[2][0])-1):
    # #     gpts.append([waypts[2][0][i], waypts[2][1][i]])
    # # # print(gpts)
    # # goal = Meters_To_WSG84(gpts, homell[2])
    # # goalx = []
    # # goaly = []
    # # for i in range(0, len(goal)-1):
    # #     goalx.append(goal[i][1])
    # #     goaly.append(goal[i][0])
    # # plt.scatter(goalx, goaly, color = 'green', marker = '*', s = 30, zorder = 50)

    # # gpts = []
    # # for i in range(0, len(waypts[1][0])-1):
    # #     gpts.append([waypts[1][0][i], waypts[1][1][i]])
    # # # print(gpts)
    # # goal = Meters_To_WSG84(gpts, homell[1])
    # # goalx = []
    # # goaly = []
    # # for i in range(0, len(goal)-1):
    # #     goalx.append(goal[i][1])
    # #     goaly.append(goal[i][0])
    # # plt.scatter(goalx, goaly, color = 'green', marker = '*', s = 30, zorder = 50)

    # # gpts = []
    # # j = 4
    # # for i in range(0, len(waypts[j][0])):
    # #     gpts.append([waypts[j][0][i], waypts[j][1][i]])
    # # # print(gpts)
    # # goal = Meters_To_WSG84(gpts, homell[j])
    # # goalx = []
    # # goaly = []
    # # for i in range(0, len(goal)):
    # #     goalx.append(goal[i][1])
    # #     goaly.append(goal[i][0])

    # # enpts = []
    # # for i in range(0, len(entry[0])):
    # #     enpts.append([entry[0][i], entry[1][i]])
    # # enp = Meters_To_WSG84(enpts, homell[j])
    # # enpx = []
    # # enpy = []
    # # for i in range(0, len(enp)):
    # #     enpx.append(enp[i][1])
    # #     enpy.append(enp[i][0])

    # # expts = []
    # # for i in range(0, len(exit[0])):
    # #     expts.append([exit[0][i], exit[1][i]])
    # # exp = Meters_To_WSG84(expts, homell[j])
    # # expx = []
    # # expy = []
    # # for i in range(0, len(exp)):
    # #     expx.append(exp[i][1])
    # #     expy.append(exp[i][0])    
    # # expx.append(goalx[-1])
    # # expy.append(goaly[-1])
    # # ax.plot(goalx, goaly, color = 'magenta', label = 'Partial Bezier Path')
    # # ax.plot(enpx, enpy,  color = 'cyan', label = 'Interception Path')
    # # ax.plot(expx, expy,  color = 'cyan')

    # # print(goal[0])
    # # print(ETA)
    # # print(res[4])
    # # print(res[5])
    # # markers = ['o', 's', '^', 'v', '>', '<', 'P', 'X', 'D', '*', '+', 'x', '|', '_', '1', '2', '3', '4', 'h', 'H']
    # used_colors = []
    # # for idx, acid in enumerate(bs.traf.id):
    # #     available_colors = [color for color in range(1, 101) if color not in used_colors]
    # #     color_idx = np.random.choice(available_colors)
    # #     color = plt.cm.tab20(color_idx)  # Use a colormap to generate a color
    # #     used_colors.append(color_idx)  # Add the used color index to the list
    # #     if idx%2 == 0:
    # #         marker = 'o'
    # #     elif idx%3 == 0:
    # #         marker = '^'
    # #     else:
    # #         marker = 'x'
    # #     # marker = np.random.choice()  # Randomly select a marker style
    # #     color = np.random.rand(3,)  # Randomly select a color
    # #     if acid == 'AC4' or acid == 'EM0':
    # #         print(idx)
    # #         plt.plot(res[::25, 1, idx], res[::25, 0, idx], marker = marker, label=f'{acid}', color=color
    # # bez1_ll = Meters_To_WSG84(waypts[4], homell[4])

    # # for idx, acid in enumerate(bs.traf.id):
    # #     available_colors = [color for color in range(1, 101) if color not in used_colors]
    # #     color_idx = np.random.choice(available_colors)
    # #     color = plt.cm.tab20(color_idx)
    # #     used_colors.append(color_idx)
    # #     # plt.scatter(bez1_ll[1], bez1_ll[0])
    # #     marker = 'o' if idx % 2 == 0 else '^' if idx % 3 == 0 else 'x' if idx == 1 else '*'
    # #     color = np.random.rand(3,)


    # # if acid == 'AC4' or acid == 'EM0':

    # plt.plot(res[::10, 1, -1], res[::10, 0, -1], marker='x', label='Emergency Vehicle', color='red', zorder = 1000)
    # plt.plot(res[::10, 1, 0], res[::10, 0, 0], marker='o', label='Aircraft 0', color='blue', zorder = 10)
    # plt.plot(res[::10, 1, 1], res[::10, 0, 1], marker='o', label='Aircraft 1', color='green', zorder = 10)
    # # plt.plot(res[::10, 1, 2], res[::10, 0, 2], marker='o', label='Aircraft 2', color='orange', zorder = 10)
    # # plt.plot(res[::10, 1, 3], res[::10, 0, 3], marker='o', label='Aircraft 3', color='purple', zorder = 10)
    # # plt.plot(res[::10, 1, 4], res[::10, 0, 4], marker='o', label='Aircraft 4', color='brown', zorder = 10)


    # h, toi_dist = qdrdist(bs.traf.lat[0], bs.traf.lon[0], bs.traf.lat[-1], bs.traf.lon[-1])

    # # ax.text(bs.traf.lon[-1]-0.005, bs.traf.lat[-1]+0.001, f'Emergency Vehicle TOA To Goal {ev_TOA:.3g}', fontsize = 7)
    # # ax.text(bs.traf.lon[0]+0.0005, bs.traf.lat[0]+0.0005, f'AC0 Resumes Flight On The Nominal Path At t = {k/100}', fontsize = 10)
    # # ax.plot([bs.traf.lon[0]+0.001, bs.traf.lon[0]], [bs.traf.lat[0], bs.traf.lat[0]], linestyle = 'dashed', color = 'black', linewidth = 2)
    # # ax.plot([bs.traf.lon[0]+0.001, bs.traf.lon[0]], [bs.traf.lat[-1], bs.traf.lat[-1]], linestyle = 'dashed', color = 'black', linewidth = 2)
    # # ax.text(bs.traf.lon[0]+0.0001, (bs.traf.lat[0]+bs.traf.lat[-1])/2, f'Time Separation: {TOI[0][-1]:.3g}s, Distance: {toi_dist:.3g}m', fontsize = 11)
    # # ax.text(bs.traf.lon[2]-0.0025, bs.traf.lat[2]-0.001, f'AC2 Begins Return To Nominal Path At t = {k/100}, After AC1 Has Passed It', fontsize = 10)
    # # ax.text(bs.traf.lon[3]-0.0025, bs.traf.lat[3]+0.001, f'AC3 Begins Holding Pattern At t = 33.81, After It Has Passed AC2', fontsize = 10)

    # # # ax.scatter()\size = 7)
    # # # ax.scatter()\
    # # # ax.plot([bs.traf.lon[0]-0.001, bs.traf.lon[1]+0.001], [bs.traf.lat[0], bs.traf.lat[1]], color = 'blue', linestyle = 'dashed', linewidth = 2)
    # # ax.plot([bs.traf.lon[0]-0.005, bs.traf.lon[0]], [39.4287, 39.4287], color = 'black', linestyle = 'dashed', linewidth = 2)
    # # ax.plot([bs.traf.lon[0]-0.005, bs.traf.lon[0]], [bs.traf.lat[1], bs.traf.lat[1]], color = 'black', linestyle = 'dashed', linewidth = 2)
    # # ax.text(bs.traf.lon[0]-0.005, (bs.traf.lat[2]+39.4287)/2, f'{k2[2]:.3g} Second Travel Time', fontsize = 10)
    # # ax.text(bs.traf.lon[1]+0.0005, bs.traf.lat[1], f'AC1 Is Passing AC0', fontsize = 10)
    # # ax.text(bs.traf.lon[1]+0.0005, bs.traf.lat[1]-0.0005, f'AC1 Initiates Holding Pattern At t = {k/100}', fontsize = 10)
    # plt.axis('equal')
    # # shape4PlotAirport('COJEZ', 'NIKOE')
    # # plotWPTCoords(wpts)
    # plt.title(f'Flight Snapshot At t = {t_max/100}')
    # # shape4PlotAirport('COJEZ', 'NIKOE')
    # # plotWPTCoords(wpts)
    # plt.xlabel('Longitude')
    # plt.ylabel('Latitude')
    # plt.grid()
    # plt.legend()

    # plt.show()
    print('FINAL TOI:', TOI[1][-1])

    # data = {
    #     'AC0': np.array([res[::10, 1, 0], res[::10, 0, 0], res[::10, 3, 0], res[::10, 4, 0]]),
    #     'AC1': np.array([res[::10, 1, 1], res[::10, 0, 1], res[::10, 3, 1], res[::10, 4, 1]]),
    #     'AC2': np.array([res[::10, 1, 2], res[::10, 0, 2], res[::10, 3, 2], res[::10, 4, 2]]),
    #     'AC3': np.array([res[::10, 1, 3], res[::10, 0, 3], res[::10, 3, 3], res[::10, 4, 3]]),
    #     'AC4': np.array([res[::10, 1, 4], res[::10, 0, 4], res[::10, 3, 4], res[::10, 4, 4]]),
    #     'EM0': np.array([res[::10, 1, 5], res[::10, 0, 5], res[::10, 3, 5], res[::10, 4, 5]]),
    #     # 'Waypts': np.array([goalx, goaly]),
    #     # 'Entry': np.array([enpx, enpy]),
    #     # 'Exit': np.array([expx, expy])
    # }
    # scipy.io.savemat('Scen1Data2Hold.mat', data)
    # print(res[::1, 0, 0], res[::1, 1, 1]) 
    # #lat, lon, alt, tas, hdg
    stateData = {
        'AC0':{
            'lat': np.array([res[::1, 0, 0]]),
            'lon': np.array([res[::1, 1, 0]]),
            'Heading': np.array([res[::1, 4, 0]]),
            'V_meters': np.array([res[::1, 3, 0]]),
            'V_knots': np.array([res[::1, 3, 0]])*1.9438445,
            'timeStamp': np.array(range(0, n_steps)),
            'ACID': 'AX0',
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        },
        'AC1':{
            'lat': np.array([res[::1, 0, 1]]),
            'lon': np.array([res[::1, 1, 1]]),
            'Heading': np.array([res[::1, 4, 1]]),
            'V_meters': np.array([res[::1, 3, 1]]),
            'V_knots': np.array([res[::1, 3, 1]])*1.9438445,
            'timeStamp': np.array(range(0, n_steps)),
            'ACID': 'AX1',
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        },
        'EV0':{
            'lat': np.array([res[::1, 0, 2]]),
            'lon': np.array([res[::1, 1, 2]]),
            'Heading': np.array([res[::1, 4, 2]]),
            'V_meters': np.array([res[::1, 3, 2]]),
            'V_knots': np.array([res[::1, 3, 2]])*1.9438445,
            'timeStamp': np.array(range(0, n_steps)),
            'ACID': 'EX0',
            'ExpNum': expnum,
            'Category': 'Emergency Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        }

    }

    # note_events = pd.DataFrame(events)
    # output_file = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\NotableEvents_{spacing}_Apart.json'
    # note_events.to_json(output_file, orient='records', indent=4)
    # print(f"Note Events data saved to {output_file}")

    # ev_specific = pd.DataFrame(ev_stuff)
    # output_file = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\EVSpecific_{spacing}_Apart.json'
    # ev_specific.to_json(output_file, orient='records', indent=4)
    # print(f"EVSPecific data saved to {output_file}")

    # output_file = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Aircraft_{spacing}_Apart.json'
    # ac_df.to_json(output_file, orient='records', indent=4)
    # print(f"Aircraft data saved to {output_file}")


    # path = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs'
    # to_json(stateData, f'State_{spacing}_Apart', path)
    # to_json(delayData, f'Delay_{spacing}_Apart', path)
    # if 'Bez' in exptype or 'Bez' in scen:
    #     bez = pd.DataFrame(bez_tot)
    #     output_file = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Bez_{spacing}_Apart.json'
    #     bez.to_json(output_file, orient='records', indent=4)
    #     print(f"Bez data saved to {output_file}")

    #     dub = pd.DataFrame(dub_tot)
    #     output_file = f'C:\\Users\\Michael\\Desktop\\BlueSkyData\\{exptype}JSONs\\Dubins_{spacing}_Apart.json'
    #     dub.to_json(output_file, orient='records', indent=4)
    #     print(f"Dubins data saved to {output_file}")
        # to_json(bezData, f'Bez_{spacing}_Apart', path)
        # to_json(dubData, f'Dubins_{spacing}_Apart', path)
    note_events = pd.DataFrame(events)
    output_file = f'~/bluesky/BlueSkyData/{exptype}JSONs/NotableEvents_{spacing}_Apart.json'
    note_events.to_json(output_file, orient='records', indent=4)
    print(f"Note Events data saved to {output_file}")

    ev_specific = pd.DataFrame(ev_stuff)
    output_file = f'~/bluesky/BlueSkyData/{exptype}JSONs/EVSpecific_{spacing}_Apart.json'
    ev_specific.to_json(output_file, orient='records', indent=4)
    print(f"EVSPecific data saved to {output_file}")

    output_file = f'~/bluesky/BlueSkyData/{exptype}JSONs/Aircraft_{spacing}_Apart.json'
    ac_df.to_json(output_file, orient='records', indent=4)
    print(f"Aircraft data saved to {output_file}")


    path = f'~/bluesky/BlueSkyData/{exptype}JSONs'
    to_json(stateData, f'State_{spacing}_Apart', path)
    to_json(delayData, f'Delay_{spacing}_Apart', path)
    # print(exptype, scen)
    if 'Bez' in exptype or 'Bez' in scen:
        bez = pd.DataFrame(bez_tot)
        output_file = f'~/bluesky/BlueSkyData/{exptype}JSONs/Bez_{spacing}_Apart.json'
        bez.to_json(output_file, orient='records', indent=4)
        print(f"Bez data saved to {output_file}")

        dub = pd.DataFrame(dub_tot)
        output_file = f'~/bluesky/BlueSkyData/{exptype}JSONs/Dubins_{spacing}_Apart.json'
        dub.to_json(output_file, orient='records', indent=4)
        print(f"Dubins data saved to {output_file}")
    #     to_json(bezData, f'Bez_{spacing}_Apart', path)
    #     to_json(dubData, f'Dubins_{spacing}_Apart', path)

    
    
    

if __name__ == '__main__':
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
    import AircraftSpacingOptimize as MBO
    import scipy.io
    import pandas as pd
    import json
    import argparse
    import scipy
    windows = False
    if not windows:
        parser = argparse.ArgumentParser(description='Apply different spacing between fleet aircraft')
        parser.add_argument('-s1', '--scenario')
        parser.add_argument('-s2', '--subscenario')
        parser.add_argument('-d', '--spacing')
        parser.add_argument('-t', '--time')
        parser.add_argument('-en', '--expnum')
        parser.add_argument('-et', '--exptype')
        args = parser.parse_args()
        run_sim(args.scenario, args.subscenario, args.spacing , args.time, args.expnum, args.exptype)
    else:
        run_sim('BezAM', "NotSingle", 311 , 17000, 0, "SuperTest")