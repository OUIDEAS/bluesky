import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy.io
import sympy as sp
from scipy.optimize import brentq
from scipy.optimize import fsolve
from pyproj import Proj
from utm import from_latlon as llutm
from utm import to_latlon as utmll
import pandas as pd

def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

def manual_bez_xy(P0, P1, P2, points):
    t = np.linspace(0, 1, points)
    x = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    y = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    return np.array([x, y])

def find_bez_xy(P0, P1, P2, t):
    x = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    y = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    return np.array([x, y])

def manual_bez_partial(P0, P1, P2, points, t_start):
    t = np.linspace(t_start,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

def manual_bez_partialExit(P0, P1, P2, points, t_end):
    t = np.linspace(0,t_end,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

def curvature(P0,P1,P2):
    mx = np.mean([P0[0], P2[0]])
    my = np.mean([P0[1], P2[1]])
    m = [mx, my]
    center1x = np.mean([mx, P0[0]])
    center1y = np.mean([my, P0[1]])
    r1 = np.sqrt(
        (center1x-mx)**2 + (center1y-my)**2
    )
 
    center2x = np.mean([mx, P2[0]])
    center2y = np.mean([my, P2[1]])
    r2 = np.sqrt(
        (center2x-mx)**2 + (center2y-my)**2
    )
    circle1 = [center1x, center1y, r1]
    circle2 = [center2x, center2y, r2]
    p1_from_c1 = np.sqrt( (circle1[0] - P1[0])**2 + (circle1[1] - P1[1])**2)
    p1_from_c2 = np.sqrt( (circle2[0] - P1[0])**2 + (circle2[1] - P1[1])**2)
    area = np.abs(
        P0[0]*P1[1] + P1[0]*P2[1] + P2[0]*P0[1] - P0[1]*P1[0] - P1[1]*P2[0] - P2[0]*P0[0]
    )/2
    if p1_from_c1 <= circle1[2]:
        # print("IN C1")
        kmax = area / (np.sqrt((P0[0]-P1[0])**2+(P0[1]-P1[1])**2))**3
 
    elif p1_from_c2 <= circle2[2]:
        # print("IN C2")
        kmax = area / (np.sqrt((P1[0]-P2[0])**2+(P1[1]-P2[1])**2))**3
 
    else:
        # print("NOT IN C1 or C2")
        kmax = (np.sqrt((m[0]-P1[0])**2+(m[1]-P1[1])**2))**3 / (area**2)
 
    roc = 1/kmax
    return roc
 
def path_length(P1, P0, P2, t):
    ax = P0[0] - 2*P1[0] + P2[0]
    ay = P0[1] - 2*P1[1] + P2[1]
    bx = 2*P1[0] - 2*P0[0]
    by = 2*P1[1] - 2*P0[1]
    A = 4 * (ax**2 + ay**2)
    B = 4 * (ax*bx + ay*by)
    C = bx**2 + by **2
    b=B/(2.0*A)
    c=C/A
    u=t+b
    k=c-(b*b)
 
    L=0.5*np.sqrt(A)*((u*np.sqrt((u*u)+k))
        -(b*np.sqrt((b*b)+k))
        +(k*np.log(np.abs((u+np.sqrt((u*u)+k))/(b+np.sqrt((b*b)+k)))))
    )
    return L

 
def solve_optim1(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr, line):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    cons = (
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
            {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
            {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
            ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
 
    return val.x , curvature(P0,val.x,P2)

def solve_optim2(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr, line):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))

    cons = (
            
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
            # {'type': 'ineq', 'fun': lambda x: x[1] - P0[1]},
            {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
        ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

    return val.x , curvature(P0,val.x,P2)

def solve_optimEntry(guess, max_bank,  min_bank, nodes, velocity, pos, lr):
    def path_cost(ba, nodes, velocity, lr):
        diff = find_diff_entry(ba, nodes, velocity, pos, lr)

        return np.abs(diff[0])
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}  
    )
    val= minimize(path_cost,guess,(nodes, velocity, lr), method='SLSQP', tol=1E-10, constraints=cons)
    t_val = find_diff_entry(val.x, nodes, velocity, pos, lr)[1]
    return val.x, t_val

def solve_optimExit(guess, max_bank, min_bank, min_t, nodes, velocity, lr):
    #t is also being minimized, so part of constraints
    def path_cost(guess, nodes, velocity, lr):
        diff = find_diff_exit(guess, nodes, velocity, lr)
        # print('DIFF FOR EXIT:', diff[0])
        return np.abs(diff[0]) 

    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank},
            {'type': 'ineq', 'fun': lambda x: x[1] - min_t},
            {'type': 'ineq', 'fun': lambda x: 1-x[1]},
    )
    val =  minimize(path_cost, guess, (nodes, velocity, lr), method = 'SLSQP', tol = 1E-10, constraints=cons)
    return val.x

def central_angle(center, point1, point2):
    # Calculate vectors
    v1 = (point1[0] - center[0], point1[1] - center[1])
    v2 = (point2[0] - center[0], point2[1] - center[1])
    print(v1,v2)
    # Compute the angle using atan2
    angle_radians = np.arctan2(v2[1]-v1[1], v2[0]-v1[0])# - math.atan2(v1[1], v1[0])
    print(angle_radians)
    # Normalize the angle to [0, 2*pi] if necessary
    angle_radians = (angle_radians + 2 * math.pi) % (2 * math.pi)
    
    # Convert to degrees (optional)
    angle_degrees = np.rad2deg(angle_radians)
    print('CA:', angle_radians, angle_degrees)
    
    return angle_radians, angle_degrees

def find_diff_entry(ba, nodes, velocity, pos, lr):
    pi = np.pi
    t_guess = 0.5
    mindiff = 50
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    tr = velocity**2 / (11.26*math.tan(ba[0]))
    tr*=0.3048
    h = lr*tr
    k = pos[1]
    circle_eq = lambda t: (Bx(t)-h)**2+(By(t)-k)**2 - tr**2

    S = fsolve(circle_eq, t_guess)
    # print(S)
    t_final = 0

    for i in S:

        if i>=0 and i <= 1:
            index = int(np.round(float(i*200)))
            if index == 200:
                index = 199
            bez_angle = np.arctan2(By(i+.01)-By(i), Bx(i+0.01)-Bx(i))

            x_l = [i for i in np.linspace(0, Bx(i), 200)] #Space between nominal path and bez
            y = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            # plt.plot(path[0], path[1])
            # plt.plot(x_l, y)
            # plt.axis('equal')
            # plt.show()

            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            diff = np.abs(bez_angle-int_angle)
            if diff < mindiff and np.abs(y[-1] - By(i))<=.025 and np.abs(x_l[-1]-Bx(i))<=0.025:
                mindiff = diff
                # plt.plot(path[0], path[1])
                # plt.plot(x_l, y)
                # plt.scatter(path[0][index], path[1][index], color = 'yellow', marker = '*', s = 100)
                # plt.scatter(x_l[-1], y[-1], color = 'purple', marker = '*', s = 100)
                # plt.axis('equal')
                # plt.show()
                t_final = i
            else:
                t_final = 0
    # print(mindiff, t_final)
    return [mindiff, t_final]

def find_diff_exit(guess, nodes, velocity, lr):
    ba = guess[0]
    t_guess = guess[1]
    # print(guess, ba, t_guess)
    mindiff = 100
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    bezHead = np.arctan2(By(t_guess)-By(t_guess-0.01), Bx(t_guess)-Bx(t_guess-0.01))

    tr = (velocity*1.94384)**2 / (11.26*math.tan(ba))
    tr*=0.3048
    h = lambda t: Bx(t) - tr*math.cos(bezHead)
    k = lambda t: By(t) + tr*math.sin(bezHead)

    circle_eq = lambda y: (0-h(t_guess))**2+(y-k(t_guess))**2 - tr**2
    if lr == -1:
        y_guess = 1200
    else:
        y_guess = 1400
    S = fsolve(circle_eq, y_guess) #gives y intersection
    t_final = .5

    y_l = [i for i in np.linspace(274, S, 200)]
    x = [0 for i in y_l]

    # print(S)
    for i in S:
        if i > 0: 
            y = [k for k in np.linspace(S, By(t_guess))]
            if lr == 1:
                x_l = [h(t_guess) - np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
                # print(x_l)
            else:
                x_l = [h(t_guess) + np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            int_angle = np.arctan2(y[0]-y[1], x_l[0]-x_l[1])
            if x_l[0] == 0 and y[0] > By(t_guess):# and np.isclose(int_angle, math.pi/2, atol = 2.5):# <= 0.01:
                # x_l = [i for i in np.linspace(750, Bx(t_guess))]
                # y = [k(t_guess)-np.sqrt(tr**2 - (x-h(t_guess))**2) for x in x_l]
                int_angle = np.arctan2(y[0]-y[1], x_l[0]-x_l[1])
                diff = np.abs((np.pi/2) - int_angle)
                guess_deg = np.rad2deg(ba)
                # if lr == 1:
                    # plt.title(f'{guess_deg} {t_guess}')
                    # plt.plot(x, y_l, linestyle = 'dashed')
                    # plt.plot(path[0], path[1])
                    # plt.plot(x_l, y)
                    # plt.axis('equal')
                    # plt.show()
                if diff < mindiff:
                    mindiff = diff
                    # print('BLARGY', diff)
                    t_final = t_guess
                

    return [mindiff, t_final]

def bez_to_wp(bez1, bez2, num):
    path_wpts = []
    point_marker = len(bez1[0])/num
    point_marker = np.round(point_marker)
    # for i in range(len(bez3[0])):
    #     if i%point_marker == 0:
    #         path_wpts.append([bez3[0][i], bez3[1][i]]) 
    for i in range(len(bez1[0])):
        if i%point_marker == 0:
            path_wpts.append([bez1[0][i], bez1[1][i]]) 
    for i in range(len(bez2[0])):
        if i%point_marker == 0:
            path_wpts.append([bez2[0][i], bez2[1][i]]) 
    wpts_x = [x[0] for x in path_wpts]
    wpts_y = [y[1] for y in path_wpts]
    return path_wpts, wpts_x, wpts_y

def paths_to_wp(paths, num):
    path_wpts = []
    num_paths = len(paths)
    for i in range(0, num_paths):
        path = paths[i]
        marker = len(path[0])/num
        # print(marker)
        for j in range(0, len(path[0])):
        #     print(j)
            if j%marker == 0:
                path_wpts.append([path[0][j], path[1][j]])
    wpts_x = [x[0] for x in path_wpts]
    wpts_y = [y[1] for y in path_wpts]
    return path_wpts, wpts_x, wpts_y

def bez_to_wp_single(bez1, num):
    path_wpts = []
    point_marker = len(bez1[0])/num
    point_marker = np.round(point_marker)
    for i in range(len(bez1[0])):
        if i%1 == 0:
            path_wpts.append([bez1[0][i], bez1[1][i]]) #building list of waypoints
    wpts_x = [x[0] for x in path_wpts]
    wpts_y = [y[1] for y in path_wpts]
    return path_wpts, wpts_x, wpts_y

def validity_check(x, y1, y2, points, corridor, curve, lr):
    # points = np.array(points)
    if lr == 1:
        for i in range(len(points)):
            # print('P',points[i][0],'X',x,'Y1',y1,'Y2',y2)
            if points[i][0]<=x and points[i][1]>=y1 and points[i][1]<=y2 or points[i][0] >= corridor:
                # print('In Positive')
                # print(curve)
                # print(points[i][0], x, points[i][0]<=x)
                # print(points[i][1], y1, points[i][1]>=y1)
                # print(points[i][1], y2, points[i][1]<=y2)
                # print(points[i][0], corridor, points[i][0] >= corridor)
                # plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    elif lr == -1:
        for i in range(len(points)):
            if points[i][0]>=x and points[i][1]>=y1 and points[i][1]<=y2 or points[i][0] <= corridor:
                # print('In Negative')
                # print(curve)
                # print(points[i][0], x, points[i][0]>=x)
                # print(points[i][1], y1, points[i][1]>=y1)
                # print(points[i][1], y2, points[i][1]<=y2)
                # print(points[i][0], corridor, points[i][0] <= corridor)
                # plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    return True

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

def __WSG84_To_Meters_Single(waypoint, home):
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    homeX, homeY = p(home[1], home[0])

    lon = waypoint[0]
    lat = waypoint[1]
    x,y = p(lat, lon)
    x = (x - homeX)
    y = (y - homeY)
        
    return [x, y]

def to_UTM(lat, lon):
    utm = llutm(lat, lon)
    x = utm[0]
    y = utm[1]
    return [x, y]

def to_LL(east, north):
    ll = utmll(east, north, 17, 'S')
    lat = ll[0]
    lon = ll[1]
    return [lat, lon]
    
def toCallOutside(velocity, turn_rate, target_toa1, target_toa2, uav_head, nodes1, nodes2, koz_x, koz_bot, koz_top, lr, latlon, corridor, timeStamp, acid, expnum, exptype, spacing, homell):
    # velocity  = 188 #ft/s
    start_time = time.time()
    # turn_rate = np.deg2rad(4.50) # RAD/s
    # turn_radius = velocity / turn_rate
    turn_radius = (velocity*1.94384)**2/(11.26*math.tan(np.deg2rad(73)))
    turn_radius*=0.3048
    valid1 = False
    valid2 = False
    orig_toa1 = target_toa1
    orig_toa2 = target_toa2
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    #LON = NODES[0]
    #LAT = NODES[1]

# UTM CONVERSIONS

    # nodes1= to_UTM(nodes1[1], nodes1[0]) # gives [East, North] Should index the same as before
    # nodes2 = to_UTM(nodes2[1], nodes2[0])
    # koz_bot = to_UTM(koz_bot[0], koz_bot[1])
    # koz_top = to_UTM(koz_top[0], koz_top[1])
    # koz_x = to_UTM(koz_x[0], koz_x[1])
    # corridor = to_UTM(corridor[0], corridor[1])
    # print('IMPORTANT', corridor[0], nodes1[0][1])
    # plt.scatter(nodes1[0], nodes1[1])
    # plt.scatter(nodes2[0], nodes2[1])
    # plt.plot([koz_x, koz_x], [koz_bot, koz_top], color = 'red', linestyle = '--')
    # plt.plot([nodes1[0][0], koz_x], [koz_bot, koz_bot], color = 'red', linestyle = '--')
    # plt.plot([nodes1[0][0], koz_x], [koz_top, koz_top], color = 'red', linestyle = '--')
    # plt.plot([corridor, corridor], [nodes1[1][0], nodes2[1][2]], color = 'orange', linestyle = 'dashed')
    # plt.show()

    
    m_p21 = (nodes2[1][1]-nodes1[1][2])/(nodes2[0][1]-nodes1[0][2])
    b2 = nodes2[1][1] - m_p21*nodes2[0][1]

    lines_coeffs2 = [m_p21, b2]
    
    c = 0
    while not valid1 or not valid2:
        optim_sol1, curv1 = solve_optim1(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=target_toa1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs2)
        
        m_p12 = (nodes1[1][2]-optim_sol1[1])/(nodes1[0][2]-optim_sol1[0])


        b = optim_sol1[1] - m_p12*optim_sol1[0]
        lines_coeffs = [m_p12, b]
        # print('M12',m_p12, 'B', b)

        x_slope = [i for i in np.linspace(optim_sol1[0]-10, lr*245, 50)]
        y_slope = [m_p12*x + b for x in x_slope]

        # print(lines_coeffs)

        optim_sol2, curv2 = solve_optim2(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=target_toa2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs)
        m_p21 = (nodes1[1][2]-optim_sol2[1])/(nodes1[0][2]-optim_sol2[0])
        b2 = optim_sol2[1] - m_p21*optim_sol2[0]
        lines_coeffs2 = [m_p21, b2]
        # if mark == True:

        # x_slope2 = [i for i in np.linspace(optim_sol2[0], nodes2[0][1], 50)]
        # y_slope2 = [m_p21*x + b2 for x in x_slope2]

        optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
        
        solved_heading = np.arctan2((optim_sol1[1]-nodes2[1][0]),(optim_sol1[1]-nodes2[0][0]))
        optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
       

        wpts, x_wpts, y_wpts = bez_to_wp(optimal_bez1, optimal_bez2, 15) 
        wpts_all1, wpts_all1x, wpts_all1y = bez_to_wp_single(optimal_bez1, 100) 
        wpts_all2, wpts_all2x, wpts_all2y = bez_to_wp_single(optimal_bez2, 100) 
        # print('WPTS', wpts_all1[0][0])
        valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 1', lr)
        # print("OPTIM SOL: ", optim_sol1[0], optim_sol1[1])
        # print(target_toa1, target_toa2)
        # plt.plot(optimal_bez1[0], optimal_bez1[1])
        # plt.plot(optimal_bez2[0], optimal_bez2[1])
        # plt.axis('equal')
        # plt.show()
        # print(valid1, valid2)
        # print('CHECK CHECK',m_p21*optim_sol1[0]+b2, optim_sol1[1])
        # print('CHECK CHECK',m_p12*optim_sol2[0]+b, optim_sol2[1])
        
        # if not np.isclose(m_p21*optim_sol1[0]+b2, optim_sol1[1], atol=0.001):
        #     valid1 = False
        # if not np.isclose(m_p12*optim_sol2[0]+b,optim_sol2[1], atol=0.001):
        #     valid2 = False
        # print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.01
            nodes1[1][1]+=-0.01
            nodes1[0][1]+=0.01*lr
            c+=1
        if valid2 == False:
            target_toa2+=0.01
            nodes2[0][1]+=-0.01*lr
            nodes2[1][1]+=0.01
            c+=1
        if c>1000:
            print('BEZ TOOK TOO LONG AND FAILED')
            break
        # if c>30:
       
        # fig = plt.figure(1)
        # ax = fig.add_subplot(111)
        # ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
        # ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
        # ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
        # ax.plot(x_slope, y_slope, color = 'purple', linestyle = 'dashed', label = 'G1 Continuity Line')
        # # ax.plot(x_slope, y_slope2, color = 'green', linestyle = 'dotted', label = 'G1 Continuity Line', zorder = 100)
        # # ax.plot(optimal_bezPre[0], optimal_bezPre[1], c = 'black')
        # # valid1 = validity_check(koz_x[0], koz_bot[1], koz_top[1], wpts_all1, corridor[0], 'Curve 1', lr)
        # # valid2 = validity_check(koz_x[0], koz_bot[1], koz_top[1], wpts_all2, corridor[0], 'Curve 1', lr)
        # # if not np.isclose(m_p21*optim_sol1[0]+b2, optim_sol1[1], atol=0.001):
        # #     valid1 = False
        # # if not np.isclose(m_p12*optim_sol2[0]+b,optim_sol2[1], atol=0.001):
        # #     valid2 = False
        # print('CHECK CHECK',m_p21*optim_sol1[0]+b2, optim_sol1[1])
        # print('CHECK CHECK',m_p12*optim_sol2[0]+b, optim_sol2[1])
        # print(valid1, valid2)
        # # print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
        # ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
        # ax.scatter(optim_sol2[0], optim_sol2[1])
        # ax.plot([corridor, corridor, corridor], [nodes1[1][0], nodes1[1][2], nodes2[1][2]], linestyle = '--')
        # ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
        # ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
        # ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
        # # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
        # ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
        # ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
        # ax.text(nodes1[0][2]-250/2,nodes1[1][2]-100/2,  r'Curve 1', fontsize = 12)
        # ax.text(nodes1[0][2]-250/2,nodes1[1][2]+100/2,  r'Curve 2', fontsize=12)
        
        # # y = [i for i in range(koz_top)]
        # # bx = [500 for i in range(koz_top)]
        # # bx2 = [1000 for i in range(koz_top)]
        # # ybot = [koz_bot for i in range(500)]
        # # bxbot = [i for i in range(500, 1000)]
        # # ytop = [koz_top for i in range(500)]
        # # y2 = [i for i in range(-50, 950)]
        # # ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
        # plt.plot([koz_x, koz_x], [koz_bot, koz_top], color = 'red', linestyle = '--')
        # plt.plot([nodes1[0][0], koz_x], [koz_bot, koz_bot], color = 'red', linestyle = '--')
        # plt.plot([nodes1[0][0], koz_x], [koz_top, koz_top], color = 'red', linestyle = '--', label = 'Emergency Vehicle Clearance Area')
        # plt.plot([corridor, corridor], [nodes1[1][0], nodes2[1][2]], color = 'orange', linestyle = 'dashed', label = 'Flight Corridor Bound')
        # ax.plot([-50, 100, 250], [711.62, 711.62, 711.62], linestyle = 'dashdot', alpha = 0.5, color = 'cyan')

        
        
        # ax.set_xlabel('X (m)')
        # ax.set_ylabel('Y (m)')
        # ax.grid(True)
        # ax.axis('equal')
        
        # ax.legend(loc = 'center left', fontsize = '8')
        # # plt.show()
        # # plt.pause(.5)
        # plt.show()

    # print(valid1, valid2)
    # latlong_XY = __WSG84_To_Meters_Single(latlon, [nodes1[0][0], nodes1[1][0]], Proj("EPSG:32667"))
    # print('OPTIMAL POINT 1', optim_sol1, 'OPTIMAL POINT 2', optim_sol2)
    
    
    nodes1 = [np.array([nodes1[0][0], optim_sol1[0], nodes1[0][2]]).flatten(),np.array([nodes1[1][0], optim_sol1[1], nodes1[1][2]]).flatten()]
    
    nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]
    end_time = time.time()
    print('TOTAL BEZ OPTIMIZATION TIME:', end_time-start_time)
    # plt.plot(optimal_bez1[0], optimal_bez1[1])
    # plt.plot(optimal_bez2[0], optimal_bez2[1])
    # plt.axis('equal')
    # plt.scatter(nodes1[0], nodes1[1])
    # plt.scatter(nodes2[0], nodes2[1])
    paths = [optimal_bez1, optimal_bez2]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    # print(x_wpts)
    # plt.scatter(x_wpts, y_wpts)

    print(f'\n\n\n TOA OF B1: {path_length(optim_sol1, [nodes1[0][0], nodes1[1][0]], [nodes1[0][2],nodes1[1][2]], 1)/velocity} TARGET TOA: {target_toa1}, RATIO: {(path_length(optim_sol1, [nodes1[0][0], nodes1[1][0]], [nodes1[0][2],nodes1[1][2]], 1)/velocity)/target_toa1}')
    print(f'\n\n\n TOA OF B2: {path_length(optim_sol2, [nodes2[0][0], nodes2[1][0]], [nodes2[0][2],nodes2[1][2]], 1)/velocity} TARGET TOA: {target_toa2}, RATIO: {(path_length(optim_sol2, [nodes2[0][0], nodes2[1][0]], [nodes2[0][2],nodes2[1][2]], 1)/velocity/target_toa2)}')
    
    optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
    optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)
    bez_data = {
        'Bez1X': optimal_bez1[0],
        'Bez1Y': optimal_bez1[1],
        'B1_NodeX': nodes1[0],
        'B1_NodeY': nodes1[1],
        'Bez2X': optimal_bez2[0],
        'Bez2Y': optimal_bez2[1],
        'B2_NodeX': nodes2[0],
        'B2_NodeY': nodes2[1],
    }
    # other_data = {
    #     'koz_x': koz_x,
    #     'koz_top': koz_top,
    #     'koz_bot': koz_bot,
    #     'corridor': corridor,
    #     'op_corridor': corridor-228,
    # }
    # scipy.io.savemat('4XBez.mat', bez_data)

    # print('SAVED')
    BezierData = {
        'B1':{
            'p0x': nodes1[0][0],
            'p0y': nodes1[1][0],
            'p1x': optim_sol1[0],
            'p1y': optim_sol1[1],
            'p2x': nodes1[0][2],
            'p2y': nodes1[1][2],
            'length': optim1_length,
            'Bez_ID': 'Optimal Bez 1',
            'TOA': target_toa1,
            'orig_TOA': orig_toa1,
            'travel_time': optim1_length/velocity,
            'timeStamp': timeStamp,
            'velocity': velocity,
            
            'ACID': acid,
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        },
        'B2':{
            'p0x': nodes2[0][0],
            'p0y': nodes2[1][0],
            'p1x': optim_sol2[0],
            'p1y': optim_sol2[1],
            'p2x': nodes2[0][2],
            'p2y': nodes2[1][2],
            'length': optim2_length,
            'Bez_ID': 'Optimal Bez 2',
            'TOA': target_toa2,
            'orig_TOA': orig_toa2,
            'travel_time': optim2_length/velocity,
            'timeStamp': timeStamp,
            'velocity': velocity,
            'home_lat': homell[0],
            'home_lon': homell[1],
            'ACID': acid,
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        }
    }

    return wpts, x_wpts, y_wpts, optimal_bez1, optimal_bez2, nodes1, nodes2, BezierData 

def EntryExitOutside(nodes1, nodes2, pos, velocity, lr, id, timeStamp, expnum, exptype, spacing):
    start = time.time()
    vel_knots = velocity*1.94384

    optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
    optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)
    
    # if id!= 0:
    ba, t_entry = solve_optimEntry(np.deg2rad(25), np.deg2rad(73), np.deg2rad(15), nodes1, vel_knots, pos, lr)
    # else:
    #     ba, t_entry = solve_optimEntry(np.deg2rad(15), np.deg2rad(73), np.deg2rad(15), nodes1, vel_knots, pos, lr)

    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry+0.0025)
    
    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))

    
    x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en], pos, lr)
    # plt.plot(x_entry, y_entry,color='cyan')
    optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
    # plt.plot(optimal_bez1[0], optimal_bez1[1], color='black')

    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], 200, t_entry)
    # plt.plot(partial_bez1[0], partial_bez1[1],color='magenta')
    # plt.show()


    t_start2 = 0.35678391959798994

    if lr == -1:
        baEx, exit_t = solve_optimExit([np.deg2rad(25), 0.45], np.deg2rad(73), np.deg2rad(15), t_start2, nodes2, velocity, lr)
    else:
        baEx, exit_t = solve_optimExit([np.deg2rad(28), 0.45], np.deg2rad(73), np.deg2rad(15), t_start2, nodes2, velocity, lr)

    # else:
    #     baEx, exit_t = solve_optimExit([np.deg2rad(60), 0.55], np.deg2rad(60), np.deg2rad(20), t_start2, nodes2, velocity, lr)
    print('EXIT BANK', np.rad2deg(baEx), 'EXIT T:', exit_t)

    x_int_ex, y_int_ex = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t)
    
    x_int_ex2, y_int_ex2 = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t-0.0025)
    
    req_ex = np.rad2deg(np.arctan2(y_int_ex-y_int_ex2,x_int_ex-x_int_ex2))

    x_exit, y_exit, central_angle_ex, exitLength, exitTOA, h_ex, k_ex = exitPath(velocity=57.412, t_exit=exit_t,
                                                                                    ba = baEx, intersect=[x_int_ex, y_int_ex],
                                                                                    nodes = [[nodes2[0][0],nodes2[1][0]],
                                                                                            [nodes2[0][1],nodes2[1][1]],
                                                                                            [nodes2[0][2],nodes2[1][2]]], lr = lr)
    # x_exit[0] = 750
    # y_exit[0] = 2350

    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], 200, exit_t+.025)
    
    pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_entry)
    pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]], P2=[nodes2[0][2],nodes2[1][2]], t=exit_t)
    total_toa = entryTOA+(pb1_length/velocity)+(pb2_length/velocity)+exitTOA
    print('TOTAL TOA:', total_toa)

    paths = [partial_bez1, partial_bez2]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    # plt.figure()
    # plt.plot(partial_bez1[0], partial_bez1[1], color = 'magenta')
    # plt.plot(partial_bez2[0], partial_bez2[1], color = 'magenta')
    # plt.plot(x_entry, y_entry, color = 'cyan')
    # plt.plot(x_exit, y_exit, color = 'orange')
    # plt.axis('equal')
    # plt.scatter(x_wpts, y_wpts)
    # plt.scatter(h_ex, k_ex)
    # plt.show()
    end = time.time()
    print('TOTAL INTERCEPTION PATH GENERATION TIME:', end-start)

    DubinsData = {
        'Entry':{
            'path_type': 'Entry',
            'bez_intx': x_int_en,
            'bez_inty': y_int_en,
            'nom_intx': pos[0],
            'nom_inty': pos[1],
            'bez_t': t_entry,
            'intersect_heading': req_ent,
            'h': h_c,
            'k': k_c,
            'bank_angle': np.rad2deg(ba[0]),
            'tr': (velocity*1.9438445)**2/(11.26*math.tan(np.deg2rad(ba))),
            'timeStamp': timeStamp,
            'ACID': id,
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        },
        'Exit':{
            'path_type': 'Exit',
            'bez_intx': x_int_ex,
            'bez_inty': y_int_ex,
            'nom_intx': x_exit[0],
            'nom_inty': y_exit[0],
            'bez_t': exit_t,
            'intersect_heading': req_ex,
            'h': h_ex,
            'k': k_ex,
            'bank_angle': np.rad2deg(baEx),
            'tr': (velocity*1.9438445)**2/(11.26*math.tan(np.deg2rad(ba))),
            'timeStamp': timeStamp,
            'ACID': id,
            'ExpNum': expnum,
            'Category': 'Fleet Aircraft',
            'ExpType': exptype,
            'ax_spacing': spacing
        }
    }

    return [x_entry, y_entry, req_ent, np.rad2deg(ba)], [x_exit, y_exit, req_ex, np.rad2deg(baEx)], total_toa, x_wpts, y_wpts, DubinsData

def entryPath(velocity, ba, intersect, pos, lr):

    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
    tr*=0.3048
    print('TR', tr)
    h, k  = lr*tr, pos[1]

    x_entry = [i for i in np.linspace(0, intersect[0], 200)]

    y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
    # y_entry[-1] = intersect[1]
    ar, ad = central_angle([h, k], [x_entry[0], y_entry[0]], [x_entry[-1], y_entry[-1]])
    # entryLength = tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
    # entryTOA = entryLength/velocity
    entryLength = 2*pi*tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
    entryTOA = entryLength/(velocity**2/(9.81*np.tan(np.deg2rad(ba))))
    print('ENTRY ToA: ', entryTOA)

    return x_entry, y_entry, ar, entryLength, entryTOA, h, k

def exitPath(velocity, t_exit, ba, intersect, nodes, lr):
    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(ba))
    tr*=0.3048

    Bx = lambda t: nodes[1][0] + (nodes[0][0] - nodes[1][0]) * (1 - t)**2 + (nodes[2][0] - nodes[1][0]) * t**2
    By = lambda t: nodes[1][1] + (nodes[0][1] - nodes[1][1]) * (1 - t)**2 + (nodes[2][1] - nodes[1][1]) * t**2

    bezHead = np.arctan2(By(t_exit)-By(t_exit-0.01), Bx(t_exit)-Bx(t_exit-0.01))

    h, k = intersect[0] - tr*math.cos(bezHead), intersect[1]+tr*math.sin(bezHead)
    
    print('H AND K', h, k)
    # print(h)
    # print(k)
    x_exit = [i for i in np.linspace(0, intersect[0], 200)]
    y_exit = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_exit]
    # y_exit[0] = 650
    # y_exit[0] = k
    # print(y_exit)
    c = 0
    # for i in range(0, len(y_exit)):
    #     if math.isnan(y_exit[i]):
    #         y_exit[i] = k-(i*2.5)
    # print(y_exit)
    ar, ad = central_angle(center = [h, k], point1=[x_exit[-1], y_exit[-1]], point2=[x_exit[0], y_exit[1]])

    # exitLength = tr*ar

    # exitTOA = exitLength/velocity
    ar, ad = central_angle(center = [h, k], point1=[x_exit[-1], y_exit[-1]], point2=[x_exit[0], y_exit[0]])
    # print(ar, ad, ba)
    exitLength = 2*pi*tr*ar

    exitTOA = exitLength/(velocity**2/(9.81*np.tan(ba)))
    
    print('EXIT ToA:', exitTOA)
    
    return x_exit, y_exit, ar, exitLength, exitTOA, h, k



if __name__ == "__main__":
 
    velocity  = 57.3024 #m/s
    turn_rate = np.deg2rad(4.50) # RAD/s
    turn_radius = 111.6**2/(11.26*math.tan(np.deg2rad(30)))
    turn_radius*=0.3048
    h = 275
    koz_bot = 15
    koz_top = h-15

    
    ev_toa = h/64.008
    
    c=0
    target_toa1 = 1.5*ev_toa
    target_toa2 = 1.5*ev_toa
    uav_head = np.deg2rad(90)
    lr = 1
    if lr == -1:
        corridor = 0
        koz_x = 500
        nodes1 = [np.array([750, 1400, 25]).flatten(),np.array([-50, 10, 450]).flatten()]
        nodes2 = [np.array([25, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()]
    else:
        corridor = 457
        koz_x = 305
        nodes1 = [np.array([229, 450, 450]).flatten(),np.array([0, -10, h/2]).flatten()]
        nodes2 = [np.array([450, 450, 229]).flatten(),np.array([h/2, h, h]).flatten()]
    
 
    print("VEHICLE VELOCITY:", velocity)
    # print("VEHICLE TURN RADIUS:", turn_radius)
    # print("TARGET LENGTH: ", target_length)
    valid1 = False
    valid2 = False
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    m_p21 = (nodes2[1][1] - nodes1[1][2])/(nodes2[0][1] - nodes1[0][2])
    print(m_p21)
    b2 = nodes2[1][1] - m_p21*nodes2[0][1]
    lines_coeffs2 = [m_p21, b2]
    bezier_variny_1 = manual_bez([nodes1[0][0],nodes1[1][0]], [nodes1[0][1],nodes1[1][1]],[nodes1[0][2],nodes1[1][2]], 20)
    bezier_variny_2 = manual_bez([nodes2[0][0],nodes2[1][0]], [nodes2[0][1],nodes2[1][1]],[nodes2[0][2],nodes2[1][2]], 20)
    #OPTIMIZE TEST:
    start = time.time()
    while not valid1 or not valid2:

        optim_sol1, curv1 = solve_optim1(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=target_toa1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs2)
        
        m_p12 = (nodes1[1][2]-optim_sol1[1])/(nodes1[0][2]-optim_sol1[0])
        # print(m_p12)
        mark = True
        #y = mx+b
        if np.isclose(nodes1[0][2], optim_sol1[0], atol = 0.01) == False:
            mark = False

        b = optim_sol1[1] - m_p12*optim_sol1[0]
        lines_coeffs = [m_p12, b]
        x_slope = [i for i in np.linspace(427, 472, 50)]
        y_slope = [m_p12*x + b for x in x_slope]
        
        
        # print(lines_coeffs)

        optim_sol2, curv2 = solve_optim2(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=target_toa2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs)
        m_p21 = (optim_sol2[1] - nodes1[1][2])/(optim_sol2[0] - nodes1[0][2])
        b2 = optim_sol2[1] - m_p21*optim_sol2[0]
        lines_coeffs2 = [m_p21, b2]
        if mark == True:
            x_slope = [nodes1[0][2] for i in np.linspace(427, 472, 50)]
            y_slope = [i for i in np.linspace(optim_sol1[1]-50, optim_sol2[1]+50, 50)]

        
        optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
        # plt.plot(optimal_bez1[0], optimal_bez1[1])
        # plt.show()
        solved_heading = np.arctan2((optim_sol1[1]-nodes2[1][0]),(optim_sol1[1]-nodes2[0][0]))
        optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
       

        wpts, x_wpts, y_wpts = bez_to_wp(optimal_bez1, optimal_bez2, 15) 
        wpts_all1, wpts_all1x, wpts_all1y = bez_to_wp_single(optimal_bez1, 100) 
        wpts_all2, wpts_all2x, wpts_all2y = bez_to_wp_single(optimal_bez2, 100) 
        valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 1', lr)
        # print("OPTIM SOL: ", optim_sol1[0], optim_sol1[1])
        # print(target_toa1, target_toa2)

        # print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.1
            nodes1[1][1]+=-.5
            nodes1[0][1]+=.5*lr
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            nodes2[0][1]+=-.5*lr
            nodes2[1][1]+=.5
            c+=1
        if c>30:
       
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
            ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
            ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
            # ax.plot(optimal_bezPre[0], optimal_bezPre[1], c = 'black')
            # valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
            # valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
            # print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
            ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
            # ax.scatter(optim_sol2[0], optim_sol2[1], marker = '*', s = 1000, zorder = 1000)
            ax.plot([corridor, corridor, corridor], [nodes1[1][0], nodes1[1][2], nodes2[1][2]], linestyle = '--')
            ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
            ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
            ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
            # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
            ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
            ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
            ax.text(nodes1[0][2]-250,nodes1[1][2]-100,  r'Curve 1')
            ax.text(nodes1[0][2]-250,nodes1[1][2]+100,  r'Curve 2')
            y = [i for i in range(koz_top)]
            bx = [153 for i in range(koz_top)]
            bx2 = [305 for i in range(koz_top)]
            ybot = [koz_bot for i in range(153)]
            bxbot = [i for i in range(152, 305)]
            ytop = [koz_top for i in range(153)]
            y2 = [i for i in range(-15, 290)]
            ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
            ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')

            ax.grid(True)
            ax.axis('equal')

            ax.legend(loc = 'center left', fontsize = '8')
            plt.show()
        
    end = time.time()-start
        
    
    
    solved_heading1 = np.arctan2((optim_sol1[1]-nodes1[1][0]),(optim_sol1[1]-nodes1[0][0]))
    solved_heading2 = np.arctan2((optim_sol2[1]-nodes2[1][0]),(optim_sol2[1]-nodes2[0][0]))
    print(curv1, path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)/velocity)
    print("HEADING TO P1", solved_heading)
    print("REQUESTED HEADING: ", uav_head)
    print('OPTIMAL POINT 1', optim_sol1)
    print('OPTIMAL POINT 2:', optim_sol2)
    optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
    optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)
    
    print("SOLVED LENGTH WITHOUT ENTRY: ", optim1_length + optim2_length)

 
    print(target_toa1, target_toa2)

    print("ToA WITHOUT ENTRY:", (optim1_length + optim2_length)/ velocity)
    
    vel_knots = 111.6
    nodes1 = [np.array([nodes1[0][0], optim_sol1[0], nodes1[0][2]]).flatten(),np.array([nodes1[1][0], optim_sol1[1], nodes1[1][2]]).flatten()]
    ba, t_entry = solve_optimEntry(np.deg2rad(25), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots)
    print(np.rad2deg(ba), t_entry)


    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    
    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
    print('REQUIRED HEADING AT ENTRY:', req_ent)
    
    
    x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en])
    act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[198], x_entry[-1]-x_entry[198]))
    print('ACTUAL HEADING AT INTERSECT:', act_ent)
    
    print('DIFF AT ENTRY:', req_ent-act_ent)

    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200, t_entry)
    
    t_start2 = 0.35678391959798994

    nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]
    baEx, exit_t = solve_optimExit([np.deg2rad(25), t_start2], np.deg2rad(30), np.deg2rad(20), t_start2, nodes2, velocity)
    print('EXIT BANK', np.rad2deg(baEx) )

    x_int_ex, y_int_ex = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t)
    
    x_exit, y_exit, central_angle_ex, exitLength, exitTOA, h_ex, k_ex = exitPath(velocity=velocity, t_exit=exit_t,
                                                                                 ba = baEx, intersect=[x_int_ex, y_int_ex],
                                                                                 nodes = [[nodes2[0][0],nodes2[1][0]],
                                                                                            [optim_sol2[0],optim_sol2[1]],
                                                                                            [nodes2[0][2],nodes2[1][2]]])
    x_exit[0] = 229
    # y_exit[0] = 487
    # print(y_exit[0])
    # y_exit[0] = 2090
    head_ex = np.rad2deg(np.arctan2(y_exit[0]-y_exit[1], x_exit[0]- x_exit[1]))
    print('ACTUAL HEADING AT EXIT:', head_ex)
    # x_exit, y_exit, exit_int, exit_t, exitTOA, exit_bank, exitLength = exitPath2(velocity, t_exit, optimal_bez2)
    print('EXIT T:', exit_t)
    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200, exit_t)
    # print(x_exit)
    # print(x_entry[0])
    entry_path = [x_entry, y_entry]
    # print(entry_path)
    # print(entry_path[0][0])
    # exit_path = [x_exit, y_exit]
    # paths = [entry_path, partial_bez1, partial_bez2, exit_path]
    # print(len(paths))
    # print(paths[0][0])
    # wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    # print(wpts)


    pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_entry)
    pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]], P2=[nodes2[0][2],nodes2[1][2]], t=t_start2)

    y = [i for i in range(15, koz_top)]
    bx = [153 for i in range(15, koz_top)]
    bx2 = [305 for i in range(15, koz_top)]
    ybot = [koz_bot for i in range(152, 305)]
    bxbot = [i for i in range(152, 305)]
    ytop = [koz_top for i in range(153)]
    y2 = [i for i in range(-396, 640)]
    xwall = [0 for i in range(-396, 640)]
    xwall2 = [457 for i in range(-396, 640)]
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    # ax.scatter(x_wpts, y_wpts, zorder = 100, marker = '^')
    ax.plot(x_exit, y_exit, color = 'orange', label = 'Exit Path')
    # ax.scatter(h_ex, k_ex, marker = 's', color = 'green', label = 'Exit Circle Center')
    ax.scatter(x_exit[-1], y_exit[-1], color = 'red', marker = '*', s = 100, zorder = 100, label = 'Exit Point')
    # print('EXIT COORDS:', x_exit, y_exit)
    ax.scatter(x_exit[0], y_exit[0], color = 'purple', marker = '*', s = 100, zorder = 100, label = 'Interception Point')
    # ax.scatter(pb2[0], pb2[1], c = 'magenta', zorder = 100)
    ax.plot([122, 305, 488], [h/2, h/2, h/2], linestyle = 'dashdot', alpha = 0.5)
    # ax.plot(optimal_bez2[0], optimal_bez2[1])
    ax.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial QBC')
    ax.plot(partial_bez2[0], partial_bez2[1], c = 'magenta')#, label = 'Partial QBC')
    # ax.plot(optimal_bez1[0],optimal_bez1[1], c='black', label='Quadratic Bezier curve')
    # ax.scatter(x_wpts, y_wpts, label = 'Waypoints', marker = '^', color = 'green')
    # print(wpts_all1[0])
    # ax.scatter(wpts_all1x, wpts_all1y, zorder = 30)
    print('ENTRY ANGLE:', np.rad2deg(np.arctan2(partial_bez1[1][1] - y_entry[-1], partial_bez1[0][1] - x_entry[-1])))
    ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
    
    # ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')

    # ax.scatter(h_c, k_c, color = 'green', marker = 's', label = 'Entry Circle Center')
    ax.scatter(x_entry[0], y_entry[0], marker = '^', color = 'black', label = 'Fleet Aircraft Start Position')


    ax.plot(x_entry, y_entry, label = 'Entry Arc', color = 'cyan')
    ax.scatter(x_entry[-1], y_entry[-1], marker = '*', color = 'purple', label = 'Intersection Point', s = 100, zorder = 100)
    # print(y_circ[0])
    

    print("ENTRY TOA:", entryTOA)
    print("EXIT TOA:", exitTOA)
    # print("SOLVED BEZIER LEGNTH: ", pb1_length + pb2_length)
    # print("PARTIAL BEZIER TRAVEL TIME:", (pb1_length+pb2_length)/velocity)

    # print("SOLVED LENGTH WITH ENTRY/EXIT: ", pb1_length + optim2_length + entryLength+exitLength)
    print("TOTAL TOA:", entryTOA + exitTOA + (pb1_length+pb2_length)/velocity)


    ax.plot(x_slope, y_slope, color = 'purple', linestyle = 'dashed', label = 'G1 Continuity Line')
    
    ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
    ax.text(nodes1[0][0]+0.25,nodes1[1][0]-10,  r'$\bf{p_{01}}$')
    ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
    ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
    # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
    ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
    ax.text(nodes2[0][2]+0.25,nodes2[1][2]+6.5,  r'$\bf{p_{22}}$')
    ax.text(nodes1[0][2]-137,nodes1[1][2]-20,  r'Curve 1', fontsize = 6)
    ax.text(nodes1[0][2]-137,nodes1[1][2]+45.7,  r'Curve 2', fontsize = 6)

    # ax.text(mx+0.25, my, r'$\bf{m}$')
    # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
    # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
    ax.plot(bx, y, color = 'red', linestyle = '--')
    ax.plot(bx2, y, color = 'red', linestyle = '--')

    ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
    ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
    ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':', color = 'orange')
    ax.plot(xwall2, y2, linestyle = ':', color = 'orange')
    # y_nom = [i for i in range(-1282, int(y_exit[0]))]
    # x_nom = [750 for i in y_nom]
    # ax.plot(x_nom, y_nom, color = 'black', label = 'Nominal Path')

    ax.plot([229, 229], [-396, y_exit[0]+50], color = 'green', label = 'Nominal Path')

    # ax.scatter(mx,my)
    # ax.scatter(center1x,center1y)
    # ax.scatter(center2x,center2y)
    # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
    # ax.add_patch(Circle((center2x, center2y), r2, color='black', fill=False))
 
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.legend(loc = 'center right', fontsize = '7')
    # plt.ylim(optim_sol1[1]-150, optim_sol2[1]+150)
    # plt.axes.set_ylim(-150, 600)
    plt.show()

#line types for, dashed, solid, dots
#