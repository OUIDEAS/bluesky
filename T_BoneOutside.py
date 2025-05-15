import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy.io
import sympy as sp
from scipy.optimize import fsolve
from pyproj import Proj
from utm import from_latlon as llutm
from utm import to_latlon as utmll
import sys
sys.path.append("C:/Users/Michael/Desktop/bluesky")
import optimize
from collections import defaultdict
import matplotlib.patches as patches
import matplotlib.lines as mlines
from optimize import get_line_o, get_travel_xy, total_travel_o, curvature_o, path_length_o, toa_diff_o, find_bez_xy_o, manual_bez_partial_o, manual_bez_o, entryPath_o, central_angle_o

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

def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

def manual_bez_xy(P0, P1, P2, points):
    t = np.linspace(0, 1, points)
    x = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    y = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    return np.array([x, y])

def find_bez_xy(P0, P1, P2, t):
    # x = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    # y = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    Bx = lambda z: P1[0] + (P0[0] - P1[0]) * (1 - z)**2 + (P2[0] - P1[0]) * z**2
    By = lambda z: P1[1] + (P0[1] - P1[1]) * (1 - z)**2 + (P2[1] - P1[1]) * z**2
    x = Bx(t)
    y = By(t)
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
    C = bx**2 + by**2
    b=B/(2.0*A)
    c=C/A
    u=t+b
    k=c-(b*b)
 
    L=0.5*np.sqrt(A)*((u*np.sqrt((u*u)+k))
        -(b*np.sqrt((b*b)+k))
        +(k*np.log(np.abs((u+np.sqrt((u*u)+k))/(b+np.sqrt((b*b)+k)))))
    )
    
    return L

def toa_diff(P1, P0, P2, t, toa, v):
    ax = P0[0] - 2*P1[0] + P2[0]
    ay = P0[1] - 2*P1[1] + P2[1]
    bx = 2*P1[0] - 2*P0[0]
    by = 2*P1[1] - 2*P0[1]
    A = 4 * (ax**2 + ay**2)
    B = 4 * (ax*bx + ay*by)
    C = bx**2 + by**2
    b=B/(2.0*A)
    c=C/A
    u=t+b
    k=c-(b*b)
 
    L=0.5*np.sqrt(A)*((u*np.sqrt((u*u)+k))
        -(b*np.sqrt((b*b)+k))
        +(k*np.log(np.abs((u+np.sqrt((u*u)+k))/(b+np.sqrt((b*b)+k)))))
    )
    diff = np.abs(L - toa*v)
    # print(f'TOA DIFF: {diff}')
    return diff


    # return L

 
def solve_optim1(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr, line):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    # if g == 'III' or g == 'IV':
    #     cons = (
    #             {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
    #             {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
    #             # {'type': 'ineq', 'fun': lambda x: x[0]-P2[0]},
    #             {'type': 'ineq', 'fun': lambda x: x[1]-P2[1]},
    #             # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
    #             {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
    #             ) 
    # else:
    cons = (
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
            # {'type': 'ineq', 'fun': lambda x: x[0]-P2[0]},
            {'type': 'ineq', 'fun': lambda x: P0[0]-x[0]-1000},
            {'type': 'ineq', 'fun': lambda x: P2[1]-x[1]-100},
            # {'type': 'ineq', 'fun': lambda x: np.sqrt((x[0] - P2[0])**2 + (x[1] - P2[1])**2) - 500},
            # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
            {'type': 'eq', 'fun': lambda x: (line[0]*x[0] + line[1])-x[1]}
            ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
 
    return val.x , curvature(P0,val.x,P2)


def solve_optim2(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr, line):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))

    cons = (
            
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
            {'type': 'ineq', 'fun': lambda x: x[1] - P0[1]},
            {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
        ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

    return val.x , curvature(P0,val.x,P2)

def solve_optimEntry(guess, max_bank,  min_bank, nodes, velocity, pos, lr, h, show):
    def path_cost(ba, nodes, velocity, lr, h, show, diff):
        diff = find_diff_entry(ba, nodes, velocity, pos, lr, h, show, diff[0])
        # if diff[1] == 0.01:  # No valid intersection
        #     return 1e6
        return np.abs(diff[0])
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}  
    )
    # val= minimize(path_cost,guess,(nodes, velocity, lr, h, show), method='SLSQP', tol=1E-3, constraints=cons)
    # print(f"Optimization result: Bank angle = {np.rad2deg(val.x[0])} degrees")
    diff = [1e6]
    val = minimize(path_cost, guess, args=(nodes, velocity, lr, h, show, diff),
                      method='SLSQP', bounds=[(min_bank, max_bank)], constraints=cons,
                      options={'xtol': 1e-3, 'gtol': 1e-3, 'maxiter': 100})

    t_val = find_diff_entry(val.x, nodes, velocity, pos, lr, h, show, diff[0])[1]
    if t_val <0:
        t_val = 0
    return val.x, t_val

def solve_optimExit(guess, max_bank, min_bank, min_t, nodes, velocity, lr):
    #t is also being minimized, so part of constraints
    def path_cost(guess, nodes, velocity, lr, diff):
        diff = find_diff_exit(guess, nodes, velocity, lr, diff[0])
        # print('DIFF FOR EXIT:', diff[0])
        return np.abs(diff[0]) 
    diff = [1e6]
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank},
            {'type': 'ineq', 'fun': lambda x: x[1] - min_t},
            {'type': 'ineq', 'fun': lambda x: 1-x[1]},
    )
    val =  minimize(path_cost, guess, (nodes, velocity, lr, diff), method = 'SLSQP', tol = 1E-10, constraints=cons)
    y_val = find_diff_exit(val.x, nodes, velocity, lr, diff[0])[2]
    # print(val)
    return val.x[0], val.x[1], y_val

def central_angle(center, point1, point2):
    """
    Computes the central angle between two points on a circle.

    Parameters:
        center (tuple): (x, y) coordinates of circle center.
        point1 (tuple): First point (x, y) on the circle.
        point2 (tuple): Second point (x, y) on the circle.

    Returns:
        angle_radians (float): Angle in radians.
        angle_degrees (float): Angle in degrees.
    """
    # Convert to NumPy arrays
    center = np.array(center)
    point1 = np.array(point1)
    point2 = np.array(point2)

    # Compute vectors from center to points
    v1 = point1 - center
    v2 = point2 - center

    # Compute central angle using the dot product formula
    dot_product = np.dot(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle_radians = np.arccos(np.clip(dot_product / norms, -1.0, 1.0))

    # Convert to degrees
    angle_degrees = np.rad2deg(angle_radians)

    return angle_radians, angle_degrees

def find_diff_entry(ba, nodes, velocity, pos, lr, head, show, mindiff):
    pi = np.pi
    t_guess = 0.0075
    # mindiff = 1e6
    # print(np.rad2deg(head))
    
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    tr = velocity**2 / (11.26*math.tan(ba[0]))
    tr*=0.3048

    # h = pos[0] + lr*tr*np.cos(head-180)
    # k = pos[1] - lr*tr*np.sin(head-180)

    h = pos[0]+ lr*tr
    k = pos[1]
    # print(pos)
    # if head == 0 or head == 2*np.pi:
    #     h = pos[0] - lr*tr*np.sin(head)
    #     k = pos[1] - lr*tr*np.cos(head)
    circle_eq = lambda t: (Bx(t)-h)**2+(By(t)-k)**2 - tr**2
    # print(f'H:{h}, k:{k}, tr:{tr}')
    
    
    if show:
        path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
        fig, ax = plt.subplots()
        plt.scatter(h, k, color = 'blue', marker = 's')
        plt.scatter(pos[0], pos[1])
        plt.axis('equal')
        circ = plt.Circle((h, k), tr, fill = False)
        ax.add_patch(circ)
        plt.plot(path[0], path[1], color = 'black', linewidth = 2)
        # plt.plot(x_l, y)
        plt.axis('equal')
        plt.arrow(pos[0], pos[1], 2500*np.cos(head), 2500*np.sin(head), width=50,
                    length_includes_head=False,
                    head_width=1000,
                    head_length=1000,
                    head_starts_at_zero=False,
                    facecolor='black',
                    edgecolor='black',
                    label = f'Heading of {np.rad2deg(head)-180}')
        plt.show()
    # dist_to_bezier = np.sqrt((Bx(0) - h) ** 2 + (By(0) - k) ** 2)

    # Use this distance to get an initial guess closer to the intersection
    # t_guess = np.clip(dist_to_bezier / 2, 0.1, 0.9)
    S = fsolve(circle_eq, t_guess)
    # result = root(circle_eq, t_guess, method = 'hybr')
    # print(S)
    

    for i in S:

        if i>=0 and i <= 1:
            index = int(np.round(float(i*200)))
            # print(index)
            if index == 200:
                index = 199
            if index == 0:
                index = 5
            bez_angle = np.arctan2(By(i+.01)-By(i), Bx(i+0.01)-Bx(i))
            # print(f'BEZ ANGLE: {np.rad2deg(bez_angle)}')
            if np.rad2deg(bez_angle) < 0:
                bez_angle+=2*pi
            # print('BEZ X:', Bx(i))
            if lr == 1:
                # print('LR + 1 LAND')
                # y = np.linspace(pos[1], By(i), 200)
                x_l = np.linspace(pos[0], Bx(i), 200)
                y = k+np.sqrt(tr**2 - (x_l-h)**2)
                # x_l = h-np.sqrt(tr**2 - (y-k)**2)
                # x_2 = h+np.sqrt(tr**2 - (y-k)**2)
                # print(x_l)
                # print(x_2)
            else:
                x_l = np.linspace(pos[0], Bx(i), 200) #Space between nominal path and bez
            # print(x_l)
                y = k+np.sqrt(tr**2 - (x_l-h)**2) #arc created by turn

                # y[0] = pos[1]
            # print(x[0], x[0]-h, x[0]-h**2, )
            # print(y)
            if show:
                fig, ax = plt.subplots()
                plt.scatter(h, k)
                plt.scatter(pos[0], pos[1])
                plt.axis('equal')
                circ = plt.Circle((h, k), tr, fill = False)
                ax.add_patch(circ)
                plt.plot(path[0], path[1])
                plt.plot(x_l, y, color = 'orange')
                # plt.plot(x_2, y, color = 'red')
                plt.axis('equal')
                plt.show()

            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            if bez_angle >= 3*pi/2 and int_angle <= pi/2:
                int_angle+=2*pi
            # print(f'INT ANGLE {np.rad2deg(int_angle)}')
            if math.isnan(int_angle):
                int_angle = np.arctan2(y[198]-y[197], x_l[198] - x_l[197])
            if np.deg2rad(int_angle) <0:
                int_angle+=2*pi
            # print(f'BEZ ANGLE {np.rad2deg(bez_angle)}, INT ANGLE {np.rad2deg(int_angle)}, DIFFERENCE {np.abs(np.rad2deg(bez_angle-int_angle))}')
            diff = np.abs(bez_angle-int_angle)
            # print(np.rad2deg(int_angle))
            if diff < mindiff and np.abs(y[-1] - By(i))<=5 and np.abs(x_l[-1]-Bx(i))<=5:
                mindiff = diff
                if show:
                    plt.plot(path[0], path[1])
                    plt.plot(x_l, y)
                    plt.scatter(path[0][index], path[1][index], color = 'yellow', marker = '*', s = 100)
                    plt.scatter(x_l[-1], y[-1], color = 'purple', marker = '*', s = 100)
                    plt.axis('equal')
                    plt.show()
                t_final = i
            else:
                # If no valid intersection found, dynamically penalize based on distance from Bézier curve
                dist_to_curve = np.sqrt((Bx(i) - x_l[-1]) ** 2 + (By(i) - y[-1]) ** 2)
                mindiff = max(mindiff, dist_to_curve * 100)  # Apply penalty
                t_final = i
        else:
            t_final = i 
    # print(f'BANK: {np.rad2deg(ba)}, T: {t_final}')
    # print('MIN DIFF:', mindiff)
    return [mindiff, t_final]

def find_diff_exit(guess, nodes, velocity, lr, mindiff):
    ba = guess[0]
    t_guess = guess[1]
    vk = velocity * 1.94384
    # print(guess, ba, t_guess)
    # mindiff = 100
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2 #Bezier x
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2 #Bezier y
    bezHead = np.arctan2(By(t_guess)-By(t_guess-0.01), Bx(t_guess)-Bx(t_guess-0.01))

    nom_y = np.linspace(1000, 4000, 200)
    nom_x = [0 for i in nom_y]
    p = []

    angle_to_rotate = 90-np.rad2deg(bezHead)
    
    for i in range(0, 200):
        p.append([nom_x[i], nom_y[i]])
    p = np.array(p) - np.array([Bx(t_guess), By(t_guess)])  # Shift to origin
    
    rotated_pts = rotate_bez(np.array([nom_x, nom_y]), angle_to_rotate, np.zeros(2))  # Rotated nominal path
    # rotated_pts += np.array([Bx(t_guess), By(t_guess)])
    m = (rotated_pts[1][100] - rotated_pts[1][0])/(rotated_pts[0][100] - rotated_pts[0][0])
    # rotated_pts = rotate_points(p, angle_to_rotate, np.array([0, 0]))
    # print(rotated_pts)
    rotated_nodes = rotate_bez(nodes, angle_to_rotate, np.array([0, 0]))
    Bx = lambda t: rotated_nodes[0][1] + (rotated_nodes[0][0] - rotated_nodes[0][1]) * (1 - t)**2 + (rotated_nodes[0][2] - rotated_nodes[0][1]) * t**2 #Rotated Bez
    By = lambda t: rotated_nodes[1][1] + (rotated_nodes[1][0] - rotated_nodes[1][1]) * (1 - t)**2 + (rotated_nodes[1][2] - rotated_nodes[1][1]) * t**2
    
    # tr = 111.6**2 / (11.26*math.tan(ba))
    tr = vk**2 / (11.26*math.tan(ba))
    tr*=0.3048
    h = lambda t: Bx(t) + lr*tr#*math.cos(bezHead)
    k = lambda t: By(t) #+ tr*math.sin(bezHead)

    path = manual_bez(P0 = [rotated_nodes[0][0], rotated_nodes[1][0]], #Rotated Bez
                      P1 = [rotated_nodes[0][1], rotated_nodes[1][1]],
                      P2 = [rotated_nodes[0][2], rotated_nodes[1][2]], 
                      points = 200)

    x_guess = Bx(1)-(tr/3)          
    y_guess = m*x_guess
    
    # print(x_guess, y_guess, m)
    circle_eq = lambda y: (x_guess-h(t_guess))**2+(y-k(t_guess))**2 - tr**2
    # y_guess = rotated_pts[1][100]
    S = fsolve(circle_eq, y_guess) #gives y intersection
    t_final = .5
    # print(vk)
    # y_l = [i for i in np.linspace(274, S, 200)]
    # x = [0 for i in y_l]

    # print(S)
    for i in S:
        if i > 0: 

            y = [z for z in np.linspace(By(t_guess), i, 200)]
            if lr == 1:
                x_l = [h(t_guess) - np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            else:
                x_l = [h(t_guess) + np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            # print(x_l, y)
            int_angle = np.arctan2(y[-1]-y[198], x_l[-1]-x_l[198])
            diff = np.abs(np.arctan2(rotated_pts[1][-1]-rotated_pts[1][0], rotated_pts[0][-1]-rotated_pts[0][0]) - int_angle)
            guess_deg = np.rad2deg(ba)
            # plt.title(f'{guess_deg} {t_guess}')
            # plt.plot(x, y_l, linestyle = 'dashed')
            # plt.arrow(Bx(t_guess), By(t_guess)+50, 75*np.cos(np.pi/2), 75*np.sin(np.pi/2), width=15,
            #     length_includes_head=False,
            #     head_width=50,
            #     head_length=50,
            #     head_starts_at_zero=False,
            #     facecolor='black',
            #     edgecolor='black', zorder = 100)
            # plt.text(-3000, 1750, 'Nominal Path', fontweight = 'bold', fontsize = 12)
                # label = f'Heading of {np.rad2deg(0*np.pi/2):0.2f}')
            # plt.plot(rotated_pts[0], rotated_pts[1], linestyle = '--', linewidth = 2, label = 'Nominal Path', color = 'blue')
            # plt.scatter(Bx(t_guess), By(t_guess), marker = '^', zorder = 100, label = 'Exit Point Guess', color = 'black')
            # print(np.hypot(h(t_guess)- Bx(t_guess), k(t_guess) - By(t_guess)), tr)
            # plt.plot(path[0], path[1], linewidth = 2, label = 'Rotated Bezier Curve', color = 'orange')
            # plt.plot(x_l, y, linewidth = 2, label = 'Exit Path', color = 'green')
            # plt.scatter(h(t_guess), k(t_guess), label = 'Circle Center', color = 'red')
            # plt.text(h(t_guess), k(t_guess)+100, f'{np.rad2deg(ba)}{chr(176)} Bank', fontweight = 'bold', fontsize = 13)
            # plt.text(Bx(t_guess)+50, By(t_guess), r'$\gamma$ = 0.75', fontsize = 13, fontweight = 'bold')
            # plt.legend(loc = 'upper right')
            # plt.axis('equal')
            # plt.xlim((-3500, -1000))
            # plt.xlabel('X (m)')
            # plt.ylabel('Y (m)')
            # plt.grid()
            # plt.show()
            h_x = h(t_guess)  # Compute h(t_guess) once for clarity
            index_match = np.argmin(np.abs(rotated_pts[0] - h_x))  # Find index of closest x value

            # print(f"Index of closest match: {index_match}, x value: {rotated_pts[0][index_match]}")
            if diff < mindiff and np.isclose(y[-1], y_guess, atol = 30) and np.isclose(x_l[-1], x_guess, atol = 30 ): #==0:# <= 0.01:
                # x_l = [i for i in np.linspace(750, Bx(t_guess))]
                # y = [k(t_guess)-np.sqrt(tr**2 - (x-h(t_guess))**2) for x in x_l]
                int_angle = np.arctan2(y[-1]-y[198], x_l[-1]-x_l[198])
                diff = np.abs(np.arctan2(rotated_pts[1][-1]-rotated_pts[1][0], rotated_pts[0][-1]-rotated_pts[0][0]) - int_angle)
                guess_deg = np.rad2deg(ba)
                # plt.title(f'{guess_deg} {t_guess}')
                # # plt.plot(x, y_l, linestyle = 'dashed')
                # plt.plot(path[0], path[1])
                # plt.plot(x_l, y)
                # plt.scatter(h(t_guess), k(t_guess))
                # plt.axis('equal')
                # plt.show()
                if diff < mindiff:
                    mindiff = diff
                    # print('BLARGY', diff)
                    t_final = t_guess             
                else:
                    dist_to_curve = np.sqrt((Bx(i) - x_l[-1]) ** 2 + (By(i) - y[-1]) ** 2)
                    mindiff = max(mindiff, dist_to_curve * 100)  # Apply penalty
                    t_final = i
                

    return [mindiff, t_final, i]

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

def LongLat_To_WSG84_Meters(waypoints, home):
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    waypoints = np.array(waypoints)
    asize = waypoints.shape
    
    waypoints_meters = []
    if(len(asize) > 1): #m x 3 matrix       
        for pt in waypoints:
            if len(pt) > 2:
                newpt = __WSG84_To_Meters_Single(pt[0], home)
                if(asize[1]==3):#altitude
                    newpt=[newpt[0],newpt[1],pt[1], pt[2]]
            else:
                newpt = __WSG84_To_Meters_Single(pt, home)
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

def entryPath(velocity, ba, intersect, pos, lr, head):

    pi = np.pi
    # print(lr)
    '''Entry Into Bezier Curve'''
    # print(f'BANK ANGLE: {ba}')
    if ba[0] > 0:
        tr = (velocity*1.94384)**2/(11.26*math.tan(np.deg2rad(ba[0])))
    else:
        tr = (velocity*1.94384)**2/(11.26*math.tan(np.deg2rad(ba)))
    tr*=0.3048
    # tr = velocity**2 / (11.26*math.tan(ba[0]))
    # tr*=0.3048
    # print('TR', tr)
    # h = pos[0] + lr*tr*np.cos(head-180)
    # k = pos[1] - lr*tr*np.sin(head-180)

    h = pos[0] + lr*tr
    k = pos[1] 

   
    x_entry = np.linspace(pos[0], intersect[0], 200)
    y_entry = k+np.sqrt(tr**2 - (x_entry-h)**2)

    if math.isnan(y_entry[0]) or math.isnan(x_entry[0]):
        # print('NANS PRESENT')
        y_entry[0] = pos[1]
        x_entry[0] = pos[0]
    if math.isnan(y_entry[-1]) or math.isnan(x_entry[-1]):
        # print('NANS PRESENT')
        y_entry[-1] = intersect[1]
        x_entry[-1] = intersect[0]
            # print(x_entry)
            # print(y_entry)]

    x = np.array([x_entry])
    y = np.array([y_entry])
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_entry = list(x[mask])
    y_entry = list(y[mask])
    
    # y_entry[-1] = intersect[1]
    ar, ad = central_angle([h, k], [x_entry[0], y_entry[0]], [x_entry[-1], y_entry[-1]])
    entryLength = tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
    entryTOA = entryLength/velocity
    # print('ENTRY ToA: ', entryTOA)

    return x_entry, y_entry, ar, entryLength, entryTOA, h, k

def exitPath(velocity, t_exit, ba, intersect, nodes, lr):
    pi = np.pi
    lr = -1
    vk = velocity*1.94384
    '''Entry Into Bezier Curve'''
    # tr = 111.6**2/(11.26*math.tan(ba))
    tr = vk**2/(11.26*math.tan(ba))
    tr*=0.3048

    # Bx = lambda t: nodes[1][0] + (nodes[0][0] - nodes[1][0]) * (1 - t)**2 + (nodes[2][0] - nodes[1][0]) * t**2
    # By = lambda t: nodes[1][1] + (nodes[0][1] - nodes[1][1]) * (1 - t)**2 + (nodes[2][1] - nodes[1][1]) * t**2

    # bezHead = np.arctan2(By(t_exit)-By(t_exit-0.01), Bx(t_exit)-Bx(t_exit-0.01))
    # print(nodes)
    # h, k = intersect[0] - tr*math.cos(bezHead), intersect[1]+tr*math.sin(bezHead)
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    bezHead = np.arctan2(By(t_exit)-By(t_exit-0.01), Bx(t_exit)-Bx(t_exit-0.01))

    nom_y = np.linspace(0, 5000, 200)
    nom_x = [0 for i in nom_y]
    p = []

    angle_to_rotate = 90-np.rad2deg(bezHead)
    for i in range(0, 200):
        p.append([nom_x[i], nom_y[i]])
    p = np.array(p) - np.array([Bx(t_exit), By(t_exit)])  # Shift to origin
   
    rotated_nodes = rotate_bez(nodes, angle_to_rotate, np.array([0, 0]))
    Bx = lambda t: rotated_nodes[0][1] + (rotated_nodes[0][0] - rotated_nodes[0][1]) * (1 - t)**2 + (rotated_nodes[0][2] - rotated_nodes[0][1]) * t**2
    By = lambda t: rotated_nodes[1][1] + (rotated_nodes[1][0] - rotated_nodes[1][1]) * (1 - t)**2 + (rotated_nodes[1][2] - rotated_nodes[1][1]) * t**2

    h = lambda t: Bx(t) + lr*tr#*math.cos(bezHead)
    k = lambda t: By(t) #+ tr*math.sin(bezHead)

    path = manual_bez(P0 = [rotated_nodes[0][0], rotated_nodes[1][0]],
                      P1 = [rotated_nodes[0][1], rotated_nodes[1][1]],
                      P2 = [rotated_nodes[0][2], rotated_nodes[1][2]], 
                      points = 200)
    # print(h, k)
    y_exit = [z for z in np.linspace(By(t_exit), intersect, 200)]
    if lr == 1:
        x_exit = [h(t_exit) - np.sqrt(tr**2 - (y_y - k(t_exit))**2) for y_y in y_exit]
    else:
        x_exit = [h(t_exit) + np.sqrt(tr**2 - (y_y - k(t_exit))**2) for y_y in y_exit]

    rp = rotate_bez(np.array([x_exit, y_exit]), -angle_to_rotate, np.array([0, 0]))
    # if math.isnan(y_exit[0]) or math.isnan(x_exit[0]):
    #     x_entry[0] = 
    # y_exit[0] = 650
    # y_exit[0] = k
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    x_exit = rp[0]
    y_exit = rp[1]
    a = np.arctan2(y_exit[100]-y_exit[99], x_exit[100]-x_exit[99])
    h_exit =  (x_exit[100]+x_exit[99])/2 - tr*np.cos(a)
    k_exit =  (y_exit[100]+y_exit[99])/2 + tr*np.sin(a)
    x = np.linspace(Bx(t_exit), 0, 200)
    y = k_exit-np.sqrt(tr**2 - (x-h_exit)**2)
    # y = np.linspace(By(t_exit), k_exit, 200)
    # x = h_exit + np.sqrt(tr**2 - (y-k_exit)**2)
    # x = np.array([x_exit])
    # y = np.array([y_exit])
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_exit = list(x[mask])
    y_exit = list(y[mask])
    # h_exit = np.mean(x_exit)  # Approximate center x
    # k_exit = np.mean(y_exit)  # Approximate center y

    # print(x_exit)
    # # h = lambda t: Bx(t) + lr * tr * np.cos(bezHead)  
    # # k = lambda t: By(t) + tr * np.sin(bezHead)
    # rot_center = rotate_bez(np.array([h(t_exit), k(t_exit)]), -90, np.array([[Bx(t_exit), By(t_exit)]]))
    # print(rot_center)
    ar, ad = central_angle(center = [h_exit, k_exit], point1=[x_exit[-1], y_exit[-1]], point2=[x_exit[0], y_exit[1]])

    exitLength = tr*ar

    exitTOA = exitLength/velocity
    
    print('EXIT ToA:', exitTOA)
    # y_exit[0] = k
    # print(y_exit[0])
    return x_exit, y_exit, ar, exitLength, exitTOA,h_exit, k_exit#rot_center[0], rot_center[1]


def interp_wps(waypts, n, v):
    interps = np.zeros((len(waypts)-1, 2, n))
    dists = []
    travel = []

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


def getTravelTime(interps, gi, travel, v):
    new_time = []
    total = 0
    # print(travel)
    for i in range(0, len(travel)):
        # print(new_time)
        if gi[0] > i:
            new_time.append(0)
            total+= 0
        elif gi[0] == i and i+1 < len(interps):
            dx = interps[i+1][0][0] - interps[i][0][gi[1]]
            dy = interps[i+1][1][0] - interps[i][1][gi[1]]
            d = np.hypot(dx, dy)
            new_time.append(d/v)
            total+=d/v
        elif gi[0] == i and i+1>= len(interps):
            dx = interps[i][0][-1] - interps[i][0][gi[1]]
            dy = interps[i][1][-1] - interps[i][1][gi[1]]
            d = np.hypot(dx, dy)
            new_time.append(d/v)
            total+=d/v
        elif gi[0] < i:
            new_time.append(travel[i])
            total+=travel[i]
    # print(f'STAR PATH Travel Time: {total}')
    return new_time, total

def getTravelXY(wps, point, gi, travel, v):
    new_time = []
    total = 0
    
    for i in range(len(travel)):
        if gi > i:
            new_time.append(0)
            total+=0
        elif gi == i and i+1 < len(wps):
           
            d = np.hypot(wps[i+1][0] - point[0] , wps[i+1][1] - point[1])
            new_time.append(d/v)
            total+=d/v
        elif gi == i and i+1>= len(wps):

            d = np.hypot(point[0] - wps[i][0], point[1] - wps[i][1])
            new_time.append(d/v)
            total+=d/v
        elif gi < i:
            new_time.append(travel[i])
            total+=travel[i]
    return new_time, total

def total_travel(waypts, v):
    dists = []
    travel = []

    total = 0
    for i in range(0, len(waypts)-1):
        d = np.hypot(waypts[i+1][0] - waypts[i][0], waypts[i+1][1] - waypts[i][1])
        dists.append(d)
        travel.append(d/v)
        total+= d/v

    print(f'Total STAR Route Travel Time: {total}')
    return dists, travel, total

def get_guess_point(travel_times, dists, required_toa, required_len, index):
    # print(travel_times, len(travel_times))
    total_time = 0
    total_dist = 0
    for i in range(len(travel_times)):
        if i > index:
            total_time+=travel_times[i]
            total_dist+=dists[i]

    segment_toa = required_toa-total_time
    # print(f'SEGMENT TOA: {segment_toa}')
    segment_dist = required_len-total_dist
    # print(f'REQUIRED DIST: {required_len}, REMAINING: {total_dist}, SEGMETNT DIST: {segment_dist}, CURRENT SEGMENT LEN: {dists[index]}')
    if segment_dist > 0:
        gp = 1 - (segment_dist/dists[index])
    else:
        gp = 0
    # print(f'GUESS: {gp}, SEGMENT DIST/TOTAL SEGEMTN: {segment_dist/dists[index]}')

    return gp

def rotate_points(points, angle_deg, origin):
    """Rotate a set of points around a given origin by angle_deg."""
    angle_rad = np.radians(angle_deg)
    R = np.array([[np.cos(angle_rad), -np.sin(angle_rad)], 
                  [np.sin(angle_rad), np.cos(angle_rad)]])

    points = np.array(points)  # Ensure points are NumPy array
    origin = np.array(origin)

    # Center points at origin, rotate, and shift back
    return (R @ (points.T - origin[:, None])).T + origin


def rotate_bez(points, angle, origin):
    """Rotate a set of points around a given origin by angle_deg."""
    angle = np.deg2rad(angle)
    R = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]
    ])

    # Ensure origin has the correct shape (2, 1) for broadcasting
    origin = origin.reshape(2, 1)

    # Apply rotation: shift to origin, rotate, then shift back
    rotated_points = R @ (points - origin) + origin  # Using @ for matrix multiplication
    return rotated_points

def rotate_interpolated_points(interps, angle_deg, origin):
    """Rotate interpolated waypoints stored in a 3D NumPy array."""
    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    
    # Reshape interps for rotation
    num_segments, _, num_points = interps.shape  # (len(waypts)-1, 2, n)
    interps_reshaped = interps.reshape(-1, 2)  # Flatten to (num_segments * n, 2)

    # Apply rotation
    rotated_points = np.dot(interps_reshaped - origin, R.T) + origin

    # Reshape back to original format
    rotated_interps = rotated_points.reshape(num_segments, 2, num_points)

    return rotated_interps

def get_line(p1, p2, t):
    if t<0:
        t=0
    elif t>1:
        t=1
    # m = (p2[1]-p1[1])/(p2[0]-p1[0])
    x = p1[0] + (p2[0]-p1[0])*t
    y = p1[1] + (p2[1]-p1[1])*t
    return [x,y]

def get_line_p1(p1, p2, t, h):
    if t>1:
        t=1
    if p2[1] <= p1[1]:
        # p2[0]+=1000
        p2[1] +=100*np.sin(h)
        if p2[0]>p1[0]:
            p2[0]+=100*np.cos(h)
        else:
            p2[0]-=100*np.cos(h)
    # m = (p2[1]-p1[1])/(p2[0]-p1[0])
    x = p1[0] + (p2[0]-p1[0])*t
    y = p1[1] + (p2[1]-p1[1])*t
    return [x,y]

def toCallOutside(velocity, nodes1, nodes2, t1, t2, lr, p1, p2):

    vel_knots = velocity * 1.94384
    turn_radius = vel_knots**2 / (11.26*np.tan(np.deg2rad(73))) #Minimum turn radius
    valid = False
    v1 = False
    v2 = False

    if nodes1[0][1] != nodes2[0][1]:
        m_p21 = (nodes2[1][1]-nodes1[1][2])/(nodes2[0][1]-nodes1[0][2])
        b2 = nodes2[1][1] - m_p21*nodes2[0][1]

        lines_coeffs2 = [m_p21, b2]
    else:
        m_p21 = 0
        b2 = nodes2[1][1] - m_p21*nodes2[0][1]

        lines_coeffs2 = [m_p21, b2]

    b1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
    b2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
    ly = np.linspace(0, 2400, 200)
    lx = 0*ly
    # 
    
    # plt.text(-750, 1500, 'Bezier 2', fontweight = 'bold', fontsize = 13)
    # plt.text(-750, 750, 'Bezier 1', fontweight = 'bold', fontsize = 13)
    # plt.plot(b1[0], b1[1], color = 'blue' , linewidth = 2)
    # plt.plot(b2[0], b2[1], color = 'orange' , linewidth = 2)
    # plt.plot(lx, ly, linewidth = 2, color = 'purple', linestyle = '--', alpha = 0.5, label = 'Nominal Path')
    # plt.scatter(p1[0], p2[1], marker = 'o', color = 'red', s = 50, label = 'Collision Point')
    # plt.scatter(nodes1[0], nodes1[1], color = 'blue', label = 'Bezier 1 Control Points')
    # plt.scatter(nodes2[0], nodes2[1], color = 'orange', label = 'Bezier 2 Control Points')
    # plt.scatter(p1[0], p1[1], marker = '^', color = 'black', label =  r'$a_0$', s = 75)
    # plt.scatter(p2[0], p2[1], marker = '>', color = 'green', label = r'$a_1$', s = 75)
    # plt.arrow(p1[0], p1[1]+50, 75*np.cos(np.pi/2), 75*np.sin(np.pi/2), width=25,
    #             length_includes_head=False,
    #             head_width=50,
    #             head_length=50,
    #             head_starts_at_zero=False,
    #             facecolor='black',
    #             edgecolor='black',)
    #             # label = f'Heading of {np.rad2deg(np.pi/2):0.2f}')
    # plt.arrow(p2[0]+50, p2[1], 75*np.cos(0*np.pi/2), 75*np.sin(0*np.pi/2), width=25,
    #             length_includes_head=False,
    #             head_width=50,
    #             head_length=50,
    #             head_starts_at_zero=False,
    #             facecolor='green',
    #             edgecolor='green',)
    #             # label = f'Heading of {np.rad2deg(0*np.pi/2):0.2f}')
    # plt.xlim((-1500, 1500))
    # plt.grid('on')
    # # plt.axis('equal')
    # plt.xlabel('X (m)')
    # plt.ylabel('Y (m)')
    # plt.legend(loc = 'center right')
    # plt.show()
    while not valid:
        optim_sol1, curv1 = solve_optim1(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=t1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs2)
        if nodes1[0][2]!=optim_sol1[0]:
            m_p12 = (nodes1[1][2]-optim_sol1[1])/(nodes1[0][2]-optim_sol1[0])
            b = optim_sol1[1] - m_p12*optim_sol1[0]
            lines_coeffs = [m_p12, b]
        else:
            m_p12 = 0
            b = optim_sol1[1] - m_p12*optim_sol1[0]
            lines_coeffs = [m_p12, b]
        # print('M12',m_p12, 'B', b)

        x_slope = [i for i in np.linspace(optim_sol1[0]-10, lr*245, 50)]
        y_slope = [m_p12*x + b for x in x_slope]

        # print(lines_coeffs)

        optim_sol2, curv2 = solve_optim2(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=t2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr, line = lines_coeffs)
        
        if nodes1[0][2] != optim_sol2[0]:
            m_p21 = (nodes1[1][2]-optim_sol2[1])/(nodes1[0][2]-optim_sol2[0])
            b2 = optim_sol2[1] - m_p21*optim_sol2[0]
            lines_coeffs2 = [m_p21, b2]
        else:
            m_p21 = 0
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
        
        nodes1 = [np.array([nodes1[0][0], optim_sol1[0], nodes1[0][2]]).flatten(),np.array([nodes1[1][0], optim_sol1[1], nodes1[1][2]]).flatten()]
        nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]

        optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
        optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)

        if np.abs((optim1_length+optim2_length)/velocity - (t1+t2)) <= 2:
            valid = True
        elif valid == False:
            nodes1[1][1] -= 1
            nodes1[0][1] -= 1
            nodes2[1][1] += 1
            nodes2[0][1] += 1
        print((optim1_length+optim2_length)/velocity)
        x = [0 for i in range(200)]
        y = [i*20 for i in range(len(x))]


    plt.plot(optimal_bez1[0], optimal_bez1[1], color = 'black')
    plt.text(-750, 750, 'Curve 1', fontweight = 'bold')
    plt.plot(optimal_bez2[0], optimal_bez2[1], color = 'black')
    plt.text(-750, 1750, 'Curve 2', fontweight = 'bold')
    plt.scatter(nodes1[0], nodes1[1], color = 'blue', label = 'Bezier 1 Control Points')
    plt.scatter(nodes2[0], nodes2[1], color = 'orange', label = 'Bezier 2 Control Points')
    plt.axis('equal')
    plt.grid()
    plt.scatter(p1[0], p1[1], color = 'black', marker = '^', label = r'$a_0$', zorder = 200, s = 75)
    plt.scatter(p2[0], p2[1], marker = '>', color = 'green', label = r'$a_1$', s = 75)
    plt.arrow(p2[0]+50, p2[1], 75*np.cos(0*np.pi/2), 75*np.sin(0*np.pi/2), width=25,
                length_includes_head=False,
                head_width=50,
                head_length=50,
                head_starts_at_zero=False,
                facecolor='green',
                edgecolor='green',)
                # label = f'Heading of {np.rad2deg(0*np.pi/2):0.2f}')
    plt.arrow(p1[0], p1[1]+50, 75*np.cos(np.pi/2), 75*np.sin(np.pi/2), width=25,
                length_includes_head=False,
                head_width=50,
                head_length=50,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                zorder = 100)
    # plt.xlim((-1500, 1500))
    # plt.grid('on')
    # # plt.axis('equal')
    # plt.xlabel('X (m)')
    # plt.ylabel('Y (m)')
    # plt.legend(loc = 'center right')
    # plt.show()
    plt.plot(x, y, linestyle = '--', linewidth = 2, color = 'blue', label = 'Nominal Path')
    plt.plot([nodes1[0][1], nodes2[0][1]], [nodes1[1][1], nodes2[1][1]], color = 'purple', linewidth = 2, label = 'G1 Continuity Line')
    
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.legend(loc = 'upper right')
    plt.show()

    return optimal_bez1, optimal_bez2, nodes1, nodes2

def EntryExitOutside(nodes1, nodes2, pos, velocity, lr, p2):
    start = time.time()
    vel_knots = velocity*1.94384

    optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
    optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)
    
    # if id!= 0:
    # if lr == 1: 
    #     ba, t_entry = solve_optimEntry(np.deg2rad(10), np.deg2rad(73), np.deg2rad(5), nodes1, vel_knots, pos, lr)
    # else:
    ba, t_entry = solve_optimEntry(np.deg2rad(30), np.deg2rad(73), np.deg2rad(5), nodes1, vel_knots, pos, lr, np.pi/2, False)

    # else:
    #     ba, t_entry = solve_optimEntry(np.deg2rad(15), np.deg2rad(73), np.deg2rad(15), nodes1, vel_knots, pos, lr)

    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry+0.0025)
    
    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))

    # solve_optimEntry(guess, max_bank,  min_bank, nodes, velocity, pos, lr, h, show):
    x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en], pos, lr, np.pi/2)
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
    print('Entry BANK', np.rad2deg(ba), 'EXIT T:', t_entry)


    t_start2 = 0.35678391959798994

    # if lr == -1:
    #     baEx, exit_t = solve_optimExit([np.deg2rad(30), 0.45], np.deg2rad(60), np.deg2rad(15), t_start2, nodes2, velocity, lr)
    # else:
    baEx, exit_t, yval = solve_optimExit([np.deg2rad(22.5), 0.75], np.deg2rad(73), np.deg2rad(5), t_start2, nodes2, velocity, lr)

    # else:
    #     baEx, exit_t = solve_optimExit([np.deg2rad(60), 0.55], np.deg2rad(60), np.deg2rad(20), t_start2, nodes2, velocity, lr)
    print('EXIT BANK', np.rad2deg(baEx), 'EXIT T:', exit_t, 'Velocity', velocity)

    x_int_ex, y_int_ex = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t)
    
    x_int_ex2, y_int_ex2 = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t-0.0025)
    
    req_ex = np.rad2deg(np.arctan2(y_int_ex-y_int_ex2,x_int_ex-x_int_ex2))

    x_exit, y_exit, central_angle_ex, exitLength, exitTOA, h_ex, k_ex = exitPath(velocity=velocity, t_exit=exit_t,
                                                                                    ba = baEx, intersect=yval,
                                                                                    nodes = nodes2, lr = lr)
    # x_exit[0] = 750
    # y_exit[0] = 2350
    print(y_exit[0])
    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], 200, exit_t+.025)
    
    pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_entry)
    pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]], P2=[nodes2[0][2],nodes2[1][2]], t=exit_t)
    total_toa = entryTOA+(pb1_length/velocity)+(pb2_length/velocity)+exitTOA
    print('TOTAL TOA:', total_toa)

    paths = [partial_bez1, partial_bez2]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    plt.figure()
    plt.text(-750, 1500, 'Bezier 2', fontweight = 'bold', fontsize = 13)
    plt.text(-750, 750, 'Bezier 1', fontweight = 'bold', fontsize = 13)
    plt.plot(partial_bez1[0], partial_bez1[1], color = 'blue', linewidth = 2)
    plt.plot(partial_bez2[0], partial_bez2[1], color = 'orange', linewidth = 2)
    plt.plot(x_entry, y_entry, color = 'green', linewidth = 2, label = 'Entry/Exit Path')
    plt.plot(x_exit, y_exit, color = 'green', linewidth = 2)
    plt.plot([0, 0], [-100, 3500], linestyle = '--', linewidth = 2, color = 'blue', label = 'Nominal Path')
    plt.scatter(pos[0], pos[1], color = 'black', marker = '^', label = r'$a_0$', zorder = 200, s = 75)
    plt.scatter(p2[0], p2[1], marker = '>', color = 'green', label = r'$a_1$', s = 75)
    plt.arrow(p2[0]+100, p2[1], 75*np.cos(0*np.pi/2), 75*np.sin(0*np.pi/2), width=25,
                length_includes_head=False,
                head_width=50,
                head_length=50,
                head_starts_at_zero=False,
                facecolor='green',
                edgecolor='green',)
                # label = f'Heading of {np.rad2deg(0*np.pi/2):0.2f}')
    plt.arrow(pos[0], pos[1]+100, 75*np.cos(np.pi/2), 75*np.sin(np.pi/2), width=25,
                length_includes_head=False,
                head_width=50,
                head_length=50,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                zorder = 100)
    plt.axis('equal')
    plt.grid('on')
    # plt.scatter(x_wpts, y_wpts)
    plt.scatter(h_c, k_c, color = 'purple', label = 'Exit Center')
    plt.scatter(h_ex, k_ex, color = 'red', label = 'Entry Center')
    plt.legend(loc = 'upper right')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.show()

    # plt.plot(optimal_bez1[0], optimal_bez1[1], color = 'black')
    # plt.text(-750, 750, 'Curve 1', fontweight = 'bold')
    # plt.plot(optimal_bez2[0], optimal_bez2[1], color = 'black')
    # plt.text(-750, 1750, 'Curve 2', fontweight = 'bold')
    # plt.scatter(nodes1[0], nodes1[1], color = 'blue', label = 'Bezier 1 Control Points')
    # plt.scatter(nodes2[0], nodes2[1], color = 'orange', label = 'Bezier 2 Control Points')
    # plt.axis('equal')
    # plt.grid()
    # plt.scatter(p1[0], p1[1], color = 'black', marker = '^', label = r'$a_0$', zorder = 200, s = 75)
    # plt.scatter(p2[0], p2[1], marker = '>', color = 'green', label = r'$a_1$', s = 75)
    # plt.arrow(p2[0]+50, p2[1], 75*np.cos(0*np.pi/2), 75*np.sin(0*np.pi/2), width=25,
    #             length_includes_head=False,
    #             head_width=50,
    #             head_length=50,
    #             head_starts_at_zero=False,
    #             facecolor='green',
    #             edgecolor='green',)
    #             # label = f'Heading of {np.rad2deg(0*np.pi/2):0.2f}')
    # plt.arrow(p1[0], p1[1]+50, 75*np.cos(np.pi/2), 75*np.sin(np.pi/2), width=25,
    #             length_includes_head=False,
    #             head_width=50,
    #             head_length=50,
    #             head_starts_at_zero=False,
    #             facecolor='black',
    #             edgecolor='black',
    #             zorder = 100)
    # # plt.xlim((-1500, 1500))
    # # plt.grid('on')
    # # # plt.axis('equal')
    # # plt.xlabel('X (m)')
    # # plt.ylabel('Y (m)')
    # # plt.legend(loc = 'center right')
    # # plt.show()
    # plt.plot(x, y, linestyle = '--', linewidth = 2, color = 'blue', label = 'Nominal Path')
    # plt.plot([nodes1[0][1], nodes2[0][1]], [nodes1[1][1], nodes2[1][1]], color = 'purple', linewidth = 2, label = 'G1 Continuity Line')
    
    # plt.xlabel('X (m)')
    # plt.ylabel('Y (m)')
    # plt.legend(loc = 'upper right')
    # plt.show()



    end = time.time()
    print('TOTAL INTERCEPTION PATH GENERATION TIME:', end-start)

    # DubinsData = {
    #     'Entry':{
    #         'path_type': 'Entry',
    #         'bez_intx': x_int_en,
    #         'bez_inty': y_int_en,
    #         'nom_intx': pos[0],
    #         'nom_inty': pos[1],
    #         'bez_t': t_entry,
    #         'intersect_heading': req_ent,
    #         'h': h_c,
    #         'k': k_c,
    #         'bank_angle': np.rad2deg(ba[0]),
    #         'tr': (velocity*1.9438445)**2/(11.26*math.tan(np.deg2rad(ba))),
    #         'timeStamp': timeStamp,
    #         'ACID': id,
    #         'ExpNum': expnum,
    #         'Category': 'Fleet Aircraft',
    #         'ExpType': exptype,
    #         'ax_spacing': spacing
    #     },
    #     'Exit':{
    #         'path_type': 'Exit',
    #         'bez_intx': x_int_ex,
    #         'bez_inty': y_int_ex,
    #         'nom_intx': x_exit[0],
    #         'nom_inty': y_exit[0],
    #         'bez_t': exit_t,
    #         'intersect_heading': req_ex,
    #         'h': h_ex,
    #         'k': k_ex,
    #         'bank_angle': np.rad2deg(baEx),
    #         'tr': (velocity*1.9438445)**2/(11.26*math.tan(np.deg2rad(ba))),
    #         'timeStamp': timeStamp,
    #         'ACID': id,
    #         'ExpNum': expnum,
    #         'Category': 'Fleet Aircraft',
    #         'ExpType': exptype,
    #         'ax_spacing': spacing
    #     }
    # }

    return [x_entry, y_entry, req_ent, np.rad2deg(ba)], [x_exit, y_exit, req_ex, exit_t, np.rad2deg(baEx)], total_toa, x_wpts, y_wpts