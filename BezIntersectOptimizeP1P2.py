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
            # {'type': 'ineq', 'fun': lambda x: x[0]-P2[0]},
            {'type': 'ineq', 'fun': lambda x: P2[1]-x[1]-1000},
            # {'type': 'ineq', 'fun': lambda x: np.sqrt((x[0] - P2[0])**2 + (x[1] - P2[1])**2) - 500},
            # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
            {'type': 'eq', 'fun': lambda x: (line[0]*x[0] + line[1])-x[1]}
            ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-6, constraints=cons)
 
    return val.x , curvature(P0,val.x,P2)

# def solve_optim12(P0, target_toa, guess, target_heading, velocity, turn_radius, lr, line, wps):#turn_radius,
#     # def path_cost(guess,target_toa, velocity):
#     #     # P1 = [guess[0], guess[1]]
#     #     # P2 = get_line(wps[0], wps[1], guess[2])
#     #     # diff = toa_diff([guess[0], guess[1]], P0, get_line(wps[0], wps[1], guess[2]), 1, target_toa, velocity)
#     #     # print(f'INPUT GUESS: {guess[2]}')
#     #     return toa_diff([guess[0], guess[1]], P0, get_line(wps[0], wps[1], guess[2]), 1, target_toa, velocity)
#     def path_cost(guess, target_toa, velocity):
#         weight_xy = 15  # Increase influence of guess[0] and guess[1]
#         weight_z = 0.1    # Keep guess[2] changes smaller
#         return weight_xy * toa_diff([guess[0], guess[1]], P0, get_line(wps[0], wps[1], guess[2]), 1, target_toa, velocity) + weight_z * np.abs(guess[2])

#     cons = (
#             {'type': 'ineq', 'fun': lambda x: curvature(P0,[x[0], x[1]],get_line(wps[0], wps[1], x[2])) - turn_radius},
#             {'type': 'ineq', 'fun': lambda x: curvature(P0,[x[0], x[1]],get_line(wps[0], wps[1], x[2]))},
#             # {'type': 'ineq', 'fun': lambda x: wps[1][1]-get_line(wps[0], wps[1], x[2])[1]},
#             # {'type': 'ineq', 'fun': lambda x: get_line(wps[0], wps[1], x[2])[1]-x[1]},
#             {'type': 'ineq', 'fun': lambda x: 1-np.abs(x[2])},
#             {'type': 'ineq', 'fun': lambda x: x[2]-0},
#             {'type': 'ineq', 'fun': lambda x: get_line(wps[0], wps[1], x[2])[1]-x[1]-100},
#             # {'type': 'ineq', 'fun': lambda x: np.sqrt((x[0] - P2[0])**2 + (x[1] - P2[1])**2) - 500},
#             # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
#             {'type': 'eq', 'fun': lambda x: (line[0]*x[0] + line[1])-x[1]}
#             ) 

#     val = minimize(path_cost,guess, (target_toa, velocity), method='SLSQP', tol=1E-3, constraints=cons)

 
#     return [val.x[0], val.x[1]], val.x[2] , curvature(P0,[val.x[0], val.x[1]],get_line(wps[0], wps[1], val.x[2]))

def solve_optim12(P0, target_toa, guess, target_heading, velocity, turn_radius, lr, line, wps):#turn_radius,
    def path_cost(guess,target_toa, velocity):
        weight_p1 = 7.5  # Increase influence of guess[0] and guess[1]
        weight_p2 = 1   # Keep guess[2] changes smaller
        P2 = get_line(wps[0], wps[1], guess[1])
        return weight_p1*toa_diff( get_line_p1(wps[0], P2 , guess[0], target_heading),
                                   P0, P2, 1, target_toa, velocity) + weight_p2 * np.abs(guess[1])

   
    
    cons = (
            {'type': 'ineq', 'fun': lambda x: curvature(P0, get_line_p1(wps[0], get_line(wps[0], wps[1], x[1]), x[0], target_heading),get_line(wps[0], wps[1], x[1])) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0, get_line_p1(wps[0], get_line(wps[0], wps[1], x[1]), x[0], target_heading),get_line(wps[0], wps[1], x[1]))},
            # {'type': 'ineq', 'fun': lambda x: wps[1][1]-get_line(wps[0], wps[1], x[2])[1]},
            # {'type': 'ineq', 'fun': lambda x: get_line(wps[0], wps[1], x[2])[1]-x[1]},
            {'type': 'ineq', 'fun': lambda x: 1-np.abs(x[1])},
            {'type': 'ineq', 'fun': lambda x: x[1]-0},
            {'type': 'ineq', 'fun': lambda x: x[1]-x[0]},
            {'type': 'ineq', 'fun': lambda x: get_line(wps[0], wps[1], x[1])[1]- get_line_p1(wps[0],get_line(wps[0], wps[1], x[1]), x[0], target_heading)[1] -100},
            # {'type': 'ineq', 'fun': lambda x: x[1]-x[0]-0.001},
            # {'type': 'ineq', 'fun': lambda x: np.sqrt((x[0] - P2[0])**2 + (x[1] - P2[1])**2) - 500},
            # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
            {'type': 'eq', 'fun': lambda x: (line[0]*get_line_p1(wps[0], get_line(wps[0], wps[1], x[1]), x[0], target_heading)[0] + line[1])-get_line_p1(wps[0], get_line(wps[0], wps[1], x[1]), x[0], target_heading)[1]}
            ) 

    val = minimize(path_cost,guess, (target_toa, velocity), method='SLSQP', tol=1E-3, constraints=cons)

 
    return val.x[0], val.x[1] , curvature(P0,get_line_p1(wps[0], get_line(wps[0], wps[1], val.x[1]), val.x[0], target_heading),get_line(wps[0], wps[1], val.x[1]))

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

    # tr = 111.6**2 / (11.26*math.tan(ba))
    tr = 106.6**2 / (11.26*math.tan(ba))
    tr*=0.3048
    h = lambda t: Bx(t) - lr*tr*math.cos(bezHead)
    k = lambda t: By(t) + tr*math.sin(bezHead)

    circle_eq = lambda y: (0-h(t_guess))**2+(y-k(t_guess))**2 - tr**2
    y_guess = 610
    S = fsolve(circle_eq, y_guess) #gives y intersection
    t_final = .5

    y_l = [i for i in np.linspace(274, S, 200)]
    x = [0 for i in y_l]

    # print(S)
    for i in S:
        if i > 0: 
            y = [i for i in np.linspace(S, By(t_guess))]
            if lr == 1:
                x_l = [h(t_guess) - np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            else:
                x_l = [h(t_guess) + np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            if x_l[0]<=0.01: #==0:# <= 0.01:
                # x_l = [i for i in np.linspace(750, Bx(t_guess))]
                # y = [k(t_guess)-np.sqrt(tr**2 - (x-h(t_guess))**2) for x in x_l]
                int_angle = np.arctan2(y[0]-y[1], x_l[0]-x_l[1])
                diff = np.abs((np.pi/2) - int_angle)
                guess_deg = np.rad2deg(ba)
                # plt.title(f'{guess_deg} {t_guess}')
                # plt.plot(x, y_l, linestyle = 'dashed')
                # plt.plot(path[0], path[1])
                # plt.plot(x_l, y)
                # plt.scatter(h(t_guess), k(t_guess))
                # plt.axis('equal')
                # plt.show()
                if diff < mindiff:
                    mindiff = diff
                    # print('BLARGY', diff)
                    t_final = t_guess
                

    return [mindiff, t_final]

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
    '''Entry Into Bezier Curve'''
    # tr = 111.6**2/(11.26*math.tan(ba))
    tr = 106.6**2/(11.26*math.tan(ba))
    tr*=0.3048

    Bx = lambda t: nodes[1][0] + (nodes[0][0] - nodes[1][0]) * (1 - t)**2 + (nodes[2][0] - nodes[1][0]) * t**2
    By = lambda t: nodes[1][1] + (nodes[0][1] - nodes[1][1]) * (1 - t)**2 + (nodes[2][1] - nodes[1][1]) * t**2

    bezHead = np.arctan2(By(t_exit)-By(t_exit-0.01), Bx(t_exit)-Bx(t_exit-0.01))

    h, k = intersect[0] - tr*math.cos(bezHead), intersect[1]+tr*math.sin(bezHead)
    print(h, k)
    x_exit = [i for i in np.linspace(0, intersect[0], 200)]
    y_exit = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_exit]
    # y_exit[0] = 650
    # y_exit[0] = k
    ar, ad = central_angle(center = [h, k], point1=[x_exit[-1], y_exit[-1]], point2=[x_exit[0], y_exit[1]])

    exitLength = tr*ar

    exitTOA = exitLength/velocity
    
    print('EXIT ToA:', exitTOA)
    y_exit[0] = k
    print(y_exit[0])
    return x_exit, y_exit, ar, exitLength, exitTOA, h, k


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
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    
    # Shift points to origin, rotate, and shift back
    rotated_points = np.dot(points - origin, R.T) + origin
    return rotated_points

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

def GenerateTrajectoryOutside(ac_pos, target, wp_list, h, velocity, kcmh, gate, plot, plans, t_it):
    v=velocity
    t_og = target
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    print(h, gate)
    if gate == 'II':
        # diff = h-90
        # heading = 360-diff
        # # heading = h+180
        # head = np.deg2rad(heading)
        # angle_to_rotate = (90-heading)
        g_index = 5
    elif gate == 'III':
        # heading=h+90
        # head = np.deg2rad(heading)
        # angle_to_rotate = (90-heading)
        g_index = 2
    elif gate == 'IV':
        # heading = h-90
        # head = np.deg2rad(heading)
        # angle_to_rotate = (90-heading)
        g_index = 2
    elif gate == 'I':
        # heading=270-h
        # head = np.deg2rad(heading+180)
        # angle_to_rotate = (90-heading+180)
        g_index = 1
    diff = h-90
    heading = 360-diff
    # heading = h+180
    head = np.deg2rad(heading)
    angle_to_rotate = (90-heading)
    rotated_heading = np.deg2rad(h+angle_to_rotate)

    
    wps = LongLat_To_WSG84_Meters(wp_list, kcmh) #Gives waypoints in global frame
    ac_met = __WSG84_To_Meters_Single(ac_pos, kcmh) #Global aircraft position
    wps_local = LongLat_To_WSG84_Meters(wp_list, ac_pos)

    print(wps)
    x = [i[0] for i in wps]
    y = [i[1] for i in wps]

    

    origin = np.array([ac_met])
    

    rotated_wps = rotate_points(wps, angle_to_rotate, origin)
    

    x_rot = [i[0] for i in rotated_wps]
    y_rot = [i[1] for i in rotated_wps]
    
    
    dists_local, toa_local, total_local = total_travel(wps_local, v)
    # print(f'REGULAR DISTANCE TO FIRST WP: {d_reg}, ROTATED: {d_rot}')

    if rotated_wps[1][0] < ac_met[0]:
        lr = -1
    else:
        lr = 1

    bank = 5
    min = 5
    plan_ahead = 10
    
    path_startX = ac_met[0] + np.cos(head)*v*plan_ahead 
    path_startY = ac_met[1] + np.sin(head)*v*plan_ahead 

    path_startX_local = ac_met[0]
    path_startY_local = ac_met[1] + v*plan_ahead
    # plt.plot(x, y, color = 'blue')
    # plt.plot(x_rot, y_rot, 'orange')
    # plt.scatter(origin[0][0], origin[0][1])
    # plt.plot([origin[0][0], path_startX], [origin[0][1], path_startY], color = 'red')
    # plt.plot([origin[0][0], path_startX_local], [origin[0][1], path_startY_local], color = 'black')
    # plt.show()


    turn_radius = v*1.94384**2/(11.26*math.tan(np.deg2rad(73)))
    turn_radius*=0.3048
    
    # path_total = 0
    valid = False
    gi = [0,0]
    t_t = target-total_local-plan_ahead

    if plans[6]:
        buffer = 0
    else:
        if plans[1]-plans[5] > 0:
            buffer = t_it - (plans[1]-plans[5])
            target+=buffer
        else:
            buffer = t_it + (plans[1]-plans[5])
            target+=buffer

    end_point = get_line(rotated_wps[1], rotated_wps[2], 0)
    end_point_ref = get_line(rotated_wps[1], rotated_wps[2], 0.1)
    new_travel, new_total = getTravelXY(wps, end_point, gi[0], toa_local, v)
    t_h = np.arctan2(end_point_ref[1]- end_point[1],end_point_ref[0]- end_point[0])
    # t_h2 = np.arctan2(rotated_wps[g_index][1]-rotated_wps[g_index-1][1], rotated_wps[g_index][0]-rotated_wps[g_index-1][0])

    t_t = target-new_total-plan_ahead
    end_point = get_line(rotated_wps[1], rotated_wps[2], gi[1])
    entry = False
    entryTOA = plan_ahead
    all_bez = []
    all_nodes = []
    opt = []
    show = False
    path_toa = new_total+plan_ahead
    start_time = time.time()
    
    
    for i in range(len(wps)-1):

        # g = 5000
        gi[1] = 0.5
        p1g = 0.25
        valid = False
        w1 = rotated_wps[gi[0]]
        w2 = rotated_wps[gi[0]+1]
        m = (w2[1] - w1[1])/(w2[0] - w1[0])
        b = w1[1] - m*w1[0]
        coeffs = [m, b]
        t_h = np.arctan2(w2[1]- w1[1],w2[1]- w1[0])
        n = 15
        for j in range(0,n):
            # if not plans[6]:
            #     print(f'OG TARGET IS {t_og}, NEW TARGRET IS {target}, DIFF IS {(path_toa+new_total) - plans[5]}, {path_toa+new_total}, {plans[5]}')

            end_point = get_line(w1, w2, gi[1])
            p1_point = get_line_p1(w1, end_point, p1g, t_h)
            new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            t_t = target-new_total-entryTOA
            
       
            nodes = [np.array([path_startX_local, p1_point[0], end_point[0]]).flatten(),
                np.array([path_startY_local, p1_point[1], end_point[1]]).flatten()]
         
            p1, end, curv1 = solve_optim12(P0 = [nodes[0][0],nodes[1][0]],
                                        target_toa=t_t,
                                        guess = [p1g, gi[1]],
                                        target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs, wps=[w1, w2])
            # print(f'OPTIMIZATION RESULT: {gi[0], end}')
            end_point = get_line(w1, w2, end)
            optim1 = get_line_p1(w1, end_point, p1, t_h)
            new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            t_t = target-new_total-entryTOA
            p1g=p1

            if optim1[0] < ac_met[0]:
                lr = -1
            else:
                lr = 1

            optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [end_point[0], end_point[1]], 50)
            optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],1)
            path_toa = optim_b_length/v
            
            nodes1 = [np.array([nodes[0][0], optim1[0], end_point[0]]).flatten(), np.array([nodes[1][0], optim1[1], end_point[1]])]

            if math.isnan(path_toa):
                    optim1[0]-=1000*np.cos(t_h)
                    optim1[1]-=1000*np.sin(t_h)
                    optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],1)
                    path_toa = optim_b_length/v
            
            if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(15)):
                entryTOA = plan_ahead
                path_toa = optim_b_length/v + entryTOA
                entry = False
                # req_entry = np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])
                
            else:
                entry = True
                bez_ang = np.rad2deg(np.arctan2(optim_b[1][1]-optim_b[1][0], optim_b[0][1]-optim_b[0][0]))
                
                # if bez_ang<0:
                #     bez_ang+=360
                # if bez_ang>=270:
                #     bank = 30
                # else:
                #     bank = 5
                # print(f'BEZ ANGLE: {bez_ang}, GUESS {bank}')
                ba, t_entry = solve_optimEntry([np.deg2rad(bank), 0.005], np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotated_heading, show)
                


                x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [end_point[0], end_point[1]], t_entry)
               
                x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [end_point[0], end_point[1]], t_entry+0.01)
                req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
                
                
                x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotated_heading)
                act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
                

                partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [end_point[0], end_point[1]], 50, t_entry)
                pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],t_entry)
                pb_toa = pb_length/v
                path_toa = entryTOA + pb_toa

            
            '''
            Plots for step by step
            '''
            ang_y = (end_point[0] - 500) *m + b
            # plt.plot([end_point[0],end_point[0] - 5000], [end_point[1], ang_y], color = 'green', linestyle = '--')
            # plt.plot(x_rot, y_rot)
            # plt.scatter(ac_met[0], ac_met[1])
            # plt.plot(x_rot, y_rot)
            # plt.scatter(optim1[0], optim1[1], marker = '^', color = 'red')
            # plt.plot(optim_b[0], optim_b[1], color = 'black')
            # plt.scatter(end_point[0], end_point[1], marker = 's', color = 'black')
            # plt.axis('equal')
            # plt.text(rotated_wps[-1][0]-5000, rotated_wps[-1][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
            # plt.text(rotated_wps[1][0], rotated_wps[1][1], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )
            # # if math.isnan(path_toa):
            # #     print(optim_b[0])
            # #     print(optim_b[1])
            # #     print(nodes1)
            # #     print(pb_length, optim_b_length)
            # #     print(x_entry)
            # #     print(y_entry)
            # if entry:
            #     plt.plot(x_entry, y_entry, color = 'cyan')
            #     plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
            #     plt.text(optim_b[0][25], optim_b[1][25]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
            #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
            #     plt.text(x_entry[0], y_entry[0], f'BANK ANGLE: {np.rad2deg(ba[0]):0.4f}', weight = 'bold')
            # else:
            #     # print(f'Bez Travel Time: {path_toa:0.2f}s')
            #     plt.text(optim_b[0][25], optim_b[1][25], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
            
            
            
            
            # print(new_total + path_toa)
            # print(f'TARGET: {target}, ACTUAL TRAVEL TIME: {path_toa+new_total}')
            # print(np.abs(target - (path_toa+new_total)))
            # print(f'Original Target: {target}, New STAR Travel Time: {new_total}, New Target Time: {t_t}, Path TOA: {path_toa}')
            # print(new_total+path_toa, target)
            # plt.show()
        
                
            if entry and not math.isnan(path_toa):
                all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, [ba, t_entry], entryTOA, [act_ent, req_ent], curv1, gi[0], end, True])
                all_nodes.append(nodes1)
            elif not entry and not math.isnan(path_toa):
                all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, curv1, gi[0], end, False])
                all_nodes.append(nodes1)
            if np.abs(target - (path_toa+new_total))<=1 and not math.isnan(path_toa):
                if entry:
                    x_close = np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5)
                    y_close = np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)
             
                    if x_close and y_close:
                        valid = True
                        break
                
                else:
                    valid = True
                    break
                    
            # if plans[6] and (path_toa+new_total) - plans[5] >= t_it:
            #     buffer = 0
            # elif not plans[6] and (path_toa+new_total) - plans[5] >= t_it:
            #     buffer = 0
            # else: 
            #     target = t_og
            #     buffer = t_it - ((path_toa+new_total)-plans[5]) + 1
            #     target+=buffer
            r_star_time = target-path_toa
            r_star_len = r_star_time*v
            gi[1] = get_guess_point(travel_times=toa_local, dists = dists_local, required_toa=r_star_time, required_len=r_star_len, index = gi[0])
            if (path_toa+new_total) > target:
                p1g+=0.25
            else:
                p1g-=0.5
            if p1g >=1:
                p1g =0.5
            

        gi[0]+=1

    stop_time = time.time()
    print(np.sum(opt))
    diffs = []
    # print(all_nodes)
    for i in all_bez:
        # if plans[6]:
        diffs.append(np.abs(target-i[0]))
        # else:
            # if i[0] - target>=-3:
                # diffs.append(i[0]-target)
        # print(np.abs(target-i[0]), i[len(i)-2])
    # if diffs:  # Check if diffs is not empty
    mindex = diffs.index(np.min(diffs))
    # else:
        # for i in all_bez:
            # diffs.append(np.abs(target-i[0]))
        # mindex = diffs.index(np.min(diffs))
    optim_b = [all_bez[mindex][1], all_bez[mindex][2]]
    nodes1 = all_nodes[mindex]
    print(f'TOTAL OPTIMIZATION TIME: {stop_time-start_time}')
    
    # print("ENTRY TOA:", entryTOA)
    # print("Bez TOA:", pb_toa)
    # print(f'Target Total is {target}\nTarget Bez TOA is: {t_t}\nBez+Entry Path TOA is: {path_toa}\nTotal Travel Time w/ STAR Route is {path_toa+new_total}\nValidity is {valid}')
    # # print(curv1, turn_radius/.3048)
    # print(f'STAR Path Intercept is {gi}')

    
    gi[0] = all_bez[mindex][len(all_bez[mindex])-3]

    end_point = get_line(wps[gi[0]], wps[gi[0]+1], all_bez[mindex][len(all_bez[mindex])-1])
    new_travel, new_total = getTravelXY(wps, end_point, gi[0], toa_local, v)
    # print(new_total)
    
    rotated_bezier = rotate_bez(np.array(optim_b), -angle_to_rotate, origin)
    rot_nodes = rotate_bez(nodes1, -angle_to_rotate, origin)
    # print(rot_nodes)
    entry = all_bez[mindex][-1]
    try:
        ep = get_line(wps[gi[0]], wps[gi[0]+1], gi[1])
        m2 = (wps[gi[0]][1] - ep[1])/(wps[gi[0]][0] - ep[0])
    except:
        ep = get_line(wps[gi[0]], wps[gi[0]+1], gi[1])
        t_h = np.arctan2(wps[gi[0]+1][1] - wps[gi[0]][1],wps[gi[0]+1][0] - wps[gi[0]][0])
        m2 = (wps[gi[0]][1] - (ep[1]+lr*1000*np.abs(np.sin(t_h))))/(wps[gi[0]][0] - (ep[0]+lr*1000*np.abs(np.cos(t_h))))
    b2 = wps[gi[0]][1] - m2*wps[gi[0]][0]
    # plt.plot(rotated_entry[0], rotated_entry[1])
    # plt.plot(rotated_bezier[0], rotated_bezier[1])
    show_local = False
    print(f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}\nDifference: {np.abs(target-all_bez[mindex][0]):0.2f}')
    
    plt.plot(x, y, linewidth = 2, label = 'STAR Path')
    
    if all_bez[mindex][-1]:
        req_ent = all_bez[mindex][11][1]
        act_ent = all_bez[mindex][11][0]

        rotated_partial = rotate_bez(np.array([all_bez[mindex][1], all_bez[mindex][2]]), -angle_to_rotate, origin)
        rotated_entry = rotate_bez(np.array([all_bez[mindex][4], all_bez[mindex][5]]), -angle_to_rotate, origin)
        x_entry = all_bez[mindex][4]
        y_entry = all_bez[mindex][5]
        print(f'Entry Bank{np.rad2deg(all_bez[mindex][9][0])}, Entry Point On Bez {all_bez[mindex][9][1]}')
        # if req_ent-angle_to_rotate>360:
        #     req_ent-=360
        #     act_ent-=360
        #     if req_ent>+360:
        #         req_ent-=360
        #         act_ent-=360
        #     print('REQUIRED HEADING AT ENTRY:', req_ent-angle_to_rotate)
        #     print('ACTUAL HEADING AT INTERSECT:', act_ent-angle_to_rotate)
        # elif req_ent-angle_to_rotate<0:
        #     req_ent+=360
        #     act_ent+=360
        #     print('REQUIRED HEADING AT ENTRY:', req_ent-angle_to_rotate)
        #     print('ACTUAL HEADING AT INTERSECT:', act_ent-angle_to_rotate)
        # else:
        print('REQUIRED HEADING AT ENTRY:', req_ent-angle_to_rotate)
        print('ACTUAL HEADING AT INTERSECT:', act_ent-angle_to_rotate)
        print('DIFF AT ENTRY:', req_ent-act_ent)
        print('BEZ CURVATURE:', all_bez[mindex][12], 'TURN RADIUS:', turn_radius)
        if show_local:
            '''
            Not Rotated -> Local frame
            '''
            # plt.plot(x_entry, y_entry, label = 'Entry Arc', color = 'cyan')
            # plt.scatter(x_entry[-1], y_entry[-1], marker = '*', color = 'yellow', label = 'Intersection Point', s = 100, zorder = 100)
            # plt.text(x_entry[-1]+1000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
            # plt.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial Bez Path')
            # plt.text(partial_bez1[0][50], partial_bez1[1][50]-2500, f'Partial Bez Travel Time: {pb_toa:0.2f}s', weight = 'bold' )


        '''
        'Rotated' -> planned trajectory in global frame
        '''
        plt.plot(rotated_entry[0], rotated_entry[1], color = 'cyan', linewidth = 2, label = 'Entry Arc')
        plt.scatter(rotated_entry[0][-1], rotated_entry[1][-1], marker = '*', color = 'yellow', s = 100, zorder = 100, label = 'Intersection Point')
        plt.text(rotated_entry[0][0]+5000,rotated_entry[1][0]-2500, f'Entry Path Travel Time: \n{all_bez[mindex][10]:0.2f}', weight = 'bold')
        plt.plot(rotated_partial[0], rotated_partial[1], c = 'magenta', linewidth = 2, label = 'Partial Bez Traj')
        plt.text(rotated_bezier[0][10], rotated_bezier[1][10]-5000, f'Partial Bez Travel Time: \n{all_bez[mindex][6]:0.2f}s', weight = 'bold' )
        plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--', linewidth = 2)
    else:
        plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--', linewidth = 2, label = 'Bez Traj')
        print(f'No Entry Path Used\nRequired Bez Heading: {np.rad2deg(np.arctan2(rotated_bezier[1][0] - rotated_bezier[1][1], rotated_bezier[0][0] - rotated_bezier[0][1]))}')
        print(f'Actual Heading: {np.rad2deg(head)}')
        print('BEZ CURVATURE:', all_bez[mindex][6], 'TURN RADIUS:', turn_radius)
        entry = False
        # plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        plt.text(rotated_bezier[0][25], rotated_bezier[1][25]-2500, f'Bez Travel Time: {all_bez[mindex][4]:0.2f}s', weight = 'bold' )
    # if head > pi == 'IV':
    #     plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head+np.pi), 2500*np.sin(head+np.pi), width=50,
    #                 length_includes_head=False,
    #                 head_width=1000,
    #                 head_length=1000,
    #                 head_starts_at_zero=False,
    #                 facecolor='black',
    #                 edgecolor='black',
    #                 label = f'Heading of {np.rad2deg(head)}')
    # else:
    plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head), 2500*np.sin(head), width=50,
                length_includes_head=False,
                head_width=1000,
                head_length=1000,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                label = f'Heading of {np.rad2deg(head):0.2f}')
    # plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(np.deg2rad(h)), 2500*np.sin(np.deg2rad(h)), width=50,
    #             length_includes_head=False,
    #             head_width=1000,
    #             head_length=1000,
    #             head_starts_at_zero=False,
    #             facecolor='purple',
    #             edgecolor='purple')#,
    #             #label = f'Heading of {np.rad2deg(head):0.2f}')
    if show_local:
        '''
        Not Rotated -> Local frame 
        '''
        # for i in range(len(wps)-1):
        #     plt.scatter(interps_local[i][0], interps_local[i][1], color = 'orange')

        # ang_y = (nodes[0][2] - 5000) *m + b
        # plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
        # plt.plot(bez[0], bez[1], color = 'black', linestyle = '--')
        # plt.plot(optim_b[0], optim_b[1], color = 'black', linestyle = '--')
        # plt.scatter(optim1[0], optim1[1], marker = '^', color = 'red')
        # plt.scatter(nodes[0][2], nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
        # plt.scatter(path_startX_local, path_startY_local, marker = '^', color = 'red')
        # plt.scatter(x_rot, y_rot, label = 'STAR Path')
        # plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

    '''
    'Rotated' -> planned trajectory in global frame
    '''
    for i in range(len(wps)-1):
        plt.scatter(wps[i][0], wps[i][1], color = 'orange', marker = 's')
    ang_y = (rotated_bezier[0][-1] + 5000) *m2 + b2
    plt.plot([rotated_bezier[0][-1],rotated_bezier[0][-1] + 5000], [rotated_bezier[1][-1], ang_y], color = 'green', linestyle = '--')
    plt.scatter(rot_nodes[0][1], rot_nodes[1][1], marker = '^', color = 'red', zorder = 10)
    # plt.scatter(rot_nodes[0][2], rot_nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
    plt.scatter(rotated_bezier[0][-1], rotated_bezier[1][-1], marker = '^', color = 'red')#, label = 'Control Points')

    plt.scatter(path_startX, path_startY, marker = '^', color = 'red')
    
    # plt.text(wps[0][0], wps[0][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    # plt.text(wps[1][0], wps[1][1], f'STAR Travel Time: {all_bez[mindex][3]:0.2f}s', weight = 'bold' )
    # plt.text(wps[-1][0], wps[-1][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    plt.text(0, 0, f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    plt.text(wps[2][0]+1000, wps[2][1], f'STAR Travel Time: \n{all_bez[mindex][3]:0.2f}s', weight = 'bold' )
    # nodes_r = rotate_bez(all_nodes[4], -angle_to_rotate, origin)
    # new_b = manual_bez([nodes_r[0][0], nodes_r[1][0]], [nodes_r[0][1], nodes_r[1][1]], [nodes_r[0][2], nodes_r[1][2]], 200)
    # plt.plot(new_b[0], new_b[1], linewidth = 3, color = 'red')
    plt.scatter(ac_met[0], ac_met[1], color = 'green')
    plt.scatter(0, 0, marker = '*', color = 'green')#, label = 'CMH')
    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    plt.grid()
    plt.title(f'Quadrant {gate[0]}, Path Optimization Time: {(stop_time-start_time):0.2f}')
    plt.legend(loc = 'upper left')
    plt.axis('equal')
    # plt.plot([0, ac_met[0]], [0, ac_met[1]], linestyle = '--')
    # plt.text(ac_met[0]/2, ac_met[1]/2, f'Angle of {(np.rad2deg(np.arctan2(ac_met[1], ac_met[0]))+180):0.1f} Heading of {heading:0.1f}')
    if plot:
        plt.show()
    else:
        plt.close()
    
    if not all_bez[mindex][-1]:
        req_ent = np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])-90
        rotated_entry = []
        ba = 0
    else:
        req_ent= req_ent-angle_to_rotate
    if gate == 'I':

        if  0<=gi[0]<= 1:
            fl = 8000
        elif gi[0] == 2:
            fl = 7000
        elif gi[0] > 2:
            fl = 6000
        else:
            fl = 10000
    if gate == 'II':
        if gi[0] == 1:
            fl = 8000
        elif gi[0] >=3:
            fl = 6000
        else:
            fl = 10000
    if gate == 'III':
        if gi[0] == 1:
            fl = 7000
        elif gi[0] == 2:
            fl = 7000
        elif gi[0] >=3:
            fl = 6000
        else:
            fl = 10000
    if gate == 'IV':
        fl = 5000
    # if gate == 'I':
    #     if ac_met[1] < 
    # if req_ent>=360:
    #     req_ent-=360
    # elif req_ent < 0:
    #     req_ent+=360
    # if gate == 'II':
    #     req_ent -= np.pi
    # elif gate == 'III':
    #     req_ent=req_ent
    # elif gate == 'IV':
    #     req_ent += np.pi
    # elif gate == 'I':
    #     req_ent=3*np.pi/2+req_ent

    # if not all_bez[mindex][-1]:
    #     req_ent = np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])
    #     rotated_entry = []
    #     ba = 0
    if entry and gate !='III':#and lr == 1 
        return wps, end_point, rot_nodes, rotated_bezier, [entry, rotated_entry, 90-req_ent, np.rad2deg(all_bez[mindex][9][0][0]), [rotated_entry[0][-1], rotated_entry[1][-1]]], fl, stop_time-start_time, [target-all_bez[mindex][0], all_bez[mindex][0]], lr
    # elif entry and lr == -1 and gate!='III':
    #     return wps, end_point, rot_nodes, rotated_bezier, [entry, rotated_entry, 90-req_ent, np.rad2deg(all_bez[mindex][9][0][0]), [rotated_entry[0][-1], rotated_entry[1][-1]]], stop_time-start_time, target-all_bez[mindex][0], lr       
    elif entry and gate == 'III':
        print(req_ent)
        req_ent =90-req_ent+360
        print(req_ent)
        return wps, end_point, rot_nodes, rotated_bezier, [entry, rotated_entry, req_ent, np.rad2deg(all_bez[mindex][9][0][0]), [rotated_entry[0][-1], rotated_entry[1][-1]]], fl,  stop_time-start_time, [target-all_bez[mindex][0], all_bez[mindex][0]], lr
    else:
        return wps, end_point, rot_nodes, rotated_bezier, [entry, rotated_entry, req_ent, np.rad2deg(ba)], fl, stop_time-start_time, [target-all_bez[mindex][0], all_bez[mindex][0]], lr


if __name__ == '__main__':

    v = 57.412
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    kcmh = [39.997472, -82.891194]
    teeze = [40.0865, -82.848111]
    # TEEZE WAYPOINTS
    melzz = [40.313417, -83.247194]
    dubln = [40.202944, -83.132306]
    trlgy = [40.167944, -83.061139]
    polrs = [40.111194, -82.953444]
    taces = [40.090111, -82.916333]
    

    teeze_path = [melzz, dubln, trlgy, polrs, taces, teeze]
    
    jaktz = [39.591028, -83.419583]
    rscot = [39.722389,  -83.286306]
    obetz = [39.787667, -83.158389]
    edwib = [39.877472, -82.984861]
    gagbe = [39.907167, -82.927278]
    jesce = [39.903556, -82.858889]
    elupy = [39.831083, -83.074778]


    elupy_path = [rscot, obetz, elupy, edwib, gagbe, jesce]


    # XAVYR WAYPOINTS
    scrlt = [39.502917, -82.350833]
    brtus = [39.730944, -82.473083]
    guber = [39.963222, -82.670889]
    bkeye = [39.982056, -82.706694]
    xavyr = [39.906167, -82.647333]
    xavyr_path = [brtus, xavyr, guber, bkeye]

    # FAVUS WAYPOINTS
    bugzz = [40.565, -82.454056]
    cbuss = [40.325306, -82.5405]
    molls = [40.132139, -82.609194]
    ordiw = [39.989306, -82.666056]
    bazel = [39.981917, -82.704556]
    favus = [40.058917, -82.639417]


    favus_path =[cbuss, molls, favus, ordiw, bazel]
    
    ac_pos = [40.10373464575462, -82.30768892311026] #g1
    h=257.97
    target = 731.8725791205976 #G1

    # ac_pos = [40.32447575882298 , -82.48594935697894]
    # h = 225
    # target = 794.6332769682364
    
    ac_pos = [40.38110069057144, -82.51492919797532]
    h = 52+180
    target = 880.3137447967732

    # ac_pos = [40.44288645743566, -82.7322246512964]
    # h = 196.86482979260154
    # target = 988.2709814198248


    # ac_pos = [40.42655576925774, -82.6248931230084]
    # h=206.9461709913814
    # target = 954.439730218743

    heading=270-h
    head = np.deg2rad(heading-180)
    angle_to_rotate = (90-heading+180)
    gate = ['I']
    wps = LongLat_To_WSG84_Meters(favus_path, kcmh) #Gives waypoints in global frame
    wps_local = LongLat_To_WSG84_Meters(favus_path, ac_pos)

    # ac_pos = [40.303627035156936 , -83.37719169326613] #g2
    # h =131
    # target = 902.1233635800432 #G2

    # ac_pos = [40.45110099161172, -83.07614134313081]
    # h = 164
    # target = 965.4033554973976

    # ac_pos = [40.16910942656836, -83.46660171897712]
    # h = 113
    # target = 978.4033554973976

    # heading = h+180
    # head = np.deg2rad(heading)
    # angle_to_rotate = (90-heading)
    # gate = ['II']
    # wps = LongLat_To_WSG84_Meters(teeze_path, kcmh) #Gives waypoints in global frame
    # wps_local = LongLat_To_WSG84_Meters(teeze_path, ac_pos)

    # ac_pos = [39.58268448211492, -83.14955754333681] #g3
    # h=27
    # target = 857.1282220067823 #G3
    # heading = h

    # ac_pos = [39.59273116598905, -83.15156205183662]
    # h=28.201610030646435
    # target=838.973276722232


    # heading=h
    # head = np.deg2rad(heading)
    # angle_to_rotate = (90-heading)
    # gate = ['III']
    # wps = LongLat_To_WSG84_Meters(elupy_path, kcmh) #Gives waypoints in global frame
    # wps_local = LongLat_To_WSG84_Meters(elupy_path, ac_pos)


    # ac_pos = [39.78201288955758, -82.35136330299837] #g4
    # h=299
    # target = 679.2001084334775 #G4

    # ac_pos = [39.90215225350915, -82.3169331457393]
    # h = 283.4527388927732
    # target = 676.850662875645

    # ac_pos = [39.9711138746171, -82.30444769380334]
    # h=274.6297714128176
    # target = 694.0649856906898

    # ac_pos = [39.62515627391903, -82.50350114681574]
    # h = 322.58407819436576 
    # target=754.2018646466217


    # heading = h-180
    # head = np.deg2rad(heading)
    # angle_to_rotate = (90-heading)
    # gate = ['IV']
    # wps = LongLat_To_WSG84_Meters(xavyr_path, kcmh) #Gives waypoints in global frame
    # wps_local = LongLat_To_WSG84_Meters(xavyr_path, ac_pos)


    rotate_head = np.deg2rad(heading+angle_to_rotate)

    print(f'HEADING: {h}, ROTATE ANGLE: {angle_to_rotate}, ROTATED HEADING: {h+angle_to_rotate}')
    
    # head = np.deg2rad(h)
    # h=270-h
    # angle_to_rotate = (90-h)
    # print(transform_angle)

    
    ac_met = __WSG84_To_Meters_Single(ac_pos, kcmh) #Global aircraft position

    

    bank = 5
    min = 5
    plan_ahead = 10

    

    x = [i[0] for i in wps]
    y = [i[1] for i in wps]

 
    origin = np.array([ac_met])

    rotated_wps = rotate_points(wps, angle_to_rotate, origin)

    x_rot = [i[0] for i in rotated_wps]
    y_rot = [i[1] for i in rotated_wps]

    d_reg = np.hypot(x[0]-ac_met[0], y[0]-ac_met[1])
    d_rot = np.hypot(x_rot[0]-ac_met[0], y_rot[0]-ac_met[1])

    dists_local, toa_local, total_local = total_travel(wps_local, v)

    # for i in range(len(wps)-1):
    #     plt.plot([rotated_wps[i][0], rotated_wps[i+1][0]], [rotated_wps[i][1], rotated_wps[i+1][1]], color = 'orange')
    #     plt.plot([wps[i][0], wps[i+1][0]], [wps[i][1], wps[i+1][1]])
    # # plt.plot(x_local, y_local)
    # plt.scatter(origin[0][0], origin[0][1])
    # plt.show()
    if rotated_wps[1][0] < ac_met[0]:
        lr = -1
    else:
        lr = 1
    
    
    path_startX = ac_met[0] + np.cos(head)*v*plan_ahead 
    path_startY = ac_met[1] + np.sin(head)*v*plan_ahead 

    path_startX_local = ac_met[0]
    path_startY_local = ac_met[1] + v*plan_ahead

    
    
   
    turn_radius = v*1.94384**2/(11.26*math.tan(np.deg2rad(73)))
    turn_radius*=0.3048
    
  
    gi = [0,0]
    t_t = target-total_local-plan_ahead
   
    end_point = get_line(rotated_wps[1], rotated_wps[2], 0)
    end_point_ref = get_line(rotated_wps[1], rotated_wps[2], 0.1)
    new_travel, new_total = getTravelXY(wps, end_point, gi[0], toa_local, v)
    t_h = np.arctan2(end_point_ref[1]- end_point[1],end_point_ref[0]- end_point[0])
    t_t = target-new_total-plan_ahead
    end_point = get_line(rotated_wps[1], rotated_wps[2], gi[1])
    entry = False
    entryTOA = plan_ahead
    all_bez = []
    all_nodes = []
    opt = []
    show = False
    
    # for i in rotated_wps:
    #     if ac_met[1]-i[1] >0:
    #         gi[0]+=1
            
    start_time = time.time()

    for i in range(gi[0], len(wps)-1):
        # g = 5000
        gi[1] = 0.6
        p1g = 0.25
        valid = False
        w1 = rotated_wps[gi[0]]
        w2 = rotated_wps[gi[0]+1]
        m = (w2[1] - w1[1])/(w2[0] - w1[0])
        b = w1[1] - m*w1[0]
        coeffs = [m, b]
        t_h = np.arctan2(w2[1]- w1[1],w2[1]- w1[0])
        ang_y = (w2[0] + 5000) *m + b
        
        # while gi[1] < 1:
        for j in range(0,11):
            print(gi, p1g)

            end_point = get_line(w1, w2, gi[1])
            p1_point = get_line_p1(w1, end_point, p1g, t_h)
            new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            print(f'NEW STAR TOTAL TIME: {new_total}')
            t_t = target-new_total-entryTOA
            print(f'BEZ TARGET TIME: {t_t}')
            # nodes = [np.array([path_startX_local, end_point[0]+lr*g*np.abs(np.cos(t_h)), end_point[0]]).flatten(),
            #     np.array([path_startY_local, end_point[1]+lr*g*np.abs(np.sin(t_h)), end_point[1]]).flatten()]
            nodes = [np.array([path_startX_local, p1_point[0], end_point[0]]).flatten(),
                np.array([path_startY_local, p1_point[1], end_point[1]]).flatten()]
            # nodes = [np.array([path_startX_local, w1[0], end_point[0]]).flatten(),
            #     np.array([path_startY_local, w1[1], end_point[1]]).flatten()]
            
           
            # end_point_ref = get_line(w1, w2, 1)
            
            print(f'TARGET HEADING: {np.rad2deg(t_h)}')
            opt_time = time.time()

            # optim1, end, curv1 = solve_optim12(P0 = [nodes[0][0],nodes[1][0]],
            #                             target_toa=t_t,
            #                             guess = [nodes[0][1], nodes[1][1], gi[1]],
            #                             target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs, wps=[w1, w2])
            p1, end, curv1 = solve_optim12(P0 = [nodes[0][0],nodes[1][0]],
                                        target_toa=t_t,
                                        guess = [p1g, gi[1]],
                                        target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs, wps=[w1, w2])
            print(f'GUESS: {gi[1]}, RESULT: {end}')
            p1g = p1
            # print('END POINT', (gi[0], end))
            # optim1, end, curv1 = solve_optim12_o(
            #     nodes[0], t_t, nodes[1], t_h, v, turn_radius, lr, coeffs, [w1,w2]
            # )
            opt.append(time.time()-opt_time)
            # print(f'OPTIMIZATION TIME: {time.time()-opt_time}')
            end_point = get_line(w1, w2, end)
            # end_point = get_line_o(w1,w2,end)
            optim1 = get_line_p1(w1, end_point, p1, t_h)
            new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            # # new_travel, new_total = get_travel_xy(rotated_wps, np.array(end_point), gi[0], toa_local, v)
            t_t = target-new_total-entryTOA

            if optim1[0] < ac_met[0]:
                lr = -1
            else:
                lr = 1

            optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [end_point[0], end_point[1]], 50)
            # optim_b = manual_bez_o(np.array([nodes[0][0],nodes[1][0]]), np.array([optim1[0], optim1[1]]), np.array([end_point[0], end_point[1]]), 50)
            optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],1)
            # optim_b_length = path_length_o(np.array([optim1[0], optim1[1]]),np.array([nodes[0][0],nodes[1][0]]),np.array([end_point[0], end_point[1]]),1)
            path_toa = optim_b_length/v
            print(f'REQUESTED: {t_t}, ACTUAl: {path_toa}')
            nodes1 = [np.array([nodes[0][0], optim1[0], end_point[0]]).flatten(), np.array([nodes[1][0], optim1[1], end_point[1]])]
            if math.isnan(path_toa):
                    optim1[0]-=1000*np.cos(t_h)
                    optim1[1]-=1000*np.sin(t_h)
                    optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],1)
                    path_toa = optim_b_length/v

            # print('HEAIDNG:',np.rad2deg(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])), np.rad2deg(np.pi/2), np.deg2rad(10))
            if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(15)):
                entryTOA = plan_ahead
                path_toa = optim_b_length/v + entryTOA
                entry = False
                # all_bez.append([optim_b[0], optim_b[1], optim_b_length, path_toa])
                # print('using plan ahead TOA')
            else:
                # show = True
                entry = True
                bez_ang = np.rad2deg(np.arctan2(optim_b[1][1]-optim_b[1][0], optim_b[0][1]-optim_b[0][0]))
                # if i ==1:
                #     show = True
                # else:
                #     show = False
                # if bez_ang<0:
                #     bez_ang+=360
                # if bez_ang>=270:
                #     bank = 42.5
                # else:
                #     bank = 5
                # print(f'BEZ ANGLE: {bez_ang}, BANK {bank}')
                opt_time = time.time()
                ba, t_entry = solve_optimEntry([np.deg2rad(bank), 0.005], np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
                # ba[0] = np.deg2rad(35)
                # t_entry+=0.005
                print('T ENTRY', t_entry)
                # ba, t_entry = solve_optimEntry_o(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
                opt.append(time.time()-opt_time)
                # print(f'OPTIMIZATION TIME: {time.time()-opt_time}')

                x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [end_point[0], end_point[1]], t_entry)
                # x_int_en, y_int_en = find_bez_xy_o(np.array([nodes1[0][0],nodes1[1][0]]),
                #                             np.array([optim1[0],optim1[1]]),
                #                             np.array([end_point[0], end_point[1]]), t_entry)
                x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [end_point[0], end_point[1]], t_entry+0.001)
                # x_int_en2, y_int_en2 = find_bez_xy_o(np.array([nodes1[0][0],nodes1[1][0]]),
                #                             np.array([optim1[0],optim1[1]]),
                #                             np.array([end_point[0], end_point[1]]), t_entry+0.01)
                req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
                
                
                x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotate_head)
                # x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath_o(v, np.rad2deg(ba), np.array([x_int_en, y_int_en]), np.array(ac_met), lr, rotate_head)
                act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
                

                partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [end_point[0], end_point[1]], 50, t_entry)
                # partial_bez1 = manual_bez_partial_o(np.array([nodes1[0][0],nodes1[1][0]]),
                #                             np.array([optim1[0],optim1[1]]),
                #                             np.array([end_point[0], end_point[1]]), 50, t_entry)
                pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],t_entry)
                # pb_length = optim_b_length - path_length_o(np.array([optim1[0], optim1[1]]),np.array([nodes[0][0],nodes[1][0]]),np.array([end_point[0], end_point[1]]),t_entry)
                pb_toa = pb_length/v
                # print(f'PB TOA: {pb_toa}')
                path_toa = entryTOA + pb_toa
                print(f'TOA: {path_toa+new_total}, NEW TOTAL {new_total}, TARGET: {target}, PATH TOA: {path_toa}, DIFF: {np.abs(target - (path_toa+new_total))}')
            ang_y = (end_point[0] - 2500) *m + b
            plt.plot([end_point[0],end_point[0] - 5000], [end_point[1], ang_y], color = 'green', linestyle = '--')
            plt.plot(x_rot, y_rot)
            plt.scatter(ac_met[0], ac_met[1])
            plt.plot(x_rot, y_rot)
            plt.scatter(optim1[0], optim1[1], marker = '^', color = 'red')
            plt.plot(optim_b[0], optim_b[1], color = 'black')
            plt.scatter(end_point[0], end_point[1], marker = 's', color = 'black')
            plt.axis('equal')
            plt.text(rotated_wps[-1][0]-10000, rotated_wps[-1][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
            plt.text(rotated_wps[1][0], rotated_wps[1][1], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )
            plt.title(f'{gi[0]:0.2f}, {gi[1]:0.2f}, {end:0.2f}, P1: {p1}')
                # if math.isnan(path_toa):
                #     print(optim_b[0])
                #     print(optim_b[1])
                #     print(nodes1)
                #     print(pb_length, optim_b_length)
                #     print(x_entry)
                #     print(y_entry)
            if entry:
                plt.plot(x_entry, y_entry, color = 'cyan')
                plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
                plt.text(optim_b[0][25], optim_b[1][25]-5000, f'Partial Bez Travel Time: {pb_toa:0.2f}s', weight = 'bold' )
                print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA:0.2f}')
                plt.text(x_entry[0], y_entry[0], f'BANK ANGLE: {np.rad2deg(ba[0]):0.4f}', weight = 'bold')
            else:
                print(f'Bez Travel Time: {path_toa:0.2f}s')
                plt.text(optim_b[0][25], optim_b[1][25], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
                
            plt.show()
                
            # if 1<=np.abs(target - (path_toa+new_total))<=5 and not math.isnan(path_toa):
            #     og = gi[1]
            #     for i in range(0, 5):
            #         w1 = rotated_wps[gi[0]]
            #         w2 = rotated_wps[gi[0]+1]
            #         m = (w2[1] - w1[1])/(w2[0] - w1[0])
            #         b = w1[1] - m*w1[0]
            #         coeffs = [m, b]
            #         # gi[1]+=0.05
            #         if w2[0] <= w1[0]:
            #             if path_toa+new_total > target:
            #                 g-=2500
            #                 gi[1]+=0.05
            #             else:
            #                 g+=2500
            #                 gi[1]-=0.05
            #         else:
            #             if path_toa+new_total > target:
            #                 g+=2500
            #                 gi[1]+=0.05
            #             else:
            #                 g-=2500
            #                 gi[1]-=0.05
            #         # print(g, gi)
            #         # print('HERHE HERE HERE')

            #         end_point = get_line(rotated_wps[gi[0]], rotated_wps[gi[0]+1], gi[1])
            #         # end_point = get_line_o(rotated_wps[gi[0]], rotated_wps[gi[0]+1], gi[1])
            #         show = False
            #         new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            #         # new_travel, new_total = get_travel_xy(rotated_wps, np.array(end_point), gi[0], toa_local, v)
            #         t_t = target-new_total-entryTOA
                    
            #         nodes = [np.array([path_startX_local, end_point[0]+lr*g*np.abs(np.cos(t_h)), end_point[0]]).flatten(),
            #             np.array([path_startY_local, end_point[1]+lr*g*np.abs(np.sin(t_h)), end_point[1]]).flatten()]
                
            #         # m = (rotated_wps[gi[0]+1][1] - rotated_wps[gi[0]][1])/(rotated_wps[gi[0]+1][0] - rotated_wps[gi[0]][0])
            #         # b = rotated_wps[gi[0]][1] - m*rotated_wps[gi[0]][0]
            #         # coeffs = [m, b]

            #         # end_point_ref = get_line(rotated_wps[gi[0]], rotated_wps[gi[0]+1], 1)
            #         # end_point_ref = get_line_o(rotated_wps[gi[0]], rotated_wps[gi[0]+1], 1)

            #         # t_h = np.arctan2(end_point_ref[1]- end_point[1],end_point_ref[0]- end_point[0])
            #         t_h = np.arctan2(w2[1]- end_point[1],w2[1]- end_point[0])
            #         opt_time = time.time()
            #         optim1, end, curv1 = solve_optim12(P0 = [nodes[0][0],nodes[1][0]],
            #                                     target_toa=t_t,
            #                                     guess = [nodes[0][1], nodes[1][1], gi[1]],
            #                                     target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs, wps=[w1, w2])
            #         # optim1, end, curv1 = solve_optim12_o(P0 = [nodes[0][0],nodes[1][0]],
            #         #                             target_toa=t_t,
            #         #                             guess = [nodes[0][1], nodes[1][1], gi[1]],
            #         #                             target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs, wps=[w1, w2])
            #         opt.append(time.time()-opt_time)
            #         # print(f'OPTIMIZATION TIME: {time.time()-opt_time}')
            #         end_point = get_line(rotated_wps[gi[0]], rotated_wps[gi[0]+1], end)
            #         # end_point = get_line_o(rotated_wps[gi[0]], rotated_wps[gi[0]+1], end)
            #         new_travel, new_total = getTravelXY(rotated_wps, end_point, gi[0], toa_local, v)
            #         # new_travel, new_total = get_travel_xy(rotated_wps, np.array(end_point), gi[0], toa_local, v)


            #         if optim1[0] < ac_met[0]:
            #             lr = -1
            #         else:
            #             lr = 1
            #         optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [end_point[0], end_point[1]], 50)
            #         # optim_b = manual_bez_o(np.array([nodes[0][0],nodes[1][0]]), np.array([optim1[0], optim1[1]]), np.array([end_point[0], end_point[1]]), 50)

            #         optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],1)
            #         # optim_b_length = path_length_o(np.array([optim1[0], optim1[1]]),np.array([nodes[0][0],nodes[1][0]]),np.array([end_point[0], end_point[1]]),1)

            #         path_toa = optim_b_length/v

            #         nodes1 = [np.array([nodes[0][0], optim1[0], end_point[0]]).flatten(), np.array([nodes[1][0], optim1[1], end_point[1]])]

            #         if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(10)):
            #             entryTOA = plan_ahead
            #             path_toa = optim_b_length/v + entryTOA
            #             entry = False
            #             # all_bez.append([optim_b[0], optim_b[1], optim_b_length, path_toa])
            #             # print('using plan ahead TOA')
            #         else:
            #             entry = True
            #             bez_ang = np.rad2deg(np.arctan2(optim_b[1][1]-optim_b[1][0], optim_b[0][1]-optim_b[0][0]))
            #             # if bez_ang<0:
            #             #     bez_ang+=360
            #             # if bez_ang>=270:
            #             #     bank = 42.5
            #             # else:
            #             #     bank = 5
            #             # print(f'BEZ ANGLE: {bez_ang}, BANK {bank}')
            #             opt_time = time.time()
            #             ba, t_entry = solve_optimEntry([np.deg2rad(bank), 0.005], np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
            #             # ba[0] = np.deg2rad(35)
            #             # t_entry+=0.005
            #             # print('T ENTRY', t_entry)
            #             # t_entry+=0.01
            #             # ba, t_entry = solve_optimEntry_o(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
            #             opt.append(time.time()-opt_time)
                        

            #             x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
            #                                         [optim1[0],optim1[1]],
            #                                         [end_point[0], end_point[1]], t_entry)
                        
            #             x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
            #                                         [optim1[0],optim1[1]],
            #                                         [end_point[0], end_point[1]], t_entry+0.001)
                        
            #             # x_int_en, y_int_en = find_bez_xy_o(np.array([nodes1[0][0],nodes1[1][0]]),
            #             #                     np.array([optim1[0],optim1[1]]),
            #             #                     np.array([end_point[0], end_point[1]]), t_entry)
                        
            #             # x_int_en2, y_int_en2 = find_bez_xy_o(np.array([nodes1[0][0],nodes1[1][0]]),
            #             #                     np.array([optim1[0],optim1[1]]),
            #             #                     np.array([end_point[0], end_point[1]]), t_entry+0.01)
            #             req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
                        
                        
            #             x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotate_head)
            #             # x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath_o(v, np.rad2deg(ba), np.array([x_int_en, y_int_en]), np.array(ac_met), lr, rotate_head)
            #             act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
                        

            #             partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
            #                                     [optim1[0],optim1[1]],
            #                                     [end_point[0], end_point[1]], 50, t_entry)
            #             # partial_bez1 = manual_bez_partial_o(np.array([nodes1[0][0],nodes1[1][0]]),
            #             #                     np.array([optim1[0],optim1[1]]),
            #             #                     np.array([end_point[0], end_point[1]]), 50, t_entry)
            #             pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[end_point[0], end_point[1]],t_entry)
            #             # pb_length = optim_b_length - path_length_o(np.array([optim1[0], optim1[1]]),np.array([nodes[0][0],nodes[1][0]]),np.array([end_point[0], end_point[1]]),t_entry)

            #             pb_toa = pb_length/v
            #             path_toa = entryTOA + pb_toa
                        

            #         ang_y = (end_point[0] - 5000) *m + b
            #         # plt.plot([end_point[0],end_point[0] - 5000], [end_point[1], ang_y], color = 'green', linestyle = '--')
            #         # plt.plot(x_rot, y_rot)
            #         # plt.scatter(ac_met[0], ac_met[1])
            #         # plt.plot(x_rot, y_rot)
            #         # plt.scatter(optim1[0], optim1[1], marker = '^', color = 'red')
            #         # plt.plot(optim_b[0], optim_b[1], color = 'black')
            #         # plt.scatter(end_point[0], end_point[1], marker = 's', color = 'black')
            #         # plt.axis('equal')
            #         # plt.text(rotated_wps[-1][0], rotated_wps[-1][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
            #         # plt.text(rotated_wps[1][0], rotated_wps[1][1], f'STAR Travel Time: {new_total:0.2f}s\nIn Region', weight = 'bold' )

            #         # if entry:
            #         #     plt.plot(x_entry, y_entry, color = 'cyan')
            #         #     plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
            #         #     plt.text(optim_b[0][25], optim_b[1][25]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
            #         #     plt.text(x_entry[0], y_entry[0], f'BANK ANGLE: {np.rad2deg(ba[0]):0.2f}', weight = 'bold')
            #         # #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
            #         # else:
            #         # #     print(f'Bez Travel Time: {path_toa:0.2f}s')
            #         #     plt.text(optim_b[0][25], optim_b[1][25], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
            #         # plt.show()
            #         if entry:
            #             x_close = np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5)
            #             y_close = np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)
            #             # print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}, {x_close}')
            #             # print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}, {y_close}')

            #         if entry and not math.isnan(path_toa) and x_close and y_close:
            #             all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, ba, entryTOA, act_ent, gi[0], end, entry])
            #             all_nodes.append(nodes1)
            #         elif not entry and not math.isnan(path_toa):
            #             all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, gi[0], end, entry])
            #             all_nodes.append(nodes1)
            #         if np.abs(target - (path_toa+new_total))<=1 and not math.isnan(path_toa):
                        
            #             if entry == True and  (not x_close or not y_close):
            #                 valid = False
            #             elif entry == True and x_close and y_close:
                       
            #                 valid = True
            #                 break
            #             else:
            #                 valid = True
            #                 break
            #     gi[1] = og
                    
                    # print(new_total + path_toa)
                    # print(f'TARGET: {target}, ACTUAL TRAVEL TIME: {path_toa+new_total}')
                    
                    # print(f'TIME DIFFERENCE: {np.abs(target - (path_toa+new_total))}')
                    # print(f'Original Target: {target}, New STAR Travel Time: {new_total}, New Target Time: {t_t}, Path TOA: {path_toa}')
                
            if entry and not math.isnan(path_toa):
                all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, [ba, t_entry], entryTOA, [act_ent, req_ent], curv1, gi[0], end, entry])
                all_nodes.append(nodes1)
            elif not entry and not math.isnan(path_toa):
                all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, curv1, gi[0], end, entry])
                all_nodes.append(nodes1)
                # if np.abs(target - (path_toa+new_total))<=1 and not math.isnan(path_toa):
                # # print('PATH IS VALID, APPENDINH')
                #     # print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
                #     # print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
                #     x_close = np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5)
                #     y_close = np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)
                #     if entry == True and (not x_close or not y_close):
                #         valid = False
                #         if act_ent > req_ent:
                #             bank-=2.5
                #         else:
                #             bank+=10
                #     elif entry == True and x_close and y_close:
                #         # print(path_toa+new_total, pb_toa, new_total)
                #         all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, ba, entryTOA, act_ent, gi[0], end, entry])
                #         # all_nodes.append(nodes1)
                #         valid = True
                #         break
                #     else:
                #         valid = True
                #         all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, gi[0], end, entry])
                #         all_nodes.append(nodes1)
                #         break
            

            if np.abs(target - (path_toa+new_total))<=1 and not math.isnan(path_toa):
                
                if entry:
                    x_close = np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5)
                    y_close = np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)
                # print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}, {x_close}')
                # print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}, {y_close}')
                    if (not x_close or not y_close):
                        valid = False

                    elif x_close and y_close:

                        valid = True
                        break
                else:
                    valid = True
                    break

            

            # if valid == False:
                # gi[1]+=0.05
                # if entry:
            r_star_time = target-path_toa
            print(f'STAR REQUIRED TIME: {r_star_time}')
            r_star_len = r_star_time*v
            gi[1] = get_guess_point(travel_times=toa_local, dists = dists_local, required_toa=r_star_time, required_len=r_star_len, index = gi[0])
            if (path_toa+new_total) > target:
                p1g+=0.25
            else:
                p1g-=0.5
            if p1g >=1:
                p1g =0.5
            # else:
            #     r_star_time = target-
            # if gi[1] > 1:
            #     gi[1] = 1
                # break

            # if np.abs(target - (path_toa+new_total))<=1 and not math.isnan(path_toa):
            #     # print('PATH IS VALID, APPENDINH')
            #     # print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
            #     # print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
            #     x_close = np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5)
            #     y_close = np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)
            #     if entry == True and (not x_close or not y_close):
            #         valid = False
            #         if act_ent > req_ent:
            #             bank-=2.5
            #         else:
            #             bank+=10
            #     elif entry == True and x_close and y_close:
            #         # print(path_toa+new_total, pb_toa, new_total)
            #         all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, ba, entryTOA, act_ent, gi[0], end, entry])
            #         all_nodes.append(nodes1)
            #         valid = True
            #     else:
            #         valid = True
            #         all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, gi[0], end, entry])
            #         all_nodes.append(nodes1)
                
            # else:
            #     gi[1]+=0.25
            #     if path_toa+new_total > target: 
            #         g-=500
            #     else:
            #         g+=500
            #     if gi[1] > 1:
            #         # print(f'PATH WAS INVALID, APPENDING, {path_toa+new_total}, {path_toa}, {new_total}')
            #         valid = True
            #         if entry:
            #             all_bez.append([path_toa+new_total, partial_bez1[0], partial_bez1[1], new_total, x_entry, y_entry, pb_toa, pb_length, entryLength, ba, entryTOA, act_ent, gi[0], end, entry])
            #         else:
            #             all_bez.append([path_toa+new_total, optim_b[0], optim_b[1], new_total,  path_toa, optim_b_length, gi[0], end, entry])
            #         all_nodes.append(nodes1)

        gi[0]+=1


    # gi[0]-=1
    stop_time = time.time()
    print('TOTAL JUST OPTIMIZATION TIME', np.sum(opt))
    diffs = []
    # print(all_nodes)
    for i in all_bez:
        diffs.append(np.abs(target-i[0]))
        print(np.abs(target-i[0]), i[len(i)-2])
    
    mindex = diffs.index(np.min(diffs))
    # print(mindex)
    # print(all_bez[mindex])
    optim_b = [all_bez[mindex][1], all_bez[mindex][2]]
    nodes1 = all_nodes[mindex]
    # print(len(all_nodes), all_nodes)
    # plt.plot(optim_b[0], optim_b[1])
    # plt.scatter(nodes1[0], nodes1[1])
    # plt.show()
    print(f'TOTAL OPTIMIZATION TIME: {stop_time-start_time}')
    
    # print("ENTRY TOA:", entryTOA)
    # print("Bez TOA:", pb_toa)
    # print(f'Target Total is {target}\nTarget Bez TOA is: {t_t}\nBez+Entry Path TOA is: {path_toa}\nTotal Travel Time w/ STAR Route is {path_toa+new_total}\nValidity is {valid}')
    # # print(curv1, turn_radius/.3048)
    # print(f'STAR Path Intercept is {gi}')
    gi[0] = all_bez[mindex][len(all_bez[mindex])-3]
    # print(all_bez[mindex][0])
    # print(all_bez[0][0], all_bez[1][0], all_bez[2][0])
    end_point = get_line(wps[gi[0]], wps[gi[0]+1], all_bez[mindex][len(all_bez[mindex])-1])
    new_travel, new_total = getTravelXY(wps, end_point, gi[0], toa_local, v)
    # print(new_total)
    
    rotated_bezier = rotate_bez(np.array(optim_b), -angle_to_rotate, origin)
    rot_nodes = rotate_bez(nodes1, -angle_to_rotate, origin)
    # print(rot_nodes)
    try:
        ep = get_line(wps[gi[0]], wps[gi[0]+1], gi[1])
        m2 = (wps[gi[0]][1] - ep[1])/(wps[gi[0]][0] - ep[0])
    except:
        ep = get_line(wps[gi[0]], wps[gi[0]+1], gi[1])
        t_h = np.arctan2(wps[gi[0]+1][1] - wps[gi[0]][1],wps[gi[0]+1][0] - wps[gi[0]][0])
        m2 = (wps[gi[0]][1] - (ep[1]+lr*1000*np.abs(np.sin(t_h))))/(wps[gi[0]][0] - (ep[0]+lr*1000*np.abs(np.cos(t_h))))
    b2 = wps[gi[0]][1] - m2*wps[gi[0]][0]
    # plt.plot(rotated_entry[0], rotated_entry[1])
    # plt.plot(rotated_bezier[0], rotated_bezier[1])
    show_local = False
    print(f'TIME DIFFERENCE OF ARRIVAL: {np.abs(all_bez[mindex][0] - target)}')
    plt.plot(x, y, label = 'STAR Path', linewidth = 2)
    if all_bez[mindex][-1]:

        rotated_partial = rotate_bez(np.array([all_bez[mindex][1], all_bez[mindex][2]]), -angle_to_rotate, origin)
        rotated_entry = rotate_bez(np.array([all_bez[mindex][4], all_bez[mindex][5]]), -angle_to_rotate, origin)
        x_entry = all_bez[mindex][4]
        y_entry = all_bez[mindex][5]
        print(f'Entry Bank{np.rad2deg(all_bez[mindex][9][0])}, Entry Point On Bez {all_bez[mindex][9][1]}')
        print('REQUIRED HEADING AT ENTRY:', all_bez[mindex][11][1]-angle_to_rotate-360)
        print('ACTUAL HEADING AT INTERSECT:', all_bez[mindex][11][0]-angle_to_rotate-360)
        print('DIFF AT ENTRY:', all_bez[mindex][11][1]-all_bez[mindex][11][0])
        print('BEZ CURVATURE:', all_bez[mindex][12], 'TURN RADIUS:', turn_radius)
        if show_local:
            '''
            Not Rotated -> Local frame
            '''
            # plt.plot(x_entry, y_entry, label = 'Entry Arc', color = 'cyan')
            # plt.scatter(x_entry[-1], y_entry[-1], marker = '*', color = 'yellow', label = 'Intersection Point', s = 100, zorder = 100)
            # plt.text(x_entry[-1]+1000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
            # plt.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial Bez Path')
            # plt.text(partial_bez1[0][50], partial_bez1[1][50]-2500, f'Partial Bez Travel Time: {pb_toa:0.2f}s', weight = 'bold' )


        '''
        'Rotated' -> planned trajectory in global frame
        '''
        plt.plot(rotated_entry[0], rotated_entry[1], label = 'Entry Arc', color = 'cyan', linewidth = 2)
        plt.scatter(rotated_entry[0][-1], rotated_entry[1][-1], marker = '*', color = 'yellow', label = 'Intersection Point', s = 100, zorder = 100)
        plt.text(rotated_entry[0][0]+5000,rotated_entry[1][0]-2500, f'Entry Path Travel Time: \n{all_bez[mindex][10]:0.2f}', weight = 'bold')
        plt.plot(rotated_partial[0], rotated_partial[1], c = 'magenta', label = 'Partial Bez Traj', linewidth = 2)
        plt.text(rotated_bezier[0][10], rotated_bezier[1][10]-5000, f'Partial Bez Travel Time: \n{all_bez[mindex][6]:0.2f}s', weight = 'bold' )
        plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--', linewidth = 2)
    else:
        plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--', linewidth = 2, label = 'Bez Traj')
        print(f'No Entry Path Used\nRequired Bez Heading: {np.rad2deg(np.arctan2(rotated_bezier[1][0] - rotated_bezier[1][1], rotated_bezier[0][0] - rotated_bezier[0][1]))}')
        print(f'Actual Heading: {np.rad2deg(head)}')
        print('BEZ CURVATURE:', all_bez[mindex][6], 'TURN RADIUS:', turn_radius)
        # plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        plt.text(rotated_bezier[0][25], rotated_bezier[1][25]-2500, f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
    # if head > pi == 'IV':
    #     plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head+np.pi), 2500*np.sin(head+np.pi), width=50,
    #                 length_includes_head=False,
    #                 head_width=1000,
    #                 head_length=1000,
    #                 head_starts_at_zero=False,
    #                 facecolor='black',
    #                 edgecolor='black',
    #                 label = f'Heading of {np.rad2deg(head)}')
    # else:
    plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head), 2500*np.sin(head), width=50,
                length_includes_head=False,
                head_width=1000,
                head_length=1000,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                label = f'Heading of {np.rad2deg(head):0.2f}')
    if show_local:
        '''
        Not Rotated -> Local frame 
        '''
        # for i in range(len(wps)-1):
        #     plt.scatter(interps_local[i][0], interps_local[i][1], color = 'orange')

        # ang_y = (nodes[0][2] - 5000) *m + b
        # plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
        # plt.plot(bez[0], bez[1], color = 'black', linestyle = '--')
        # plt.plot(optim_b[0], optim_b[1], color = 'black', linestyle = '--')
        # plt.scatter(optim1[0], optim1[1], marker = '^', color = 'red')
        # plt.scatter(nodes[0][2], nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
        # plt.scatter(path_startX_local, path_startY_local, marker = '^', color = 'red')
        # plt.scatter(x_rot, y_rot, label = 'STAR Path')
        # plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

    '''
    'Rotated' -> planned trajectory in global frame
    '''
    for i in range(len(wps)-1):
        plt.scatter(wps[i][0], wps[i][1], color = 'orange', marker = 's')
    ang_y = (rotated_bezier[0][-1] + 5000) *m2 + b2
    plt.plot([rotated_bezier[0][-1],rotated_bezier[0][-1] + 5000], [rotated_bezier[1][-1], ang_y], color = 'green', linestyle = '--')
    plt.scatter(rot_nodes[0][1], rot_nodes[1][1], marker = '^', color = 'red')
    # plt.scatter(rot_nodes[0][2], rot_nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
    plt.scatter(rotated_bezier[0][-1], rotated_bezier[1][-1], marker = '^', color = 'red', label = 'Control Points')

    plt.scatter(path_startX, path_startY, marker = '^', color = 'red')
    
    # plt.text(wps[0][0], wps[0][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    # plt.text(wps[1][0], wps[1][1], f'STAR Travel Time: {all_bez[mindex][3]:0.2f}s', weight = 'bold' )
    # plt.text(wps[-1][0], wps[-1][1], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    plt.text(-20000, 5000, f'Required Travel Time: {target:0.2f}\nActual Travel Time: {all_bez[mindex][0]:0.2f}', weight = 'bold')
    plt.text(wps[2][0]+1000, wps[2][1], f'STAR Travel Time: \n{all_bez[mindex][3]:0.2f}s', weight = 'bold' )
    # nodes_r = rotate_bez(all_nodes[4], -angle_to_rotate, origin)
    # new_b = manual_bez([nodes_r[0][0], nodes_r[1][0]], [nodes_r[0][1], nodes_r[1][1]], [nodes_r[0][2], nodes_r[1][2]], 200)
    # plt.plot(new_b[0], new_b[1], linewidth = 3, color = 'red')
    plt.scatter(ac_met[0], ac_met[1], color = 'green')
    plt.scatter(0, 0, marker = '*', color = 'green', label = 'CMH')
    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    plt.grid()
    plt.title(f'Quadrant {gate[0]}, Path Optimization Time: {(stop_time-start_time):0.2f}')
    plt.legend(loc = 'upper left')
    plt.axis('equal')
    plt.show()
