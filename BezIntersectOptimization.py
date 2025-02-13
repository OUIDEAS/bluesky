import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy.io
import sympy as sp
from scipy.optimize import root
from scipy.optimize import brentq
from scipy.optimize import fsolve
from pyproj import Proj
from utm import from_latlon as llutm
from utm import to_latlon as utmll

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
    def path_cost(ba, nodes, velocity, lr, h, show):
        diff = find_diff_entry(ba, nodes, velocity, pos, lr, h, show)
        # if diff[1] == 0.01:  # No valid intersection
        #     return 1e6
        return np.abs(diff[0])
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}  
    )
    # val= minimize(path_cost,guess,(nodes, velocity, lr, h, show), method='SLSQP', tol=1E-3, constraints=cons)
    # print(f"Optimization result: Bank angle = {np.rad2deg(val.x[0])} degrees")
    val = minimize(path_cost, guess, args=(nodes, velocity, lr, h, show),
                      method='trust-constr', bounds=[(min_bank, max_bank)], constraints=cons,
                      options={'xtol': 1e-3, 'gtol': 1e-3, 'maxiter': 100})

    t_val = find_diff_entry(val.x, nodes, velocity, pos, lr, h, show)[1]
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

def find_diff_entry(ba, nodes, velocity, pos, lr, head, show):
    pi = np.pi
    t_guess = 0.01
    mindiff = 1e6
    # print(np.rad2deg(head))
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
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
    dist_to_bezier = np.sqrt((Bx(0) - h) ** 2 + (By(0) - k) ** 2)

    # Use this distance to get an initial guess closer to the intersection
    # t_guess = np.clip(dist_to_bezier / 2, 0.1, 0.9)
    S = fsolve(circle_eq, t_guess)
    # result = root(circle_eq, t_guess, method = 'hybr')
    # print(S)
    t_final = 0.01

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
                y = np.linspace(pos[1], By(i), 200)
                x_l = h-np.sqrt(tr**2 - (y-k)**2)
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
                plt.plot(x_l, y)
                plt.axis('equal')
                plt.show()

            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            # print(f'INT ANGLE {np.rad2deg(int_angle)}')
            if math.isnan(int_angle):
                int_angle = np.arctan2(y[198]-y[197], x_l[198] - x_l[197])
            if np.deg2rad(int_angle) <0:
                int_angle+=2*pi
            # print(f'BEZ ANGLE {np.rad2deg(bez_angle)}, INT ANGLE {np.rad2deg(int_angle)}, DIFFERENCE {np.abs(bez_angle-int_angle)}')
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

    # if not result.success:
    #     print("Root finding failed")
    #     return [1e6, 0.01]  # Return a large penalty and a failed solution

    # # Process valid intersection
    # t_final = result.x[0]
    # if not (0 <= t_final <= 1):
    #     return [1e6, 0.01]  # No valid intersection found

    # # Calculate corresponding values
    # index = int(np.round(t_final * 200))
    # if index == 200:
    #     index = 199
    # if index == 0:
    #     index = 5
    # bez_angle = np.arctan2(By(t_final + 0.01) - By(t_final), Bx(t_final + 0.01) - Bx(t_final))

    # if lr == 1:
    #     y = np.linspace(pos[1], By(t_final), 200)
    #     x_l = h - np.sqrt(tr ** 2 - (y - k) ** 2)
    # else:
    #     x_l = np.linspace(pos[0], Bx(t_final), 200)
    #     y = k + np.sqrt(tr ** 2 - (x_l - h) ** 2)

    # int_angle = np.arctan2(y[-1] - y[198], x_l[-1] - x_l[198])
    # diff = np.abs(bez_angle - int_angle)

    # if diff < mindiff and np.abs(y[-1] - By(t_final)) <= 5 and np.abs(x_l[-1] - Bx(t_final)) <= 5:
    #     mindiff = diff
    #     t_final = t_final

    # else:
    #     dist_to_curve = np.sqrt((Bx(t_final) - x_l[-1]) ** 2 + (By(t_final) - y[-1]) ** 2)
    #     mindiff = max(mindiff, dist_to_curve * 100)  # Apply penalty

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

    tr = 111.6**2 / (11.26*math.tan(ba))
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
    print(f'BANK ANGLE: {ba}')
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

    # if head == 0 or head == 2*np.pi:
    #     h = pos[0] - lr*tr*np.sin(head)
    #     k = pos[1] - lr*tr*np.cos(head)
    # if lr == 1:
    #             y = [j for j in  np.linspace(pos[1], By(i), 200)]
    #             x_l = [h+np.sqrt(tr**2 - (j-k)**2) for j in y]
    #         else:
    #             x_l = [j for j in np.linspace(pos[0], Bx(i), 200)] #Space between nominal path and bez
    #         # print(x_l)
            
    #             y = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
    if lr == 1:
        y_entry = list(np.linspace(pos[1], intersect[1], 200))
        x_entry = [h-lr*np.sqrt(tr**2 - (y-k)**2) for y in y_entry]
    else:
        x_entry = list(np.linspace(pos[0], intersect[0], 200))
        y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
    if math.isnan(y_entry[0]):
        x_entry = list(np.linspace(h-lr*tr, intersect[0], 200))
        y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
        # print(x_entry)
        # print(y_entry)
        if math.isnan(y_entry[0]):
            y_entry[0] = pos[1]
            x_entry[0] = pos[0]
            # print(x_entry)
            # print(y_entry)]
    if math.isnan(y_entry[-1]):
        y_entry = list(np.linspace(pos[1], intersect[1], 200))
        x_entry = [h-lr*np.sqrt(tr**2 - (y-k)**2) for y in y_entry]
    for i in range(0,len(y_entry)):
        if i<len(y_entry) and math.isnan(y_entry[i]):
            y_entry.pop(i)
            x_entry.pop(i)
        elif i<len(y_entry) and math.isnan(x_entry[i]):
            y_entry.pop(i)
            x_entry.pop(i)
            # ar, ad = central_angle([h, k], [x_entry[0], y_entry[0]], [x_entry[-1], y_entry[-1]])
            # entryLength = tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
            # entryTOA = entryLength/velocity
            # print(entryTOA)
            # y_entry.pop(0)
            # x_entry.pop(0)
            # ar, ad = central_angle([h, k], [x_entry[0], y_entry[0]], [x_entry[-1], y_entry[-1]])
            # entryLength = tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
            # entryTOA = entryLength/velocity
            # print(entryTOA)
    # y_entry[-1] = intersect[1]
    ar, ad = central_angle([h, k], [x_entry[0], y_entry[0]], [x_entry[-1], y_entry[-1]])
    entryLength = tr*ar #entryLength = 2*pi*tr * (central_angle/(2*pi))
    entryTOA = entryLength/velocity
    # print('ENTRY ToA: ', entryTOA)

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
        if gi[0] > i:
            new_time.append(0)
            total+=0
        elif gi[0] == i and i+1 < len(wps):
            dx = point[0] - wps[i+1][0]
            dy = point[1] - wps[i+1][1]
            d = np.hypoy(dx, dy)
            new_time.append(d/v)
            total+=d/v
        elif gi[0] == i and i+1>= len(wps):
            dx = point[0] - wps[i][0]
            dy = point[1] - wps[i][1]
            d = np.hypot(dx, dy)
            new_time.append(d/v)
            total+=d/v
        elif gi[0] < i:
            new_time.append(travel[i])
            total+=travel[i]

def total_travel(waypts, v):
    # interps = np.zeros((len(waypts)-1, 2, n))
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
        # x = np.linspace(waypts[i][0], waypts[i+1][0], n, endpoint=False)
        # interps[i][0][0], interps[i][1][0] = x[0], waypts[i][1]
        # for j in range(1, len(x)):
        #     y = waypts[i][1] + m * (x[j]-waypts[i][0])
        #     interps[i][0][j] = x[j]
        #     interps[i][1][j] = y
    print(f'Total STAR Route Travel Time: {total}')
    return dists, travel, total

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
    # m = (p2[1]-p1[1])/(p2[0]-p1[0])
    x = p1[0] + (p2[0]-p1[0])*t
    y = p1[1] + (p2[1]-p1[1])*t
    return [x,y]

def GenerateTrajectoryOutside(ac_pos, target, wp_list, heading, velocity, kcmh, gate):
    v=velocity
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    print(heading, gate)
    if gate == 'II':
        heading = heading+180
        head = np.deg2rad(heading)
        angle_to_rotate = (90-heading)
    elif gate == 'III':
        heading=heading
        head = np.deg2rad(heading)
        angle_to_rotate = (90-heading)
    elif gate == 'IV':
        heading = heading-180
        head = np.deg2rad(heading)
        angle_to_rotate = (90-heading)
    elif gate == 'I':
        heading=270-heading
        head = np.deg2rad(heading-180)
        angle_to_rotate = (90-heading+180)
    rotated_heading = np.deg2rad(heading+angle_to_rotate)
    plt.scatter(ac_pos[1], ac_pos[0])
    # print(transform_angle)
    
    wps = LongLat_To_WSG84_Meters(wp_list, kcmh) #Gives waypoints in global frame
    print(wps)
    
    ac_met = __WSG84_To_Meters_Single(ac_pos, kcmh) #Global aircraft position

    wps_local = LongLat_To_WSG84_Meters(wp_list, ac_pos)
    # ac_local = __WSG84_To_Meters_Single(ac_pos, ac_pos)
    # he, de = qdrdist(ac_pos[0], ac_pos[1], dubln[0], dubln[1])
    # print(de/v)

    start_time = time.time()

    x = [i[0] for i in wps]
    y = [i[1] for i in wps]

    # x_local = [i[0] for i in wps_local]
    # y_local = [i[1] for i in wps_local]



    # interp_pts = 100
    # interps, dists, toa, total = interp_wps(wps, interp_pts, velocity) #Interpolate points between waypoints
    for i in range(len(wps)-1):
        plt.scatter(interps[i][0], interps[i][1], color = 'orange', marker = 's')
    plt.show()
    origin = np.array([ac_met])
    # print(origin)

    rotated_wps = rotate_points(wps, angle_to_rotate, origin)
    interps_local, dists_local, toa_local, total_local = interp_wps(rotated_wps,  interp_pts, velocity)

    x_rot = [i[0] for i in rotated_wps]
    y_rot = [i[1] for i in rotated_wps]

    d_reg = np.hypot(x[0]-ac_met[0], y[0]-ac_met[1])
    # d_loc = np.hypot(x_local[0]-ac_local[0], y_local[0]-ac_local[1])
    d_rot = np.hypot(x_rot[0]-ac_met[0], y_rot[0]-ac_met[1])
    print(f'REGULAR DISTANCE TO FIRST WP: {d_reg}, ROTATED: {d_rot}')

    if interps_local[1][0][0] < ac_met[0]:
        lr = -1
        # g = 5
        # min = 2
        # plan_ahead = 60
    else:
        lr = 1
    bank = 5
    min = 5
    plan_ahead = 10
    
    path_startX = ac_met[0] + np.cos(head)*v*plan_ahead 
    path_startY = ac_met[1] + np.sin(head)*v*plan_ahead 

    path_startX_local = ac_met[0]
    path_startY_local = ac_met[1] + v*plan_ahead

    
    # nodes = [np.array([path_startX, (path_startX + interps[0][0][25])/2, interps[0][0][25]]).flatten(),
    #          np.array([path_startY, (path_startY + interps[0][1][25])/1.5, interps[0][1][25]]).flatten()]
    
    # nodes_local = [np.array([path_startX_local, (path_startX_local + interps_local[0][0][25])/2, interps_local[0][0][25]]).flatten(),
    #          np.array([path_startY, (path_startY + interps_local[0][1][25])/1.5, interps_local[0][1][25]]).flatten()]
    # # print(nodes[0], nodes[0][1])
    # # print(interps[0][1][25])
    turn_radius = v*1.94384**2/(11.26*math.tan(np.deg2rad(73)))
    turn_radius*=0.3048
    
    # # print(interps_local)
    # bez = manual_bez([nodes[0][0],nodes[1][0]], [nodes[0][1], nodes[1][1]], [nodes[0][2], nodes[1][2]], 200)

    # bez_local = manual_bez([nodes_local[0][0],nodes_local[1][0]], [nodes_local[0][1], nodes_local[1][1]], [nodes_local[0][2], nodes_local[1][2]], 200)
    valid = False
    gi = [0,2]
    t_t = target-total-plan_ahead
    new_travel, new_total = getTravelTime(interps, gi, toa, v)

    new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
    t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
    # print(new_total)
    t_t = target-new_total-plan_ahead
    
    entry = False
    entryTOA = plan_ahead
    while not valid:
        show = False
        new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
        t_t = target-new_total-entryTOA
        
        nodes = [np.array([path_startX_local, interps_local[gi[0]][0][gi[1]]+lr*2500*np.abs(np.cos(t_h)), interps_local[gi[0]][0][gi[1]]]).flatten(),
             np.array([path_startY_local, interps_local[gi[0]][1][gi[1]]+lr*2500*np.abs(np.sin(t_h)), interps_local[gi[0]][1][gi[1]]]).flatten()]
       
        m = (interps_local[gi[0]][1][0] - interps_local[gi[0]][1][gi[1]])/(interps_local[gi[0]][0][0] - interps_local[gi[0]][0][gi[1]])
        b = interps_local[gi[0]][1][0] - m*interps_local[gi[0]][0][0]
        coeffs = [m, b]
        
        t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
        
        optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                    target_toa=t_t,
                                    guess = [nodes[0][1], nodes[1][1]],
                                    target_heading=t_h, velocity=velocity, turn_radius=turn_radius, lr = lr, line = coeffs)
        if optim1[0] < ac_met[0]:
            lr = -1
        else:
            lr = 1
        optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
        optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)
        path_toa = optim_b_length/v
        
        nodes1 = [np.array([nodes[0][0], optim1[0], nodes[0][2]]).flatten(), np.array([nodes[1][0], optim1[1], nodes[1][2]])]
        
        if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(10)):
            entryTOA = plan_ahead
            path_toa = optim_b_length/v + entryTOA
            entry = False
            
        else:
            entry = True
            bez_ang = np.rad2deg(np.arctan2(optim_b[1][2]-optim_b[1][0], optim_b[0][2]-optim_b[0][0]))
            if bez_ang<0:
                bez_ang+=360
            if bez_ang>=270:
                bank = 30
            ba, t_entry = solve_optimEntry(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotated_heading, show)
            
            x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [nodes1[0][2],nodes1[1][2]], t_entry)
            
            x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
            req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
            
            
            x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotated_heading)
            act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
            

            partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                    [optim1[0],optim1[1]],
                                    [nodes1[0][2],nodes1[1][2]], 200, t_entry)
            pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],t_entry)
            pb_toa = pb_length/v
            path_toa = entryTOA + pb_toa

        ang_y = (nodes[0][2] - 5000) *m + b
        '''
        Plots for step by step
        '''
        # plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
        # plt.plot(x_rot, y_rot)
        # plt.scatter(ac_met[0], ac_met[1])
        # plt.plot(x_rot, y_rot)
        # plt.scatter(optim1[0], optim1[1])
        # plt.plot(optim_b[0], optim_b[1], color = 'black')
        # plt.axis('equal')
        # plt.text(interps_local[-1][0][0], interps_local[-1][1][0], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
        # plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

        # if entry:
        #     plt.plot(x_entry, y_entry, color = 'cyan')
        #     plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
        #     plt.text(optim_b[0][50], optim_b[1][50]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
        # #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
        # else:
        # #     print(f'Bez Travel Time: {path_toa:0.2f}s')
        #     plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        # plt.show()
        
        
        
        print(new_total + path_toa)
        print(f'TARGET: {target}, ACTUAL TRAVEL TIME: {path_toa+new_total}')
        print(np.abs(target - (path_toa+new_total)))
        print(f'Original Target: {target}, New STAR Travel Time: {new_total}, New Target Time: {t_t}, Path TOA: {path_toa}')
        print(new_total+path_toa, target)
        if 1<= np.abs(target - (path_toa+new_total)) <= 2.5:
            g = 5000
            show = False
            for k in range(0, 5):
                new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
                t_t = target-new_total-entryTOA
                
                nodes = [np.array([path_startX_local, interps_local[gi[0]][0][0] + lr*g*np.abs(np.cos(t_h)), interps_local[gi[0]][0][gi[1]]]).flatten(),
                    np.array([path_startY_local, interps_local[gi[0]][1][0]+lr*g*np.abs(np.sin(t_h)), interps_local[gi[0]][1][gi[1]]]).flatten()]
                
                m = (interps_local[gi[0]][1][0] - interps_local[gi[0]][1][gi[1]])/(interps_local[gi[0]][0][0] - interps_local[gi[0]][0][gi[1]])
                b = interps_local[gi[0]][1][0] - m*interps_local[gi[0]][0][0]
                coeffs = [m, b]
                
                t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
                
                optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                            target_toa=t_t,
                                            guess = [nodes[0][1], nodes[1][1]],
                                            target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs)
                if optim1[0] < ac_met[0]:
                    lr = -1
                else:
                    lr = 1
                optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
                optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)
                path_toa = optim_b_length/v
                
                nodes1 = [np.array([nodes[0][0], optim1[0], nodes[0][2]]).flatten(), np.array([nodes[1][0], optim1[1], nodes[1][2]])]
               
                if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(10)):
                    entryTOA = plan_ahead
                    path_toa = optim_b_length/v + entryTOA
                    entry = False
                    # print('using plan ahead TOA')
                else:
                    entry = True
                    bez_ang = np.rad2deg(np.arctan2(optim_b[1][2]-optim_b[1][0], optim_b[0][2]-optim_b[0][0]))
                    if bez_ang<0:
                        bez_ang+=360
                    if bez_ang>=270:
                        bank = 30
                    ba, t_entry = solve_optimEntry(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotated_heading, show)
                    
                    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                                [optim1[0],optim1[1]],
                                                [nodes1[0][2],nodes1[1][2]], t_entry)
                    
                    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                                [optim1[0],optim1[1]],
                                                [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
                    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
                    
                    
                    x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotated_heading)
                    act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
                    

                    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [nodes1[0][2],nodes1[1][2]], 200, t_entry)
                    pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],t_entry)
                    pb_toa = pb_length/v
                    path_toa = entryTOA + pb_toa
                ang_y = (nodes[0][2] - 5000) *m + b
                '''
                Plots for step by step
                '''
                # plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
                # plt.plot(x_rot, y_rot)
                # plt.scatter(ac_met[0], ac_met[1])
                # plt.plot(x_rot, y_rot)
                # plt.scatter(optim1[0], optim1[1])
                # plt.plot(optim_b[0], optim_b[1], color = 'black')
                # plt.axis('equal')
                # plt.text(interps_local[-1][0][0], interps_local[-1][1][0], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
                # plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

                # if entry:
                #     plt.plot(x_entry, y_entry, color = 'cyan')
                #     plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
                #     plt.text(optim_b[0][50], optim_b[1][50]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
                # #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
                # else:
                # #     print(f'Bez Travel Time: {path_toa:0.2f}s')
                #     plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
                # plt.show()
                g+=2500
                if bank > min:
                    bank -= 2.5
                elif bank <= min:
                    bank = min
                if np.abs(target - (path_toa+new_total))<=1:
                    print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
                    print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
                    if entry == True and (not np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5) or not np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)):
                        valid = False
                    else:
                        valid = True
                        break

        if np.abs(target - (path_toa+new_total))<=1:
            print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
            print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
            if entry == True and (not np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5) or not np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)):
                valid = False
            else:
                valid = True
        

        if valid == False:
            # print(f'Target TOA is: {target-total-plan_ahead}\nPath TOA is: {path_toa}\n Validity is {valid}')
            gi[1]+=2
            if gi[1] >= interp_pts/2 or gi[0]>=2:
                gi[1]+=3
            if gi[1]>=interp_pts:
                gi[0]+=1
                gi[1] = 0

            if gi[0] >= len(interps):
                gi = [0, 0]
                nodes = [np.array([ac_met[0], (ac_met[0] + interps[gi[0]][0][gi[1]])/2, interps[gi[0]][0][gi[1]]]).flatten(),
                np.array([ac_met[1], (ac_met[1] + interps[gi[0]][1][gi[1]])/2+10, interps[gi[0]][1][gi[1]]]).flatten()]
                m = (interps[gi[0]][1][10] - interps[gi[0]][1][0])/(interps[gi[0]][0][10] - interps[gi[0]][0][0])
                b = interps[gi[0]][1][0] - m*interps[gi[0]][0][0]
                coeffs = [m, b]
                t_h = np.arctan2(interps[gi[0]][1][10]- interps[gi[0]][1][0],interps[gi[0]][0][10]- interps[gi[0]][0][0] )
                optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                            target_toa=target-total,
                                            guess = [nodes[0][1], nodes[1][1]],
                                            target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = 1, line = coeffs)
                t_t= target-total
                optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
                path_toa = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)/v
                break

    stop_time = time.time()
    print(f'TOTAL OPTIMIZATION TIME: {stop_time-start_time}')
    
    print("ENTRY TOA:", entryTOA)
    print("Bez TOA:", pb_toa)
    print(f'Target Total is {target}\nTarget Bez TOA is: {t_t}\nBez+Entry Path TOA is: {path_toa}\nTotal Travel Time w/ STAR Route is {path_toa+new_total}\nValidity is {valid}')
    # print(curv1, turn_radius/.3048)
    print(f'STAR Path Intercept is {gi}')

    

    rotated_bezier = rotate_bez(np.array(optim_b), -angle_to_rotate, origin)
    rotated_partial = rotate_bez(np.array(partial_bez1), -angle_to_rotate, origin)
    rotated_entry = rotate_bez(np.array([x_entry, y_entry]), -angle_to_rotate, origin)
    rot_nodes = rotate_bez(nodes1, -angle_to_rotate, origin)
    m2 = (interps[gi[0]][1][0] - interps[gi[0]][1][gi[1]])/(interps[gi[0]][0][0] - interps[gi[0]][0][gi[1]])
    b2 = interps[gi[0]][1][0] - m2*interps[gi[0]][0][0]
    # plt.plot(rotated_entry[0], rotated_entry[1])
    # plt.plot(rotated_bezier[0], rotated_bezier[1])
    show_local = False

    if entry:
        print(f'Entry Bank{np.rad2deg(ba)}, Entry Point On Bez {t_entry}')
        print('REQUIRED HEADING AT ENTRY:', req_ent)
        print('ACTUAL HEADING AT INTERSECT:', act_ent)
            
        print('DIFF AT ENTRY:', req_ent-act_ent)
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
        plt.plot(rotated_entry[0], rotated_entry[1], label = 'Entry Arc', color = 'cyan')
        plt.scatter(rotated_entry[0][-1], rotated_entry[1][-1], marker = '*', color = 'yellow', label = 'Intersection Point', s = 100, zorder = 100)
        plt.text(x_entry[-1]+10000, y_entry[-1]-5000, f'Entry Path Travel Time: {entryTOA:0.2f}', weight = 'bold')
        plt.plot(rotated_partial[0], rotated_partial[1], c = 'magenta', label = 'Partial Bez Path')
        plt.text(rotated_bezier[0][50], rotated_bezier[1][50]-2500, f'Partial Bez Travel Time: {pb_toa:0.2f}s', weight = 'bold' )
    else:
        print(f'No Entry Path Used\nRequired Bez Heading: {np.rad2deg(np.arctan2(rotated_bezier[1][0] - rotated_bezier[1][1], rotated_bezier[0][0] - rotated_bezier[0][1]))+180}')
        print(f'Actual Heading: {np.rad2deg(head)}')
        # plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        plt.text(rotated_bezier[0][50], rotated_bezier[1][50]-2500, f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )

    plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head), 2500*np.sin(head), width=50,
                length_includes_head=False,
                head_width=1000,
                head_length=1000,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                label = f'Heading of {np.rad2deg(head)}')
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
        plt.scatter(interps[i][0], interps[i][1], color = 'orange', marker = 's')
    ang_y = (rot_nodes[0][2] - 5000) *m2 + b2
    plt.plot([rot_nodes[0][2],rot_nodes[0][2] - 5000], [rot_nodes[1][2], ang_y], color = 'green', linestyle = '--')
    plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--')
    plt.scatter(rot_nodes[0][1], rot_nodes[1][1], marker = '^', color = 'red')
    plt.scatter(rot_nodes[0][2], rot_nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
    plt.scatter(path_startX, path_startY, marker = '^', color = 'red')
    plt.scatter(x, y, label = 'STAR Path')
    plt.text(-50000, 35000, f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
    plt.text(interps[1][0][0], interps[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

    plt.scatter(ac_met[0], ac_met[1], color = 'green')
    plt.scatter(0, 0, marker = '*', color = 'green', label = 'CMH')
    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    plt.grid()
    plt.legend(loc = 'lower left')
    plt.axis('equal')
    plt.show()
    return interps, rot_nodes, rotated_bezier, [rotated_entry, req_ent-2*np.pi, np.rad2deg(ba)]



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
    ac_pos = [40.303627035156936 , -83.37719169326613] #g2
    h =131
    ac_pos = [39.58268448211492, -83.14955754333681] #g3
    h=27
    ac_pos = [39.78201288955758, -82.35136330299837] #g4
    h=299
    # head = np.deg2rad(h)
    # angle_to_rotate = (90-h)
    # rotate_head = np.deg2rad(h+angle_to_rotate)

    # heading = h
    # head = np.deg2rad(heading)
    # angle_to_rotate = (90-heading)
    heading = h-180
    head = np.deg2rad(heading)
    angle_to_rotate = (90-heading)
    rotate_head = np.deg2rad(heading+angle_to_rotate)

    print(f'HEADING: {h}, ROTATE ANGLE: {angle_to_rotate}, ROTATED HEADING: {h+angle_to_rotate}')
    
    # head = np.deg2rad(h)
    # h=270-h
    # angle_to_rotate = (90-h)
    # print(transform_angle)

    wps = LongLat_To_WSG84_Meters(xavyr_path, kcmh) #Gives waypoints in global frame
    print(wps)
    ac_met = __WSG84_To_Meters_Single(ac_pos, kcmh) #Global aircraft position

    wps_local = LongLat_To_WSG84_Meters(xavyr_path, ac_pos)
    print('WAYPOINTS:',wps)
    print('LOCAL WAYPOINTS:',wps_local)
    # ac_local = __WSG84_To_Meters_Single(ac_pos, ac_pos)
    # he, de = qdrdist(ac_pos[0], ac_pos[1], dubln[0], dubln[1])
    # print(de/v)

    start_time = time.time()

    x = [i[0] for i in wps]
    y = [i[1] for i in wps]

    # x_local = [i[0] for i in wps_local]
    # y_local = [i[1] for i in wps_local]



    interp_pts = 100
    interps, dists, toa, total = interp_wps(wps, interp_pts, v) #Interpolate points between waypoints

    origin = np.array([ac_met])
    print('ORIGIN', origin)

    rotated_wps = rotate_points(wps, angle_to_rotate, origin)
    print('ROTATED WAYPOINTS:',rotated_wps)

    # interps_local, dists_local, toa_local, total_local = interp_wps(rotated_wps,  interp_pts, v)

    x_rot = [i[0] for i in rotated_wps]
    y_rot = [i[1] for i in rotated_wps]

    d_reg = np.hypot(x[0]-ac_met[0], y[0]-ac_met[1])
    # d_loc = np.hypot(x_local[0]-ac_local[0], y_local[0]-ac_local[1])
    d_rot = np.hypot(x_rot[0]-ac_met[0], y_rot[0]-ac_met[1])
    print(f'REGULAR DISTANCE TO FIRST WP: {d_reg}, ROTATED: {d_rot}')
    # interps_local = rotate_interpolated_points(interps, angle_to_rotate, origin)
    # rotated_ac = rotate_points(np.array([ac_pos]), angle_to_rotate, origin)[0]
    # rotated_bezier = rotate_points(np.array(bez_path), angle_to_rotate, origin)
    # rotated_star = rotate_points(np.array(star_path), angle_to_rotate, origin)
    interps_local, dists_local, toa_local, total_local = interp_wps(wps_local, interp_pts, v)
    dists_local, toa_local, total_local = total_travel(wps_local, v)
    # print(dists, dists_local)
    # plt.plot(x_rot, y_rot)
    for i in range(len(wps)-1):
        plt.plot([rotated_wps[i][0], rotated_wps[i+1][0]], [rotated_wps[i][1], rotated_wps[i+1][1]], color = 'orange')
        plt.plot([wps[i][0], wps[i+1][0]], [wps[i][1], wps[i+1][1]])
    # plt.plot(x_local, y_local)
    plt.scatter(origin[0][0], origin[0][1])
    plt.show()
    if rotated_wps[1][0] < ac_met[0]:
        lr = -1
        # g = 5
        # min = 2
        # plan_ahead = 60
    else:
        lr = 1
    bank = 5
    min = 5
    plan_ahead = 10
    
    path_startX = ac_met[0] + np.cos(head)*v*plan_ahead 
    path_startY = ac_met[1] + np.sin(head)*v*plan_ahead 

    path_startX_local = ac_met[0]
    path_startY_local = ac_met[1] + v*plan_ahead

    
    target = 731.8725791205976 #G1
    target = 902.1233635800432 #G2
    target = 857.1282220067823 #G3
    target = 679.2001084334775 #G4
    # nodes = [np.array([path_startX, (path_startX + interps[0][0][25])/2, interps[0][0][25]]).flatten(),
    #          np.array([path_startY, (path_startY + interps[0][1][25])/1.5, interps[0][1][25]]).flatten()]
    
    # nodes_local = [np.array([path_startX_local, (path_startX_local + interps_local[0][0][25])/2, interps_local[0][0][25]]).flatten(),
    #          np.array([path_startY, (path_startY + interps_local[0][1][25])/1.5, interps_local[0][1][25]]).flatten()]
    # print(nodes[0], nodes[0][1])
    # print(interps[0][1][25])
    turn_radius = v*1.94384**2/(11.26*math.tan(np.deg2rad(73)))
    turn_radius*=0.3048
    
    # print(interps_local)
    # bez = manual_bez([nodes[0][0],nodes[1][0]], [nodes[0][1], nodes[1][1]], [nodes[0][2], nodes[1][2]], 200)

    # bez_local = manual_bez([nodes_local[0][0],nodes_local[1][0]], [nodes_local[0][1], nodes_local[1][1]], [nodes_local[0][2], nodes_local[1][2]], 200)
    valid = False
    gi = [1,2]
    t_t = target-total_local-plan_ahead
    # new_travel, new_total = getTravelXY(interps, gi, toa, v)

    new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
    t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
    # print(new_total)
    t_t = target-new_total-plan_ahead
    
    print('LR',lr)
    entry = False
    entryTOA = plan_ahead
    while not valid:
        show = False
        # new_travel, new_total = getTravelTime(interps, gi, new_travel)
        new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
        t_t = target-new_total-entryTOA
        # nodes = [np.array([path_startX, (path_startX + interps[gi[0]][0][gi[1]])/2, interps[gi[0]][0][gi[1]]]).flatten(),
        #      np.array([path_startY, (path_startY + interps[gi[0]][1][gi[1]])/1.5, interps[gi[0]][1][gi[1]]]).flatten()]
        nodes = [np.array([path_startX_local, interps_local[gi[0]][0][gi[1]]+lr*4000*np.abs(np.cos(t_h)), interps_local[gi[0]][0][gi[1]]]).flatten(),
             np.array([path_startY_local, interps_local[gi[0]][1][gi[1]]+lr*4000*np.abs(np.sin(t_h)), interps_local[gi[0]][1][gi[1]]]).flatten()]
        # nodes = [np.array([path_startX, (interps[gi[0]][0][gi[1]])-10000, interps[gi[0]][0][gi[1]]]).flatten(),
        #      np.array([path_startY, (interps[gi[0]][1][gi[1]])+10000, interps[gi[0]][1][gi[1]]]).flatten()]
        m = (interps_local[gi[0]][1][0] - interps_local[gi[0]][1][gi[1]])/(interps_local[gi[0]][0][0] - interps_local[gi[0]][0][gi[1]])
        b = interps_local[gi[0]][1][0] - m*interps_local[gi[0]][0][0]
        coeffs = [m, b]
        # t_h = np.arctan2(interps[gi[0]][1][gi[1]]- interps[gi[0]][1][gi[1]-2],interps[gi[0]][0][gi[0]]- interps[gi[0]][0][gi[1]-2] )
        t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
        # print(np.rad2deg(t_h))
        optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                    target_toa=t_t,
                                    guess = [nodes[0][1], nodes[1][1]],
                                    target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs)
        print(f'BEZ CURVATURE {curv1}, TURN RADIUS {turn_radius}')
        if optim1[0] < ac_met[0]:
            lr = -1
        else:
            lr = 1
        optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
        optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)
        path_toa = optim_b_length/v
        # print(2*np.pi+np.arctan2(optim_b[1][20] - ac_met[1],  optim_b[0][20]-ac_met[0]), head)
        nodes1 = [np.array([nodes[0][0], optim1[0], nodes[0][2]]).flatten(), np.array([nodes[1][0], optim1[1], nodes[1][2]])]
        # if np.isclose(2*np.pi+np.arctan2(optim_b[1][20] - ac_met[1],  optim_b[0][20]-ac_met[0]), head, atol = np.deg2rad(10)):
        # print('HEAIDNG:',np.rad2deg(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])), np.rad2deg(np.pi/2), np.deg2rad(10))
        if np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(10)):
            entryTOA = plan_ahead
            path_toa = optim_b_length/v + entryTOA
            entry = False
            # print('using plan ahead TOA')
        else:
            entry = True
            ba, t_entry = solve_optimEntry(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
            
            x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [nodes1[0][2],nodes1[1][2]], t_entry)
            
            x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                        [optim1[0],optim1[1]],
                                        [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
            req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
            
            
            x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotate_head)
            act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[len(y_entry)-2], x_entry[-1]-x_entry[len(x_entry)-2]))
            

            partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                    [optim1[0],optim1[1]],
                                    [nodes1[0][2],nodes1[1][2]], 200, t_entry)
            pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],t_entry)
            pb_toa = pb_length/v
            path_toa = entryTOA + pb_toa

        ang_y = (nodes[0][2] - 5000) *m + b
        plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
        plt.plot(x_rot, y_rot)
        plt.scatter(ac_met[0], ac_met[1])
        plt.plot(x_rot, y_rot)
        plt.scatter(optim1[0], optim1[1])
        plt.plot(optim_b[0], optim_b[1], color = 'black')
        plt.axis('equal')
        plt.text(interps_local[-1][0][0], interps_local[-1][1][0], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
        plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

        if entry:
            plt.plot(x_entry, y_entry, color = 'cyan')
            plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
            plt.text(optim_b[0][50], optim_b[1][50]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
        #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
        else:
        #     print(f'Bez Travel Time: {path_toa:0.2f}s')
            plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        plt.show()
        
        
        
        
        print(new_total + path_toa)
        print(f'TARGET: {target}, ACTUAL TRAVEL TIME: {path_toa+new_total}')
        print(np.abs(target - (path_toa+new_total)))
        print(f'Original Target: {target}, New STAR Travel Time: {new_total}, New Target Time: {t_t}, Path TOA: {path_toa}')
        print(new_total+path_toa, target)
        if 1<= np.abs(target - (path_toa+new_total)) <= 3:
            g = 5000
            show = False
            for k in range(0, 5):
                # new_travel, new_total = getTravelTime(interps, gi, new_travel)
                new_travel, new_total = getTravelTime(interps_local, gi, toa, v)
                t_t = target-new_total-entryTOA
                # nodes = [np.array([path_startX, (path_startX + interps[gi[0]][0][gi[1]])/2, interps[gi[0]][0][gi[1]]]).flatten(),
                #      np.array([path_startY, (path_startY + interps[gi[0]][1][gi[1]])/1.5, interps[gi[0]][1][gi[1]]]).flatten()]
                nodes = [np.array([path_startX_local, interps_local[gi[0]][0][0] + lr*g*np.abs(np.cos(t_h)), interps_local[gi[0]][0][gi[1]]]).flatten(),
                    np.array([path_startY_local, interps_local[gi[0]][1][0]+lr*g*np.abs(np.sin(t_h)), interps_local[gi[0]][1][gi[1]]]).flatten()]
                # nodes = [np.array([path_startX, (interps[gi[0]][0][gi[1]])-10000, interps[gi[0]][0][gi[1]]]).flatten(),
                #      np.array([path_startY, (interps[gi[0]][1][gi[1]])+10000, interps[gi[0]][1][gi[1]]]).flatten()]
                m = (interps_local[gi[0]][1][0] - interps_local[gi[0]][1][gi[1]])/(interps_local[gi[0]][0][0] - interps_local[gi[0]][0][gi[1]])
                b = interps_local[gi[0]][1][0] - m*interps_local[gi[0]][0][0]
                coeffs = [m, b]
                # t_h = np.arctan2(interps[gi[0]][1][gi[1]]- interps[gi[0]][1][gi[1]-2],interps[gi[0]][0][gi[0]]- interps[gi[0]][0][gi[1]-2] )
                t_h = np.arctan2(interps_local[gi[0]][1][gi[1]]- interps_local[gi[0]][1][gi[1]-2],interps_local[gi[0]][0][gi[0]]- interps_local[gi[0]][0][gi[1]-2] )
                # print(np.rad2deg(t_h))
                optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                            target_toa=t_t,
                                            guess = [nodes[0][1], nodes[1][1]],
                                            target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = lr, line = coeffs)
                if optim1[0] < ac_met[0]:
                    lr = -1
                else:
                    lr = 1
                optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
                optim_b_length = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)
                path_toa = optim_b_length/v
                # print(2*np.pi+np.arctan2(optim_b[1][20] - ac_met[1],  optim_b[0][20]-ac_met[0]), head)
                nodes1 = [np.array([nodes[0][0], optim1[0], nodes[0][2]]).flatten(), np.array([nodes[1][0], optim1[1], nodes[1][2]])]
                # if np.isclose(2*np.pi+np.arctan2(optim_b[1][20] - ac_met[1],  optim_b[0][20]-ac_met[0]), head, atol = np.deg2rad(10)):
                # print('HEAIDNG:',np.rad2deg(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0])), np.rad2deg(np.pi/2), np.deg2rad(10))
                if target<=800 and np.isclose(np.arctan2(optim_b[1][1] - optim_b[1][0],  optim_b[0][1]-optim_b[0][0]), np.pi/2, atol = np.deg2rad(10)):
                    entryTOA = plan_ahead
                    path_toa = optim_b_length/v + entryTOA
                    entry = False
                    # print('using plan ahead TOA')
                else:
                    entry = True
                    ba, t_entry = solve_optimEntry(np.deg2rad(bank), np.deg2rad(73), np.deg2rad(min), nodes1, v*1.94384, ac_met, lr, rotate_head, show)
                    
                    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                                [optim1[0],optim1[1]],
                                                [nodes1[0][2],nodes1[1][2]], t_entry)
                    
                    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                                [optim1[0],optim1[1]],
                                                [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
                    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
                    
                    
                    x_entry, y_entry, ca, entryLength, entryTOA, h_c, k_c = entryPath(v, np.rad2deg(ba), [x_int_en, y_int_en], ac_met, lr, rotate_head)
                    act_ent = np.rad2deg(np.arctan2(y_entry[-1]-y_entry[198], x_entry[-1]-x_entry[198]))
                    

                    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                            [optim1[0],optim1[1]],
                                            [nodes1[0][2],nodes1[1][2]], 200, t_entry)
                    pb_length = optim_b_length - path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],t_entry)
                    pb_toa = pb_length/v
                    path_toa = entryTOA + pb_toa
                ang_y = (nodes[0][2] - 5000) *m + b
                plt.plot([nodes[0][2],nodes[0][2] - 5000], [nodes[1][2], ang_y], color = 'green', linestyle = '--')
                plt.plot(x_rot, y_rot)
                plt.scatter(ac_met[0], ac_met[1])
                plt.plot(x_rot, y_rot)
                plt.scatter(optim1[0], optim1[1])
                plt.plot(optim_b[0], optim_b[1], color = 'black')
                plt.axis('equal')
                plt.text(interps_local[-1][0][0], interps_local[-1][1][0], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
                plt.text(interps_local[1][0][0], interps_local[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

                if entry:
                    plt.plot(x_entry, y_entry, color = 'cyan')
                    plt.text(x_entry[-1]-5000, y_entry[-1]-5000, f'Entry Path Travel Time: \n{entryTOA:0.2f}', weight = 'bold')
                    plt.text(optim_b[0][50], optim_b[1][50]-5000, f'Partial Bez Travel Time: {pb_toa}s', weight = 'bold' )
                #     print(f'Partial Bez Travel Time: {pb_toa}, Entry Path Travel Time: {entryTOA}')
                else:
                #     print(f'Bez Travel Time: {path_toa:0.2f}s')
                    plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
                plt.show()
        
                g+=2500
                if np.abs(target - (path_toa+new_total))<=1:
                    print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
                    print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
                    if entry == True and (not np.isclose(x_entry[-1],partial_bez1[0][0], atol = 1) or not np.isclose(y_entry[-1],partial_bez1[1][0], atol = 1)):
                        valid = False
                        if bank > min:
                            bank -=2.5
                        elif bank <= min:
                            bank = min
                    else:
                        valid = True
                        break

        if np.abs(target - (path_toa+new_total))<=1:
            print(f'X ENTRY, PB1 Start: {x_entry[-1]}, {partial_bez1[0][0]}')
            print(f'Y ENTRY, PB1 Start: {y_entry[-1]}, {partial_bez1[1][0]}')
            if entry == True and (not np.isclose(x_entry[-1],partial_bez1[0][0], atol = 5) or not np.isclose(y_entry[-1],partial_bez1[1][0], atol = 5)):
                valid = False
                bank-=2.5
            else:
                valid = True
        

        if valid == False:
            # print(f'Target TOA is: {target-total-plan_ahead}\nPath TOA is: {path_toa}\n Validity is {valid}')
            gi[1]+=2
            # if gi[1] >= interp_pts/2 or gi[0]>=2:
            #     gi[1]+=3
            if gi[1]>=interp_pts:
                gi[0]+=1
                gi[1] = 0

            if gi[0] >= len(interps):
                gi = [0, 0]
                nodes = [np.array([ac_met[0], (ac_met[0] + interps[gi[0]][0][gi[1]])/2, interps[gi[0]][0][gi[1]]]).flatten(),
                np.array([ac_met[1], (ac_met[1] + interps[gi[0]][1][gi[1]])/2+10, interps[gi[0]][1][gi[1]]]).flatten()]
                m = (interps[gi[0]][1][10] - interps[gi[0]][1][0])/(interps[gi[0]][0][10] - interps[gi[0]][0][0])
                b = interps[gi[0]][1][0] - m*interps[gi[0]][0][0]
                coeffs = [m, b]
                t_h = np.arctan2(interps[gi[0]][1][10]- interps[gi[0]][1][0],interps[gi[0]][0][10]- interps[gi[0]][0][0] )
                optim1, curv1 = solve_optim1(P0 = [nodes[0][0],nodes[1][0]], P2 = [nodes[0][2], nodes[1][2]],
                                            target_toa=target-total,
                                            guess = [nodes[0][1], nodes[1][1]],
                                            target_heading=t_h, velocity=v, turn_radius=turn_radius, lr = 1, line = coeffs)
                t_t= target-total
                optim_b = manual_bez([nodes[0][0],nodes[1][0]], [optim1[0], optim1[1]], [nodes[0][2], nodes[1][2]], 200)
                path_toa = path_length([optim1[0], optim1[1]],[nodes[0][0],nodes[1][0]],[nodes[0][2], nodes[1][2]],1)/v
                break
            # print(gi)

    stop_time = time.time()
    print(f'TOTAL OPTIMIZATION TIME: {stop_time-start_time}')
    
    print("ENTRY TOA:", entryTOA)
    print("Bez TOA:", pb_toa)
    print(f'Target Total is {target}\nTarget Bez TOA is: {t_t}\nBez+Entry Path TOA is: {path_toa}\nTotal Travel Time w/ STAR Route is {path_toa+new_total}\nValidity is {valid}')
    # print(curv1, turn_radius/.3048)
    print(f'STAR Path Intercept is {gi}')

    

    rotated_bezier = rotate_bez(np.array(optim_b), -angle_to_rotate, origin)
    rotated_partial = rotate_bez(np.array(partial_bez1), -angle_to_rotate, origin)
    rotated_entry = rotate_bez(np.array([x_entry, y_entry]), -angle_to_rotate, origin)
    rot_nodes = rotate_bez(nodes1, -angle_to_rotate, origin)
    m2 = (interps[gi[0]][1][0] - interps[gi[0]][1][gi[1]])/(interps[gi[0]][0][0] - interps[gi[0]][0][gi[1]])
    b2 = interps[gi[0]][1][0] - m2*interps[gi[0]][0][0]
    # plt.plot(rotated_entry[0], rotated_entry[1])
    # plt.plot(rotated_bezier[0], rotated_bezier[1])
    show_local = False

    if entry:
        print(f'Entry Bank{np.rad2deg(ba)}, Entry Point On Bez {t_entry}')
        print('REQUIRED HEADING AT ENTRY:', req_ent-angle_to_rotate-360)
        print('ACTUAL HEADING AT INTERSECT:', act_ent-angle_to_rotate-360)
            
        print('DIFF AT ENTRY:', req_ent-act_ent)
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
        plt.plot(rotated_entry[0], rotated_entry[1], label = 'Entry Arc', color = 'cyan')
        plt.scatter(rotated_entry[0][-1], rotated_entry[1][-1], marker = '*', color = 'yellow', label = 'Intersection Point', s = 100, zorder = 100)
        plt.text(x_entry[-1]+10000, y_entry[-1]-5000, f'Entry Path Travel Time: {entryTOA:0.2f}', weight = 'bold')
        plt.plot(rotated_partial[0], rotated_partial[1], c = 'magenta', label = 'Partial Bez Path')
        plt.text(rotated_bezier[0][50], rotated_bezier[1][50]-2500, f'Partial Bez Travel Time: {pb_toa:0.2f}s', weight = 'bold' )
    else:
        print(f'No Entry Path Used\nRequired Bez Heading: {np.rad2deg(np.arctan2(rotated_bezier[1][0] - rotated_bezier[1][1], rotated_bezier[0][0] - rotated_bezier[0][1]))+180}')
        print(f'Actual Heading: {np.rad2deg(head)}')
        # plt.text(optim_b[0][50], optim_b[1][50], f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )
        plt.text(rotated_bezier[0][50], rotated_bezier[1][50]-2500, f'Bez Travel Time: {path_toa:0.2f}s', weight = 'bold' )

    plt.arrow(ac_met[0], ac_met[1], 2500*np.cos(head), 2500*np.sin(head), width=50,
                length_includes_head=False,
                head_width=1000,
                head_length=1000,
                head_starts_at_zero=False,
                facecolor='black',
                edgecolor='black',
                label = f'Heading of {np.rad2deg(head)-180}')
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
        plt.scatter(interps[i][0], interps[i][1], color = 'orange', marker = 's')
    ang_y = (rot_nodes[0][2] - 5000) *m2 + b2
    plt.plot([rot_nodes[0][2],rot_nodes[0][2] - 5000], [rot_nodes[1][2], ang_y], color = 'green', linestyle = '--')
    plt.plot(rotated_bezier[0], rotated_bezier[1], color = 'black', linestyle = '--')
    plt.scatter(rot_nodes[0][1], rot_nodes[1][1], marker = '^', color = 'red')
    plt.scatter(rot_nodes[0][2], rot_nodes[1][2], marker = '^', color = 'red', label = 'Control Points')
    plt.scatter(path_startX, path_startY, marker = '^', color = 'red')
    plt.scatter(x, y, label = 'STAR Path')
    plt.text(interps[0][0][0]-10000, interps[0][1][0], f'Required Travel Time: {target:0.2f}\nActual Travel Time: {path_toa+new_total:0.2f}', weight = 'bold')
    plt.text(interps[1][0][0], interps[1][1][0], f'STAR Travel Time: {new_total:0.2f}s', weight = 'bold' )

    plt.scatter(ac_met[0], ac_met[1], color = 'green')
    plt.scatter(0, 0, marker = '*', color = 'green', label = 'CMH')
    plt.ylabel('Y (m)')
    plt.xlabel('X (m)')
    plt.grid()
    plt.legend(loc = 'lower left')
    plt.axis('equal')
    plt.show()