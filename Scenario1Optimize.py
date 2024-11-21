import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy
import sympy as sp
from scipy.optimize import brentq
from scipy.optimize import fsolve
from pyproj import Proj
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
            # {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]}
            {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
        ) 

    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

    return val.x , curvature(P0,val.x,P2)


def solve_optimEntry(guess, max_bank,  min_bank, nodes, velocity):
    def path_cost(ba, nodes, velocity):
        diff = find_diff_entry(ba, nodes, velocity)

        return np.abs(diff[0])
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}  
    )
    val= minimize(path_cost,guess,(nodes, velocity), method='SLSQP', tol=1E-10, constraints=cons)
    t_val = find_diff_entry(val.x, nodes, velocity)[1]
    return val.x, t_val

def solve_optimExit(guess, max_bank, min_bank, min_t, nodes, velocity):
    #t is also being minimized, so part of constraints
    def path_cost(guess, nodes, velocity):
        diff = find_diff_exit(guess, nodes, velocity)
        # print('DIFF FOR EXIT:', diff[0])
        return np.abs(diff[0]) 

    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank},
            {'type': 'ineq', 'fun': lambda x: x[1] - min_t},
            {'type': 'ineq', 'fun': lambda x: 1-x[1]},
    )
    val =  minimize(path_cost, guess, (nodes, velocity), method = 'SLSQP', tol = 1E-10, constraints=cons)
    return val.x

def find_diff_entry(ba, nodes, velocity):
    pi = np.pi
    t_guess = 0.75
    mindiff = 50
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    tr = velocity**2 / (11.26*math.tan(ba[0]))
    tr*=0.3048
    h = tr + 229
    k = -391
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

            x_l = [i for i in np.linspace(229, Bx(i), 200)] #Space between nominal path and bez
            y = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            # plt.plot(path[0], path[1])
            # plt.plot(x_l, y)
            # plt.axis('equal')
            # plt.show()

            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            diff = np.abs(bez_angle-int_angle)
            if diff < mindiff and np.abs(y[-1] - By(i))<=.00001 and np.abs(x_l[-1]-Bx(i))<=0.00001:
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

    return [mindiff, t_final]

def find_diff_exit(guess, nodes, velocity):
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
    h = lambda t: Bx(t) - tr*math.cos(bezHead)
    k = lambda t: By(t) + tr*math.sin(bezHead)

    circle_eq = lambda y: (229-h(t_guess))**2+(y-k(t_guess))**2 - tr**2
    y_guess = 610
    S = fsolve(circle_eq, y_guess) #gives y intersection
    t_final = .5

    y_l = [i for i in np.linspace(274, S, 200)]
    x = [229 for i in y_l]

    # print(S)
    for i in S:
        if i > 0: 
            y = [i for i in np.linspace(S, By(t_guess))]
            x_l = [h(t_guess) - np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            if x_l[0] ==229:# <= 0.01:
                # x_l = [i for i in np.linspace(750, Bx(t_guess))]
                # y = [k(t_guess)-np.sqrt(tr**2 - (x-h(t_guess))**2) for x in x_l]
                int_angle = np.arctan2(y[0]-y[1], x_l[0]-x_l[1])
                diff = np.abs((np.pi/2) - int_angle)
                guess_deg = np.rad2deg(ba)
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
    points = np.array(points)
    if lr == 1:
        for i in range(len(points)):
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
        
def toCallOutside(velocity, turn_rate, target_toa1, target_toa2, uav_head, nodes1, nodes2, koz_x, koz_bot, koz_top, lr):
    # velocity  = 188 #ft/s
    turn_rate = np.deg2rad(4.50) # RAD/s
    turn_radius = velocity / turn_rate
    valid1 = False
    valid2 = False
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    print(nodes2[1][1])
    m_p21 = (nodes2[1][1] - nodes1[1][2])/(nodes2[0][1] - nodes1[0][2])
    print(m_p21)
    b2 = nodes2[1][1] - m_p21*nodes2[0][1]
    lines_coeffs2 = [m_p21, b2]
    if lr == 1:
        # corridor = nodes1[0][0]+0.002055
        corridor = Meters_To_WSG84([0,750],[nodes1[0][0], nodes1[1][0]])[1]
    else:
        corridor = -750
        # corridor = nodes1[0][0]-0.002055
    c = 0
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
        conv14 = Meters_To_WSG84([0,1400],[nodes1[1][0], nodes1[0][0]])[1]
        conv15 = Meters_To_WSG84([0,1550],[nodes1[1][0], nodes1[0][0]])[1]
        print(conv14, conv15)
        x_slope = [i for i in np.linspace(conv14, conv15, 50)]
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
            x_slope = [nodes1[0][2] for i in np.linspace(conv14, conv15, 50)]
            y_slope = [i for i in np.linspace(optim_sol1[1]-50, optim_sol2[1]+50, 50)]

        
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
        valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 1', lr)
        # print("OPTIM SOL: ", optim_sol1[0], optim_sol1[1])
        # print(target_toa1, target_toa2)
        plt.plot(optimal_bez1[0], optimal_bez1[1])
        plt.plot(optimal_bez2[0], optimal_bez2[1])
        plt.axis('equal')
        plt.show()
        # print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.1
            # nodes1[1][1]+=-25
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            # nodes2[1][1]+=25
            c+=1
        if c>30:
       
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
            ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
            ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
            # ax.plot(optimal_bezPre[0], optimal_bezPre[1], c = 'black')
            valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
            valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
            # print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
            ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
            ax.scatter(optim_sol2[0], optim_sol2[1], marker = '*', s = 1000, zorder = 1000)
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
            bx = [500 for i in range(koz_top)]
            bx2 = [1000 for i in range(koz_top)]
            ybot = [koz_bot for i in range(500)]
            bxbot = [i for i in range(500, 1000)]
            ytop = [koz_top for i in range(500)]
            y2 = [i for i in range(-50, 950)]
            ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
            ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')

            ax.grid(True)
            ax.axis('equal')

            ax.legend(loc = 'center left', fontsize = '8')
            plt.show()

    vel_knots = 111.6
    nodes1 = [np.array([nodes1[0][0], optim_sol1[0], nodes1[0][2]]).flatten(),np.array([nodes1[1][0], optim_sol1[1], nodes1[1][2]]).flatten()]
    ba, t_entry = solve_optimEntry(np.deg2rad(25), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots)
    
    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    x_entry, y_entry, central_angle, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en])

    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200, t_entry)


    t_start2 = 0.35678391959798994

    nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]
    baEx, exit_t = solve_optimExit([np.deg2rad(27.5), t_start2], np.deg2rad(30), np.deg2rad(20), t_start2, nodes2, velocity)
    print('EXIT BANK', np.rad2deg(baEx) )

    x_int_ex, y_int_ex = find_bez_xy([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], exit_t)
    
    x_exit, y_exit, central_angle_ex, exitLength, exitTOA, h_ex, k_ex = exitPath(velocity=velocity, t_exit=exit_t,
                                                                                 ba = baEx, intersect=[x_int_ex, y_int_ex],
                                                                                 nodes = [[nodes2[0][0],nodes2[1][0]],
                                                                                            [optim_sol2[0],optim_sol2[1]],
                                                                                            [nodes2[0][2],nodes2[1][2]]])
    x_exit[0] = 750
    y_exit[0] = 2350

    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200, exit_t)

    # pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_entry)
    # pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]], P2=[nodes2[0][2],nodes2[1][2]], t=t_start2)


    paths = [[x_entry, y_entry], partial_bez1, partial_bez2, [x_exit, y_exit]]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 5)





    plt.figure()
    plt.plot(partial_bez1[0], partial_bez1[1], color = 'magenta')
    plt.plot(partial_bez2[0], partial_bez2[1], color = 'magenta')
    plt.plot(x_entry, y_entry, color = 'cyan')
    plt.plot(x_exit, y_exit, color = 'orange')
    plt.show()


    return wpts, x_wpts, y_wpts, optimal_bez1, optimal_bez2 

def entryPath(velocity, ba, intersect):
    if intersect[0] == 229:
        intersect[0] = 400
        ba = 25
    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
    tr*=0.3048
    print('TR', tr)
    h, k  = tr+229, -391
    print('INTERSECT', intersect)
    
    x_entry = [i for i in np.linspace(229, intersect[0], 200)]
    
    y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]

    # y_entry[-1] = intersect[1]
    central_angle = pi - np.arctan2(y_entry[-1]-k, x_entry[-1] - h)
    entryLength = tr*central_angle #entryLength = 2*pi*tr * (central_angle/(2*pi))
    
    entryTOA = entryLength/velocity
    print('ENTRY ToA: ', entryTOA)

    return x_entry, y_entry, central_angle, entryLength, entryTOA, h, k


def exitPath2(velocity, t_start, path):
    pi = np.pi

    t_vals = []
    central_angle_l = []
    exitLength_l = []
    exitTOA_l = []
    # sols = []
    angles = []
    exit_angle = []
    int_point = []
    center = []
    bas = []
    # y_path = [i for i in range(900, 3000)]
    # x = [750 for i in y_path]

    for t in np.linspace(t_start, .9, 75):
        for ba in np.linspace(25, 45, 100):
            index = int(np.round(t*200))
            bezPos = [path[0][index], path[1][index]]
            bezPosPrev = [path[0][index-1], path[1][index-1]]
            bezHead = np.arctan2(bezPos[1] - bezPosPrev[1], bezPos[0] - bezPosPrev[0])

            tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
            h = bezPos[0] - tr*math.cos(bezHead)
            k = bezPos[1] + tr*math.sin(bezHead)

            circ_int = np.sqrt(tr**2 - (750-h)**2)+k #Will be used to check if it intercepts
            x_l = [i for i in np.linspace(750, bezPos[0])] #Space between nominal path and bez
            y = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            int_angle = np.arctan2(y[0]-y[1], x_l[0] - x_l[1])

            if circ_int > 0 and np.isclose(pi/2, int_angle, np.deg2rad(15)):

                central_angle = pi + np.arctan2(y[-1]-k, x_l[-1] - h)
                central_angle_l.append(central_angle)

                exitLength = 2*pi*tr * (central_angle/(2*pi))
                exitLength_l.append(exitLength)

                exitTOA = exitLength/velocity
                exitTOA_l.append(exitTOA)

                int_point.append(circ_int)
                t_vals.append(t)
                angles.append(ba)
                center.append([h, k])
                exit_angle.append(int_angle)
                bas.append(ba)
    
    
    mindex = exitTOA_l.index(min(exitTOA_l))
    tFinal = int(np.round(t_vals[mindex]*200))
    print('FINAL BANK ANGLE:', bas[mindex])
    realT = t_vals[mindex]
    centerFinal = center[mindex]

    bezPosFinal = [path[0][tFinal], path[1][tFinal]]
    trFinal = 111.6**2/(11.26*np.tan(np.deg2rad(angles[mindex])))

    x_l = [i for i in np.linspace(750, bezPosFinal[0], 200)] #Space between nominal path and bez
    y = [centerFinal[1]-np.sqrt(trFinal**2 - (x-centerFinal[0])**2) for x in x_l] #arc created by turn

    print('T VALUE:', realT)
    print('TOA:', exitTOA_l[mindex])
    print('BANK ANGLE:', angles[mindex])
    print('INTERCEPTION ANGLE:', np.rad2deg(exit_angle[mindex]))
    print('INTERCEPTION POINT:', (750,int_point[mindex]))

    return x_l, y, bezPosFinal, realT, exitTOA_l[mindex], angles[mindex], exitLength_l[mindex]


def exitPath(velocity, t_exit, ba, intersect, nodes):
    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(ba))
    tr*=0.3048

    Bx = lambda t: nodes[1][0] + (nodes[0][0] - nodes[1][0]) * (1 - t)**2 + (nodes[2][0] - nodes[1][0]) * t**2
    By = lambda t: nodes[1][1] + (nodes[0][1] - nodes[1][1]) * (1 - t)**2 + (nodes[2][1] - nodes[1][1]) * t**2

    bezHead = np.arctan2(By(t_exit)-By(t_exit-0.01), Bx(t_exit)-Bx(t_exit-0.01))

    h, k = intersect[0] - tr*math.cos(bezHead), intersect[1]+tr*math.sin(bezHead)

    x_exit = [i for i in np.linspace(229, intersect[0], 200)]
    y_exit = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_exit]
    # y_exit[0] = 650
    # y_exit[0] = k
    central_angle = pi + np.arctan2(y_exit[-1]-k, x_exit[-1] - h)

    exitLength = 2*pi*tr * (central_angle/(2*pi))

    exitTOA = exitLength/velocity
    
    print('EXIT ToA:', exitTOA)

    return x_exit, y_exit, central_angle, exitLength, exitTOA, h, k





if __name__ == "__main__":
 
    velocity  = 57.412 #m/s
    knots = velocity*1.94384
    # turn_rate = np.deg2rad(20) # RAD/s
    # turn_radius = velocity/turn_rate
    # print(turn_radius)
    turn_radius = knots**2/(11.26*math.tan(np.deg2rad(60)))
    turn_radius*=0.3048
    # print(np.rad2deg(velocity/turn_radius))
    # print(np.rad2deg(np.arctan(57.412**2/(11.26*turn_radius))))

    h = 225
    koz_bot = 15
    koz_top = h-15

    
    ev_toa = 275/64.008
    
    c=0
    target_toa1 = 2*ev_toa
    target_toa2 = 2*ev_toa
    uav_head = np.deg2rad(90)
    lr = 1
    # if lr == -1:
    #     corridor = 0
    #     koz_x = 500
    #     nodes1 = [np.array([750, 1400, 25]).flatten(),np.array([-50, 10, 450]).flatten()]
    #     nodes2 = [np.array([25, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()]
    # else:
    corridor = 457
    koz_x = 305
    nodes1 = [np.array([229, 442, 450]).flatten(),np.array([0, 5/60, h/2]).flatten()]
    nodes2 = [np.array([450, 454, 229]).flatten(),np.array([h/2, h-100, h]).flatten()]
    
 
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
        
        # m_p12 = (nodes1[1][2]-optim_sol1[1])/(nodes1[0][2]-optim_sol1[0])
        # # print(m_p12)
        # mark = True
        # #y = mx+b
        # if np.isclose(nodes1[0][2], optim_sol1[0], atol = 0.01) == False:
        #     mark = False

        # b = optim_sol1[1] - m_p12*optim_sol1[0]
        # lines_coeffs = [m_p12, b]
        # x_slope = [i for i in np.linspace(427, 472, 50)]
        # y_slope = [m_p12*x + b for x in x_slope]
        m_p12 = (nodes1[1][2]-optim_sol1[1])/(nodes1[0][2]-optim_sol1[0])


        b = optim_sol1[1] - m_p12*optim_sol1[0]
        lines_coeffs = [m_p12, b]
        # print('M12',m_p12, 'B', b)

        x_slope = [i for i in np.linspace(optim_sol1[0]-10, 470, 50)]
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
        # print('CHECK CHECK',m_p21*optim_sol1[0]+b2, optim_sol1[1])
        # if not np.isclose(m_p21*optim_sol1[0]+b2, optim_sol1[1], atol=0.01):
        #     valid1 = False

        # if not np.isclose(m_p12*optim_sol2[0]+b, optim_sol2[1], atol = 0.01):
        #     valid2 = False
        # print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.1
            nodes1[1][1]+=-1
            nodes1[0][1]+=1
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            nodes2[0][1]+=-1
            nodes2[1][1]+=1
            c+=1
        # if c>30:
        # if valid1 or valid2:
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
        ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
        ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
        ax.plot(x_slope, y_slope, color = 'purple', linestyle = 'dashed', label = 'G1 Continuity Line')
        # # ax.plot(optimal_bezPre[0], optimal_bezPre[1], c = 'black')
        # valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        # valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
        # if not np.isclose(m_p21*optim_sol1[0]+b2, optim_sol1[1], atol=0.01):
        #     valid1 = False

        # if not np.isclose(m_p12*optim_sol2[0]+b, optim_sol2[1], atol = 0.01):
        #     valid2 = False
        # print('CHECK CHECK',m_p12*optim_sol2[0]+b, optim_sol2[1])
        # # if m_p12*optim_sol2[0]+b != optim_sol2[1]:
        # #     valid2 = False
        # # if m_p21*optim_sol1[0]+b2 != optim_sol1[1]:
        # #     valid1 = False
        # print(valid1, valid2)
        # print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
        ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
        ax.scatter(optim_sol2[0], optim_sol2[1])
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
    ba, t_entry = solve_optimEntry(np.deg2rad(28), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots)
    print(np.rad2deg(ba), t_entry)


    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    
    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry+0.01)
    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))
    print('REQUIRED HEADING AT ENTRY:', req_ent)
    
    
    x_entry, y_entry, central_angle, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en])
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
    # paths = [partial_bez1, partial_bez2]
    # print(len(paths))
    # print(paths[0][0])
    # wpts, x_wpts, y_wpts = paths_to_wp(paths, 10)
    # print(wpts)

    print('CURVHSADHAKSD',curvature(P0 = [nodes1[0][0],nodes1[1][0]], P1 = [optim_sol1[0], optim_sol1[1]], P2 = [nodes1[0][2], nodes1[1][2]]), turn_radius)  
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



    # ax.plot(x_exit, y_exit, color = 'orange', label = 'Exit Path')
    # ax.scatter(h_ex, k_ex, marker = 's', color = 'green', label = 'Exit Circle Center')
    # ax.scatter(x_exit[-1], y_exit[-1], color = 'red', marker = '*', s = 100, zorder = 100, label = 'Exit Point')

    # # print('EXIT COORDS:', x_exit, y_exit)
    # ax.scatter(x_exit[0], y_exit[0], color = 'purple', marker = '*', s = 100, zorder = 100, label = 'Interception Point')

    ax.plot([122, 305, 488], [h/2, h/2, h/2], linestyle = 'dashdot', alpha = 0.5)

    # ax.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial QBC')
    # ax.plot(partial_bez2[0], partial_bez2[1], c = 'magenta')#, label = 'Partial QBC')

    ax.plot(optimal_bez1[0],optimal_bez1[1], c='black', label='Quadratic Bezier curve')


    # print('ENTRY ANGLE:', np.rad2deg(np.arctan2(partial_bez1[1][1] - y_entry[-1], partial_bez1[0][1] - x_entry[-1])))
    ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
    
    ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')

    # # ax.scatter(h_c, k_c, color = 'green', marker = 's', label = 'Entry Circle Center')
    # ax.scatter(x_entry[0], y_entry[0], marker = '^', color = 'black', label = 'Fleet Aircraft Start Position')


    # ax.plot(x_entry, y_entry, label = 'Entry Arc', color = 'cyan')
    # ax.scatter(x_entry[-1], y_entry[-1], marker = '*', color = 'purple', label = 'Intersection Point', s = 100, zorder = 100)
    # # print(y_circ[0])
    


    # print("ENTRY TOA:", entryTOA)
    # print("EXIT TOA:", exitTOA)
    # print("SOLVED BEZIER LEGNTH: ", pb1_length + pb2_length)
    # print("PARTIAL BEZIER TRAVEL TIME:", (pb1_length+pb2_length)/velocity)

    # print("SOLVED LENGTH WITH ENTRY/EXIT: ", pb1_length + optim2_length + entryLength+exitLength)
    # print("TOTAL TOA:", entryTOA + exitTOA + (pb1_length+pb2_length)/velocity)


    ax.plot(x_slope, y_slope, color = 'purple', linestyle = 'dashed', label = 'G1 Continuity Line')
    
    ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
    ax.text(nodes1[0][0]+1,nodes1[1][0]-35,  r'$\bf{p_{01}}$', fontsize=16)
    ax.text(optim_sol1[0]+10,optim_sol1[1],  r'$\bf{p_{11}}$', fontsize=16)
    ax.text(nodes1[0][2]+10,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$', fontsize=16)
    # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
    ax.text(optim_sol2[0]+10,optim_sol2[1],  r'$\bf{p_{12}}$', fontsize=16)
    ax.text(nodes2[0][2]+25,nodes2[1][2]+25,  r'$\bf{p_{22}}$', fontsize=16)
    ax.text(nodes1[0][2]-137,nodes1[1][2]-35,  r'Curve 1', fontsize = 18)
    ax.text(nodes1[0][2]-137,nodes1[1][2]+20,  r'Curve 2', fontsize = 18)

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

    ax.plot([229, 229], [-396, 600+50], color = 'green', label = 'Nominal Path')

    # ax.scatter(mx,my)
    # ax.scatter(center1x,center1y)
    # ax.scatter(center2x,center2y)
    # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
    # ax.add_patch(Circle((center2x, center2y), r2, color='black', fill=False))
 
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.legend(loc = 'upper left', fontsize = '10')
    # plt.ylim(optim_sol1[1]-150, optim_sol2[1]+150)
    # plt.axes.set_ylim(-150, 600)
    plt.show()
    bezXY = {
        'Bez1X': optimal_bez1[0],
        'Bez1Y': optimal_bez1[1],
        'Bez2X': optimal_bez2[0],
        'Bez2Y': optimal_bez2[1]
    }

    bez = pd.DataFrame(bezXY)
    bez.to_json('Bez3.json', orient = 'records', lines = 'True')
    print('CPX1',nodes1[0])
    print('CPY1', nodes1[1],)
    print('OS1', optim_sol1,)
    print('CPX2', nodes2[0],)
    print('CPY2', nodes2[1],)
    print('OS2', optim_sol2,)
    print('contX', x_slope,)
    print('contY', y_slope)

    print(f'\n\n\n TOA OF B1: {path_length(optim_sol1, [nodes1[0][0], nodes1[1][0]], [nodes1[0][2],nodes1[1][2]], 1)/velocity} TARGET TOA: {target_toa1}, RATIO: {(path_length(optim_sol1, [nodes1[0][0], nodes1[1][0]], [nodes1[0][2],nodes1[1][2]], 1)/velocity)/target_toa1}')
    print(f'\n\n\n TOA OF B2: {path_length(optim_sol2, [nodes2[0][0], nodes2[1][0]], [nodes2[0][2],nodes2[1][2]], 1)/velocity} TARGET TOA: {target_toa2}, RATIO: {(path_length(optim_sol2, [nodes2[0][0], nodes2[1][0]], [nodes2[0][2],nodes2[1][2]], 1)/velocity/target_toa2)}')

#line types for, dashed, solid, dots
#