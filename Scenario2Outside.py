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
            # {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
            # {'type': 'ineq', 'fun': lambda x: 1500 - x[0]},
            # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
            # {'type': 'ineq', 'fun': lambda x: x[0] - 1000},
            # {'type': 'ineq', 'fun': lambda x: x[1] -P0[1]},
            # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
            {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))},
            # {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
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
            # {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
            # {'type': 'ineq', 'fun': lambda x: 1500 - x[0]},
            # {'type': 'ineq', 'fun': lambda x: np.deg2rad(10) - np.abs(np.arctan2(x[1]-P2[1], x[0] - P2[0]))},
            # {'type': 'ineq', 'fun': lambda x: x[0] - 1000},
            # {'type': 'ineq', 'fun': lambda x: x[1] - P0[1]},
            # {'type': 'ineq', 'fun': lambda x: P2[1]-x[1]},
            {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P0[1]-x[1]), (P0[0]-x[0]))) - target_heading}, 
            {'type': 'ineq', 'fun': lambda x: x[1] - P0[1]},
            {'type': 'eq', 'fun': lambda x: x[1] - (line[0]*x[0] + line[1])}
        ) 
    
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

    return val.x , curvature(P0,val.x,P2)


def solve_optimEntry(guess, max_bank,  min_bank, nodes, velocity, pos, j):
    def path_cost(ba, nodes, velocity):
        diff = find_diff_entry(ba, nodes, velocity, pos, j)

        return np.abs(diff[0])
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}  
    )
    val= minimize(path_cost,guess,(nodes, velocity), method='SLSQP', tol=1E-10, constraints=cons)
    t_val = find_diff_entry(val.x, nodes, velocity, pos, j)[1]
    return val.x, t_val



def find_diff_entry(ba, nodes, velocity, pos,j ):
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
    h = tr
    k = pos[1]
    circle_eq = lambda t: (Bx(t)-h)**2+(By(t)-k)**2 - tr**2

    S = fsolve(circle_eq, t_guess)
    # print(S, np.rad2deg(ba[0]))
    t_final = 0

    for i in S:

        if i>=0 and i <= 1:
            index = int(np.round(float(i*200)))
            if index == 200:
                index = 199
            bez_angle = np.arctan2(By(i+.01)-By(i), Bx(i+0.01)-Bx(i))

            x_l = [i for i in np.linspace(0, Bx(i), 200)] #Space between nominal path and bez
            y = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            # if j == 2 or j == 4:
            #     plt.plot(path[0], path[1])
            #     plt.plot(x_l, y)
            #     plt.axis('equal')
            #     plt.show()

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
        path_wpts.append([path[0][-1], path[1][-1]])
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
    for i in range(len(points)):
        if points[i][0]<=x and points[i][1]>=535 and points[i][1]<=625 or points[i][0] >= corridor:
            # print('In Positive')
            # print(curve)
            # print(points[i][0], x, points[i][0]<=x)
            # print(points[i][1], y1, points[i][1]>=y1)
            # print(points[i][1], y2, points[i][1]<=y2)
            # print(points[i][0], corridor, points[i][0] >= corridor)
            # plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
            return False
    return True


def toCallOutside(velocity, turn_rate, target_toa1, target_toa2, uav_head, nodes1, nodes2, koz_x, koz_bot, koz_top, lr, latlon, corridor, ev_TOA):
    start_time = time.time()
    # turn_rate = np.deg2rad(4.50) # RAD/s
    # turn_radius = velocity / turn_rate
    turn_radius = 111.6**2/(11.26*math.tan(np.deg2rad(30)))
    turn_radius*=0.3048
    valid1 = False
    valid2 = False

    c = 0
    m_p21 = (nodes2[1][1]-nodes1[1][2])/(nodes2[0][1]-nodes1[0][2])
    b2 = nodes2[1][1] - m_p21*nodes2[0][1]

    lines_coeffs2 = [m_p21, b2]
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
        
        optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
        
        solved_heading = np.arctan2((optim_sol1[1]-nodes2[1][0]),(optim_sol1[1]-nodes2[0][0]))
        optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200)

        wpts, x_wpts, y_wpts = bez_to_wp(optimal_bez1, optimal_bez2, 15) 
        wpts_all1, wpts_all1x, wpts_all1y = bez_to_wp_single(optimal_bez1, 30) 
        # print(wpts_all1)
        wpts_all2, wpts_all2x, wpts_all2y = bez_to_wp_single(optimal_bez2, 30) 
        valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)

        optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
        optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]],P2=[nodes2[0][2], nodes2[1][2]], t=0.75)
        optim_toa = ((optim1_length+optim2_length) / velocity)
        # print('OPTIM TOA:', optim_toa, 'EV TOA', ev_TOA)

        if valid2 == False and optim_toa >= ev_TOA:
            valid2 = True
        if optim_toa < ev_TOA:
            valid2 = False
        if valid1 == False:
            target_toa1+=0.1
            nodes1[1][1]+=-.5
            nodes1[0][1]+=.5
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            nodes2[1][1]+=.5
            nodes2[0][1]+=-.5
            c+=1
        # if c>5:
            
        # print(target_toa2)
        # #     c+=1
        # # if c>=5:
            # y = [i for i in range(10)]
            # bx = [koz_x for i in range(10)]
            # ybot = [koz_bot for i in range(250)]
            # bxbot = [i for i in range(39)]
            # ytop = [koz_top for i in range(2)]
            # y2 = [i for i in range(-50, 950)]
        #     # xwall = [1500 for i in range(-50, 950)]
        # fig = plt.figure(1)
        # ax = fig.add_subplot(111)
        # ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
        # ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
        # ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')

        # print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
        # ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
        # ax.plot([corridor, corridor, corridor], [nodes1[1][0], nodes1[1][2], nodes2[1][2]], linestyle = '--')
        # ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
        # ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
        # ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
        # # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
        # ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
        # ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
        # ax.text(nodes1[0][2]-250,nodes1[1][2]-100,  r'Curve 1')
        # ax.text(nodes1[0][2]-250,nodes1[1][2]+100,  r'Curve 2')

        # # ax.text(mx+0.25, my, r'$\bf{m}$')
        # # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
        # # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
        # # ax.plot(bx, y, color = 'red', linestyle = '--')
        # # ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
        # # ax.plot([-83.2, -83.199], ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
        # plt.plot([koz_x, koz_x], [koz_bot, koz_top], color = 'red', linestyle = '--')
        # plt.plot([nodes1[0][0], koz_x], [koz_bot, koz_bot], color = 'red', linestyle = '--')
        # plt.plot([nodes1[0][0], koz_x], [koz_top, koz_top], color = 'red', linestyle = '--')
        # plt.plot([corridor, corridor], [nodes1[1][0], nodes2[1][2]], color = 'orange', linestyle = 'dashed')
        # # ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':')
        # ax.grid(True)
        # ax.axis('equal')
        # # ax.set_xlabel('X (ft)')
        # # ax.set_ylabel('Y (ft)')
        # ax.legend(loc = 'center left', fontsize = '8')
        # # def
        # plt.show()

    nodes1 = [np.array([nodes1[0][0], optim_sol1[0], nodes1[0][2]]).flatten(),np.array([nodes1[1][0], optim_sol1[1], nodes1[1][2]]).flatten()]
    print(curvature(P0 = [nodes1[0][0],nodes1[1][0]], P1 = [nodes1[0][1],nodes1[1][1]], P2 = [nodes1[0][2], nodes1[1][2]]) - turn_radius)
    nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]
    end_time = time.time()
    print('TOTAL BEZ OPTIMIZATION TIME:', end_time-start_time)
    # plt.plot(optimal_bez1[0], optimal_bez1[1])
    # plt.plot(optimal_bez2[0], optimal_bez2[1])
    # plt.axis('equal')
    # plt.scatter(nodes1[0], nodes1[1])
    # plt.scatter(nodes2[0], nodes2[1])    
    # plt.show()
    paths = [optimal_bez1, optimal_bez2]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    # print(x_wpts)
    # plt.scatter(x_wpts, y_wpts)
    
    return wpts, x_wpts, y_wpts, optimal_bez1, optimal_bez2, nodes1, nodes2 

    # x_entry, y_entry, central_angle, entryLength, entryTOA, h, k = entryPath(velocity, 29.49, 1353.34518698047)
    # t_start1 = 0.578786727511435
    # partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
    #                             [optim_sol1[0],optim_sol1[1]],
    #                             [nodes1[0][2],nodes1[1][2]], int(np.round(200*t_start1)), 0.601342376572023)

    # paths = [[x_entry, y_entry], partial_bez1]
    # wpts, x_wpts, y_wpts = paths_to_wp(paths, 5)

def EntryExitOutside(nodes1, nodes2, pos, velocity, lr, id):
    start = time.time()
    vel_knots = 111.6

    optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
    optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[nodes2[0][1],nodes2[1][1]],P2=[nodes2[0][2], nodes2[1][2]], t=1)
    # if id == 2:
    #     ba, t_entry = solve_optimEntry(np.deg2rad(10), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots, pos)
    # else:
    ba, t_entry = solve_optimEntry(np.deg2rad(20), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots, pos, id)
    # else:
    #     ba, t_entry = solve_optimEntry(np.deg2rad(26), np.deg2rad(30), np.deg2rad(20), nodes1, vel_knots, pos, lr)

    x_int_en, y_int_en = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry)
    x_int_en2, y_int_en2 = find_bez_xy([nodes1[0][0],nodes1[1][0]],
                                [nodes1[0][1],nodes1[1][1]],
                                [nodes1[0][2],nodes1[1][2]], t_entry+0.0025)
    
    req_ent = np.rad2deg(np.arctan2(y_int_en2-y_int_en,x_int_en2-x_int_en))

    
    x_entry, y_entry, central_angle, entryLength, entryTOA, h_c, k_c = entryPath(velocity, np.rad2deg(ba), [x_int_en, y_int_en], pos)
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

    optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [nodes2[0][1],nodes2[1][1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
    
    
    pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[nodes1[0][1],nodes1[1][1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_entry)
    
    total_toa = entryTOA+(pb1_length/velocity)+(optim2_length/velocity)
    print('TOTAL TOA:', total_toa)

    paths = [partial_bez1, optimal_bez2]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    plt.figure()
    plt.plot(partial_bez1[0], partial_bez1[1], color = 'magenta')
    plt.plot(optimal_bez2[0], optimal_bez2[1], color = 'magenta')
    plt.plot(x_entry, y_entry, color = 'cyan')
    plt.axis('equal')
    plt.scatter(x_wpts, y_wpts)
    plt.show()
    end = time.time()
    print('TOTAL INTERCEPTION PATH GENERATION TIME:', end-start)
    return [x_entry, y_entry, req_ent, np.rad2deg(ba)], total_toa, x_wpts, y_wpts

def entryPath(velocity, ba, intersect, pos):

    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
    tr*=0.3048
    # print('TR', tr)
    h, k  = tr, pos[1]

    x_entry = [i for i in np.linspace(0, intersect[0], 200)]

    y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
    # print(x_entry, y_entry)
    # y_entry[-1] = intersect[1]
    central_angle = pi - np.arctan2(y_entry[-1]-k, x_entry[-1] - h)
    entryLength = tr*central_angle #entryLength = 2*pi*tr * (central_angle/(2*pi))
    
    entryTOA = entryLength/velocity
    print('ENTRY ToA: ', entryTOA)

    return x_entry, y_entry, central_angle, entryLength, entryTOA, h, k




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


if __name__ == "__main__":
 
    velocity  = 57.3024 #m/s
    turn_rate = np.deg2rad(4.50) # RAD/s
    turn_rate = np.deg2rad(4.50) # RAD/s
    turn_radius = 111.6**2/(11.26*math.tan(np.deg2rad(30)))
    turn_radius*=0.3048
    h = 275
    koz_bot = 182.88
    koz_top = 60.96

    
    ev_toa = 10
    
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
        nodes1 = [np.array([229, 442, 450]).flatten(),np.array([0, h/60, h/2]).flatten()]
        nodes2 = [np.array([450, 465, 213.36]).flatten(),np.array([h/2, h, 137.16]).flatten()]
    # nodes1 = [np.array([750, 1450, 1475]).flatten(),np.array([0, h/20, 450]).flatten()]
    # nodes2 = [np.array([1475, 1400, 700]).flatten(),np.array([450, 1500, 450]).flatten()]
 
    print("VEHICLE VELOCITY:", velocity)
    # print("VEHICLE TURN RADIUS:", turn_radius)
    # print("TARGET LENGTH: ", target_length)
    valid1 = False
    valid2 = False
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    m_p21 = (nodes2[1][1] - nodes1[1][2])/(nodes2[0][1] - nodes1[0][2])
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
                                    target_heading=np.pi/2,
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
        optim1_length = path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1)
        optim2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]],P2=[nodes2[0][2], nodes2[1][2]], t=0.75)
        optim_toa = ((optim1_length+optim2_length) / velocity)
        print(optim_toa)

        if valid2 == False and optim_toa >= ev_toa:
            valid2 = True
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

    print('Entry Bank:', np.rad2deg(ba),'Entry t:', t_entry)


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
    
    # t_start2 = 0.35678391959798994

    # nodes2 = [np.array([nodes2[0][0], optim_sol2[0], nodes2[0][2]]).flatten(),np.array([nodes2[1][0], optim_sol2[1], nodes2[1][2]]).flatten()]
    # baEx, exit_t = solve_optimExit([np.deg2rad(30), t_start2], np.deg2rad(30), np.deg2rad(20), t_start2, nodes2, 111.6)
    # print('EXIT BANK', np.rad2deg(baEx), exit_t )

    # x_int_ex, y_int_ex = find_bez_xy([nodes2[0][0],nodes2[1][0]],
    #                             [optim_sol2[0],optim_sol2[1]],
    #                             [nodes2[0][2],nodes2[1][2]], exit_t)
    
    # x_exit, y_exit, central_angle_ex, exitLength, exitTOA, h_ex, k_ex = exitPath(velocity=188, t_exit=exit_t,
    #                                                                              ba = baEx, intersect=[x_int_ex, y_int_ex],
    #                                                                              nodes = nodes2)
    # x_exit[0] = 750
    # print(y_exit[0])
    # y_exit[0] = 2090
    # print(len(optimal_bez2[0]))
    head_ex = np.rad2deg(np.arctan2(optimal_bez2[1][-1]-optimal_bez2[1][198], optimal_bez2[0][-1]- optimal_bez2[0][198]))
    print('ACTUAL HEADING AT EXIT:', head_ex)
    # x_exit, y_exit, exit_int, exit_t, exitTOA, exit_bank, exitLength = exitPath2(velocity, t_exit, optimal_bez2)
    exitTOA = optim2_length/velocity
    # print('EXIT T:', exit_t)
    # partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
    #                             [optim_sol2[0],optim_sol2[1]],
    #                             [nodes2[0][2],nodes2[1][2]], 200, exit_t)
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
    # pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]], P2=[nodes2[0][2],nodes2[1][2]], t=t_start2)


    y = [i for i in np.linspace(koz_bot, koz_top, 100)]
    bx = [0 for i in y]
    bx2 = [305 for i in y]
    ybot = [koz_bot for i in range(305)]
    bxbot = [i for i in range(0, 305)]
    ytop = [koz_top for i in range(305)]
    y2 = [i for i in range(-396, h+100)]
    xwall = [0 for i in range(-396, h+100)]
    xwall2 = [457 for i in range(-396, h+100)]

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    # ax.scatter(x_wpts, y_wpts, zorder = 100, marker = '^')
    # ax.plot(x_exit, y_exit, color = 'orange', label = 'Exit Path')
    # ax.scatter(h_ex, k_ex, marker = 's', color = 'green', label = 'Exit Circle Center')
    # ax.scatter(x_exit[-1], y_exit[-1], color = 'red', marker = '*', s = 100, zorder = 100, label = 'Exit Point')
    # print('EXIT COORDS:', x_exit, y_exit)
    # ax.scatter(x_exit[0], y_exit[0], color = 'purple', marker = '*', s = 100, zorder = 100, label = 'Interception Point')
    # ax.scatter(pb2[0], pb2[1], c = 'magenta', zorder = 100)
    ax.plot([122, 305, 488], [h/2, h/2, h/2], linestyle = 'dashdot', alpha = 0.5)
    # ax.plot(optimal_bez2[0], optimal_bez2[1])
    ax.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial QBC')
    # ax.plot(partial_bez2[0], partial_bez2[1], c = 'magenta')#, label = 'Partial QBC')
    # ax.plot(optimal_bez1[0],optimal_bez1[1], c='black', label='Quadratic Bezier curve')
    # ax.scatter(x_wpts, y_wpts, label = 'Waypoints', marker = '^', color = 'green')
    # print(wpts_all1[0])
    ax.scatter(213.36, 137.16, marker = '*', color = 'green', label = 'Landing Pad Goal Point', zorder = 30)
    # ax.scatter(wpts_all1x, wpts_all1y, zorder = 30)
    print('ENTRY ANGLE:', np.rad2deg(np.arctan2(partial_bez1[1][1] - y_entry[-1], partial_bez1[0][1] - x_entry[-1])))
    ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
    
    ax.plot(optimal_bez2[0],optimal_bez2[1], c='black', label='Quadratic Bezier curve')

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
    print("TOTAL TOA:", entryTOA + exitTOA + pb1_length/velocity)


    
    ax.plot(x_slope, y_slope, color = 'purple', linestyle = 'dashed', label = 'G1 Continuity Line')
    ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
    ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
    ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
    ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
    # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
    ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
    ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
    ax.text(nodes1[0][2]-100,nodes1[1][2]-25,  r'Curve 1', fontsize = 6)
    ax.text(nodes1[0][2]-100,nodes1[1][2]+25,  r'Curve 2', fontsize = 6)

    # ax.text(mx+0.25, my, r'$\bf{m}$')
    # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
    # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
    # print(bx2, y)
    ax.plot(bx, y, color = 'red', linestyle = '--')
    ax.plot(bx2, y, color = 'red', linestyle = '--')

    ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
    ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
    ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':', color = 'orange')
    ax.plot(xwall2, y2, linestyle = ':', color = 'orange')
    # y_nom = [i for i in range(-1282, int(y_exit[0]))]
    # x_nom = [750 for i in y_nom]
    # ax.plot(x_nom, y_nom, color = 'black', label = 'Nominal Path')

    # ax.plot([700, 700], [-1300, y_exit[0]+50], color = 'green', label = 'Nominal Path')

    # ax.scatter(mx,my)
    # ax.scatter(center1x,center1y)
    # ax.scatter(center2x,center2y)
    # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
    # ax.add_patch(Circle((center2x, center2y), r2, color='black', fill=False))
 
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.legend(loc = 'center left', fontsize = '7')
    # plt.ylim(optim_sol1[1]-150, optim_sol2[1]+150)
    # plt.axes.set_ylim(-150, 600)
    plt.show()

#line types for, dashed, solid, dots
#