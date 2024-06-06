import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy
import sympy as sp
from scipy.optimize import brentq

def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

def manual_bez_xy(P0, P1, P2, points):
    t = np.linspace(0, 1, points)
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

 
def solve_optim1(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    if lr == 1:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: 1500 - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0] - 1000},
                {'type': 'ineq', 'fun': lambda x: x[1] -P0[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
    else:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0]-P0[0]},
                {'type': 'ineq', 'fun': lambda x: 0-x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: 500 - x[0]},
                {'type': 'ineq', 'fun': lambda x: x[1] -P0[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
 
    return val.x , curvature(P0,val.x,P2)

def solve_optim2(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    if lr == 1:
        cons = (
                
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: 1500 - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(10) - np.abs(np.arctan2(x[1]-P2[1], x[0] - P2[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0] - 1000},
                {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
                {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P2[1]-x[1]), (P2[0]-x[0]))-target_heading)}
            ) 
    else:
        cons = (
                
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: 0-x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(10) - np.abs(np.arctan2(x[1]-P2[1], x[0] - P2[0]))},
                {'type': 'ineq', 'fun': lambda x: 500 - x[0]},
                {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
                {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P2[1]-x[1]), (P2[0]-x[0]))-target_heading)}
            ) 
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

    return val.x , curvature(P0,val.x,P2)



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
                print('In Positive')
                print(curve)
                print(points[i][0], x, points[i][0]<=x)
                print(points[i][1], y1, points[i][1]>=y1)
                print(points[i][1], y2, points[i][1]<=y2)
                print(points[i][0], corridor, points[i][0] >= corridor)
                # plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    elif lr == -1:
        for i in range(len(points)):
            if points[i][0]>=x and points[i][1]>=y1 and points[i][1]<=y2 or points[i][0] <= corridor:
                print('In Negative')
                print(curve)
                print(points[i][0], x, points[i][0]>=x)
                print(points[i][1], y1, points[i][1]>=y1)
                print(points[i][1], y2, points[i][1]<=y2)
                print(points[i][0], corridor, points[i][0] <= corridor)
                # plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    return True


def toCallOutside(velocity, turn_rate, target_toa1, target_toa2, uav_head, nodes1, nodes2, koz_x, koz_bot, koz_top, lr):
    turn_radius = velocity / turn_rate
    valid1 = False
    valid2 = False
    if lr == 1:
        # corridor = nodes1[0][0]+0.002055
        corridor = 1500
    else:
        corridor = -1500
        # corridor = nodes1[0][0]-0.002055
    c = 0
    while not valid1 or not valid2:
        optim_sol1, curv1 = solve_optim1(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=target_toa1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        optim_sol2, curv2 = solve_optim2(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=target_toa2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        # print(optim_sol)
        
        optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
        
        solved_heading = np.arctan2((optim_sol1[1]-nodes2[1][0]),(optim_sol1[1]-nodes2[0][0]))
        optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
        # print(len(optimal_bez1[0]))

        wpts, x_wpts, y_wpts = bez_to_wp(optimal_bez1, optimal_bez2, 15) 
        wpts_all1, wpts_all1x, wpts_all1y = bez_to_wp_single(optimal_bez1, 30) 
        # print(wpts_all1)
        wpts_all2, wpts_all2x, wpts_all2y = bez_to_wp_single(optimal_bez2, 30) 
        valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
        if valid1 == False:
            target_toa1+=0.1
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            c+=1
        if c>5:
            
        # print(target_toa2)
        # #     c+=1
        # # if c>=5:
            y = [i for i in range(10)]
            bx = [koz_x for i in range(10)]
            ybot = [koz_bot for i in range(250)]
            bxbot = [i for i in range(39)]
            ytop = [koz_top for i in range(2)]
            y2 = [i for i in range(-50, 950)]
            xwall = [1500 for i in range(-50, 950)]
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
            ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
            ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
            valid1 = validity_check(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
            valid2 = validity_check(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
            print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
            ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
            ax.plot([corridor, corridor, corridor], [nodes1[1][0], nodes1[1][2], nodes2[1][2]], linestyle = '--')
            ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
            ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
            ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
            # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
            ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
            ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
            ax.text(nodes1[0][2]-250,nodes1[1][2]-100,  r'Curve 1')
            ax.text(nodes1[0][2]-250,nodes1[1][2]+100,  r'Curve 2')

            # ax.text(mx+0.25, my, r'$\bf{m}$')
            # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
            # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
            # ax.plot(bx, y, color = 'red', linestyle = '--')
            # ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
            ax.plot([-83.2, -83.199], ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
            # ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':')
            ax.grid(True)
            ax.axis('equal')
            # ax.set_xlabel('X (ft)')
            # ax.set_ylabel('Y (ft)')
            ax.legend(loc = 'center left', fontsize = '8')
            # def
            plt.show()

    x_entry, y_entry, central_angle, entryLength, entryTOA, h, k = entryPath(velocity, 29.49, 1353.34518698047)
    t_start1 = 0.578786727511435
    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], int(np.round(200*t_start1)), 0.601342376572023)
    
    t_start2 = 0.35678391959798994
    x_exit, y_exit, exit_int, exit_t, exitTOA, exit_bank, exitLength = exitPath2(velocity, t_start2, optimal_bez2)
    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], int(np.round(200*t_start2)), t_start2)
    paths = [[x_entry, y_entry], partial_bez1, partial_bez2, [x_exit, y_exit]]
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 5)


    return wpts, x_wpts, y_wpts, optimal_bez1, optimal_bez2 

def entryPath(velocity, ba, intersect):
    pi = np.pi
    '''Entry Into Bezier Curve'''
    tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
    h, k  = tr+750, -1282

    x_entry = [i for i in np.linspace(750, intersect, 200)]

    y_entry = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_entry]
    
    central_angle = pi - np.arctan2(y_entry[-1]-k, x_entry[-1] - h)
    entryLength = 2*pi*tr * (central_angle/(2*pi))
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
    
    
    mindex = exitTOA_l.index(min(exitTOA_l))
    tFinal = int(np.round(t_vals[mindex]*200))

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

if __name__ == "__main__":
 
    velocity  = 188 #ft/s
    turn_rate = np.deg2rad(20) # RAD/s
    turn_radius = velocity / turn_rate
    koz_bot = 0
    koz_top = 850
    
    c=0
    target_toa1 = 1.5*3.71
    target_toa2 = 1.5*3.71
    uav_head = np.deg2rad(90)
    lr = 1
    if lr == -1:
        corridor = 0
        koz_x = 500
        nodes1 = [np.array([750, 100, 25]).flatten(),np.array([-50, 0, 450]).flatten()]
        nodes2 = [np.array([25, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()]
    else:
        corridor = 1500
        koz_x = 1000
        nodes1 = [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()]
        nodes2 = [np.array([1475, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()]
    # target_length = target_toa/velocity
 
    print("VEHICLE VELOCITY:", velocity)
    # print("VEHICLE TURN RADIUS:", turn_radius)
    # print("TARGET LENGTH: ", target_length)
    valid1 = False
    valid2 = False
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    
    bezier_variny_1 = manual_bez([nodes1[0][0],nodes1[1][0]], [nodes1[0][1],nodes1[1][1]],[nodes1[0][2],nodes1[1][2]], 20)
    bezier_variny_2 = manual_bez([nodes2[0][0],nodes2[1][0]], [nodes2[0][1],nodes2[1][1]],[nodes2[0][2],nodes2[1][2]], 20)
    #OPTIMIZE TEST:
    start = time.time()
    while not valid1 or not valid2:

        optim_sol1, curv1 = solve_optim1(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=target_toa1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        optim_sol2, curv2 = solve_optim2(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=target_toa2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        
        
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
        print("OPTIM SOL: ", optim_sol1[0], optim_sol1[1])
        print(target_toa1)

        print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.1
            c+=1
        if valid2 == False:
            target_toa2+=0.1
            c+=1
        if c>5:
       
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
            ax.plot([corridor, corridor, corridor], [nodes1[1][0], nodes1[1][2], nodes2[1][2]], linestyle = '--')
            ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
            ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
            ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
            # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
            ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
            ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
            ax.text(nodes1[0][2]-250,nodes1[1][2]-100,  r'Curve 1')
            ax.text(nodes1[0][2]-250,nodes1[1][2]+100,  r'Curve 2')

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

    x_entry, y_entry, central_angle, entryLength, entryTOA, h, k = entryPath(velocity, 29.49, 1353.34518698047)
    t_start1 = 0.578786727511435
    partial_bez1 = manual_bez_partial([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200, 0.601342376572023)
    
    t_start2 = 0.35678391959798994
    x_exit, y_exit, exit_int, exit_t, exitTOA, exit_bank, exitLength = exitPath2(velocity, t_start2, optimal_bez2)
    partial_bez2 = manual_bez_partialExit([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200, t_start2)
    # print(x_exit)
    # print(x_entry[0])
    entry_path = [x_entry, y_entry]
    # print(entry_path[0][0])
    exit_path = [x_exit, y_exit]
    paths = [entry_path, partial_bez1, partial_bez2, exit_path]
    # print(len(paths))
    # print(paths[0][0])
    wpts, x_wpts, y_wpts = paths_to_wp(paths, 20)
    # print(wpts)


    pb1_length = optim1_length - path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol1[0],optim_sol1[1]],P2=[nodes1[0][2], nodes1[1][2]], t=t_start1)
    pb2_length = path_length(P0=[nodes2[0][0],nodes2[1][0]], P1=[optim_sol2[0],optim_sol2[1]], P2=[nodes2[0][2],nodes2[1][2]], t=t_start2)

    y = [i for i in range(850)]
    bx = [500 for i in range(850)]
    bx2 = [1000 for i in range(850)]
    ybot = [0 for i in range(500)]
    bxbot = [i for i in range(500, 1000)]
    ytop = [850 for i in range(500)]
    y2 = [i for i in range(-50, 950)]
    xwall = [0 for i in range(-50, 950)]
    xwall2 = [1500 for i in range(-50, 950)]
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    # ax.scatter(x_wpts, y_wpts, zorder = 100, marker = '^')
    ax.plot(x_exit, y_exit, color = 'orange', label = 'Exit Path')
    ax.scatter(x_exit[-1], y_exit[-1], color = 'purple', marker = '*', s = 100, zorder = 100)
    # ax.scatter(pb2[0], pb2[1], c = 'magenta', zorder = 100)
    ax.plot([400, 1000, 1600], [450, 450, 450], linestyle = 'dashdot', alpha = 0.5)
    ax.plot(partial_bez1[0], partial_bez1[1], c = 'magenta', label = 'Partial QBC')
    ax.plot(partial_bez2[0], partial_bez2[1], c = 'magenta')#, label = 'Partial QBC')
    # ax.plot(optimal_bez1[0],optimal_bez1[1], c='black', label='Quadratic Bezier curve')
    # ax.scatter(x_wpts, y_wpts, label = 'Waypoints', marker = '^', color = 'green')
    # print(wpts_all1[0])
    # ax.scatter(wpts_all1x, wpts_all1y, zorder = 30)
    print('ENTRY ANGLE:', np.rad2deg(np.arctan2(partial_bez1[1][1] - y_entry[-1], partial_bez1[0][1] - x_entry[-1])))
    ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
    
    # ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')

    # plt.scatter(h, k, label = 'Entry Turn Radius Center')


    ax.plot(x_entry, y_entry, label = 'Entry Arc', color = 'cyan')
    ax.scatter(x_entry[-1], y_entry[-1], marker = '*', color = 'purple', label = 'Intersection Point', s = 100, zorder = 100)
    # print(y_circ[0])
    

    print("ENTRY TOA:", entryTOA)
    print("EXIT TOA:", exitTOA)
    print("SOLVED BEZIER LEGNTH: ", pb1_length + pb2_length)
    print("PARTIAL BEZIER TRAVEL TIME:", (pb1_length+pb2_length)/velocity)

    print("SOLVED LENGTH WITH ENTRY/EXIT: ", pb1_length + optim2_length + entryLength+exitLength)
    print("TOTAL TOA:", entryTOA + exitTOA + (pb1_length+pb2_length)/velocity)


    

    ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
    ax.text(nodes1[0][0]+0.25,nodes1[1][0]-30,  r'$\bf{p_{01}}$')
    ax.text(optim_sol1[0]+0.25,optim_sol1[1],  r'$\bf{p_{11}}$')
    ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_{21}/p_{02}}$')
    # ax.text(nodes2[0][0]+0.25,nodes2[1][0],  r'$\bf{p_0}$')
    ax.text(optim_sol2[0]+0.25,optim_sol2[1],  r'$\bf{p_{12}}$')
    ax.text(nodes2[0][2]+0.25,nodes2[1][2]+20,  r'$\bf{p_{22}}$')
    ax.text(nodes1[0][2]-250,nodes1[1][2]-100,  r'Curve 1')
    ax.text(nodes1[0][2]-250,nodes1[1][2]+100,  r'Curve 2')

    # ax.text(mx+0.25, my, r'$\bf{m}$')
    # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
    # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
    ax.plot(bx, y, color = 'red', linestyle = '--')
    ax.plot(bx2, y, color = 'red', linestyle = '--')

    ax.plot(bxbot, ybot, color = 'red', linestyle = '--')
    ax.plot(bxbot, ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
    ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':', color = 'orange')
    ax.plot(xwall2, y2, linestyle = ':', color = 'orange')
    y_nom = [i for i in range(-1282, int(y_exit[0]))]
    x_nom = [750 for i in y_nom]
    ax.plot(x_nom, y_nom, color = 'black', label = 'Nominal Path')
    # ax.scatter(mx,my)
    # ax.scatter(center1x,center1y)
    # ax.scatter(center2x,center2y)
    # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
    # ax.add_patch(Circle((center2x, center2y), r2, color='black', fill=False))
 
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('X (ft)')
    ax.set_ylabel('Y (ft)')
    ax.legend(loc = 'lower left', fontsize = '8')
    # def
    plt.show()

#line types for, dashed, solid, dots
#