import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
 
def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
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



def solve_optim1latlon(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    if lr == 1:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: (P0[0]+0.002055) - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0] - (P0[0]+0.000686)},
                {'type': 'ineq', 'fun': lambda x: x[1] -P0[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
    else:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: (P0[0]-0.002055) - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0]+(P0[0]-0.000686)},
                {'type': 'ineq', 'fun': lambda x: x[1] -P0[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
    # print(val)
 
    return val.x , curvature(P0,val.x,P2)

def solve_optim2latlon(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr):#turn_radius,
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
  
    if lr == 1:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: (P0[0]+0.0002055) - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0] - (P0[0]+0.000686)},
                {'type': 'ineq', 'fun': lambda x: x[1] - P0[1]},
                {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
    else:
        cons = (
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
                {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
                {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
                {'type': 'ineq', 'fun': lambda x: (P0[0]-0.002055) - x[0]},
                # {'type': 'ineq', 'fun': lambda x: np.deg2rad(20) - np.abs(np.arctan2(x[1]-P0[1], x[0] - P0[0]))},
                {'type': 'ineq', 'fun': lambda x: x[0]+(P0[0]-0.000686)},
                {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
                # {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P0[1]-x[1]), (P0[0]-x[0])))},
                {'type': 'ineq', 'fun': lambda x: np.abs(target_heading-np.arctan2((P2[1]-x[1]), (P2[0]-x[0])))}
                ) 
   #
    
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
 
    # print(val)
    # enumerate()
    return val.x , curvature(P0,val.x,P2)

# def solve_optimPrelatlon(P0, P2, target_toa,  guess, target_heading, velocity, turn_radius, lr):#turn_radius,
#     def path_cost(P1):
#         return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
#     if lr == 1:
#         cons = (
                
#                 {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
#                 {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
#                 {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
#                 # {'type': 'ineq', 'fun': lambda x: P0[0]+x[0]},
#                 {'type': 'ineq', 'fun': lambda x: x[1]-P0[1]},
#                 # {'type': 'ineq', 'fun': lambda x: np.deg2rad(10) - np.abs(np.arctan2(x[1]-P2[1], x[0] - P2[0]))},
#                 # {'type': 'ineq', 'fun': lambda x: x[0] - 1000},
#                 {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
#                 {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P2[1]-x[1]), (P2[0]-x[0]))-target_heading)}
#             ) 
#     else:
#         cons = (
                
#                 {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
#                 {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
#                 {'type': 'ineq', 'fun': lambda x: P0[0]-x[0]},
#                 # {'type': 'ineq', 'fun': lambda x: P0[0]+x[0]},

#                 # {'type': 'ineq', 'fun': lambda x: 0+x[0]},
#                 # {'type': 'ineq', 'fun': lambda x: np.deg2rad(10) - np.abs(np.arctan2(x[1]-P2[1], x[0] - P2[0]))},
#                 # {'type': 'ineq', 'fun': lambda x: 500 - x[0]},
#                 {'type': 'ineq', 'fun': lambda x: x[1]-P0[1]},
#                 {'type': 'ineq', 'fun': lambda x: P2[1] - x[1]},
#                 {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P2[1]-x[1]), (P2[0]-x[0]))-target_heading)}
#             ) 
#     val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)

#     return val.x , curvature(P0,val.x,P2)

def bez_to_wp(bez1, bez2, num):
    path_wpts = []
    point_marker = len(bez1[0])/num
    point_marker = np.round(point_marker)
    # for i in range(len(bez3[0])):
    #     if i%point_marker == 0:
    #         path_wpts.append([bez3[1][i], bez3[0][i]]) 
    for i in range(len(bez1[0])):
        if i%point_marker == 0:
            path_wpts.append([bez1[1][i], bez1[0][i]]) #building list of waypoints FLIPPED FOR LAT LON
    for i in range(len(bez2[0])):
        if i%point_marker == 0:
            path_wpts.append([bez2[1][i], bez2[0][i]]) #building list of waypoints FLIPPED FOR LAT LON
    wpts_x = [x[1] for x in path_wpts]
    wpts_y = [y[0] for y in path_wpts]
    return path_wpts, wpts_x, wpts_y

def bez_to_wp_single(bez1, num):
    path_wpts = []
    point_marker = len(bez1[0])/num
    point_marker = np.round(point_marker)
    for i in range(len(bez1[0])):
        if i%point_marker == 0:
            path_wpts.append([bez1[0][i], bez1[1][i]]) #building list of waypoints
    wpts_x = [x[0] for x in path_wpts]
    wpts_y = [y[1] for y in path_wpts]
    return path_wpts, wpts_x, wpts_y


def validity_checklatlon(x, y1, y2, points, corridor, curve, lr):
    points = np.array(points)
    print('LR IS: ', lr)
    if lr == 1:
        for i in range(len(points)):
            if points[i][0]<=x and points[i][1]>=y1 and points[i][1]<=y2 and points[i][0] >= corridor:
                print('In Positive')
                print(curve)
                print(points[i][0], x, points[i][0]<=x)
                print(points[i][1], y1, points[i][1]>=y1)
                print(points[i][1], y2, points[i][1]<=y2)
                print(points[i][0], corridor, points[i][0] >= corridor)
                plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    elif lr == -1:
        for i in range(len(points)):
            if points[i][0]>=x and points[i][1]>=y1 and points[i][1]<=y2 and points[i][0] <= corridor:
                print('In Negative')
                print(curve)
                print(points[i][0], x, points[i][0]>=x)
                print(points[i][1], y1, points[i][1]>=y1)
                print(points[i][1], y2, points[i][1]<=y2)
                print(points[i][0], corridor, points[i][0] <= corridor)
                plt.scatter(points[i][0], points[i][1], marker = '*', s = 150, zorder = 50)
                return False
    return True

def toCallOutside(velocity, turn_rate, target_toa1, target_toa2, uav_head, nodes1, nodes2, koz_x, koz_bot, koz_top, lr):
    turn_radius = velocity / turn_rate
    valid1 = False
    valid2 = False
    if lr == 1:
        corridor = nodes1[0][0]+0.002055
    else:
        corridor = nodes1[0][0]-0.002055
    c = 0
    while not valid1 or not valid2:
        optim_sol1, curv1 = solve_optim1latlon(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                    target_toa=target_toa1,
                                    guess=[nodes1[0][1], nodes1[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        optim_sol2, curv2 = solve_optim2latlon(P0=[nodes2[0][0],nodes2[1][0]],P2=[nodes2[0][2], nodes2[1][2]],
                                    target_toa=target_toa2,
                                    guess=[nodes2[0][1], nodes2[1][1]],
                                    target_heading=np.pi/2,
                                    velocity=velocity, turn_radius=turn_radius, lr=lr)
        # optim_solPre, curvPre = solve_optimPrelatlon(P0=[nodesPre[0][0],nodesPre[1][0]],P2=[nodesPre[0][2], nodesPre[1][2]],
        #                             target_toa=1.91,#target_toa1/target_toa1,
        #                             guess=[nodesPre[0][1], nodesPre[1][1]],
        #                             target_heading=np.pi/4,
        #                             velocity=velocity, turn_radius=turn_radius, lr=lr)
        
        # print(optim_sol)
        
        optimal_bez1 = manual_bez([nodes1[0][0],nodes1[1][0]],
                                [optim_sol1[0],optim_sol1[1]],
                                [nodes1[0][2],nodes1[1][2]], 200)
        
        solved_heading = np.arctan2((optim_sol1[1]-nodes2[1][0]),(optim_sol1[1]-nodes2[0][0]))
        optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [optim_sol2[0],optim_sol2[1]],
                                [nodes2[0][2],nodes2[1][2]], 200)
        # optimal_bezPre = manual_bez([nodesPre[0][0],nodesPre[1][0]],
        #                         [optim_solPre[0],optim_solPre[1]],
        #                         [nodesPre[0][2],nodesPre[1][2]], 200)
        # print(len(optimal_bez1[0]))

        wpts, x_wpts, y_wpts = bez_to_wp(optimal_bez1, optimal_bez2, 15) 
        wpts_all1, wpts_all1x, wpts_all1y = bez_to_wp_single(optimal_bez1, 30) 
        # print(wpts_all1)
        wpts_all2, wpts_all2x, wpts_all2y = bez_to_wp_single(optimal_bez2, 30) 
        valid1 = validity_checklatlon(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
        valid2 = validity_checklatlon(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
        # print(optim_sol1[0])
        # print(target_toa1)

        # print(valid1, valid2)
        if valid1 == False:
            target_toa1+=0.1
        if valid2 == False:
            target_toa2+=0.1
        # print(target_toa2)
        # #     c+=1
        # # if c>=5:
        # # y = [i for i in range(10)]
        # # bx = [koz_x for i in range(10)]
        # # ybot = [koz_bot for i in range(250)]
        # # bxbot = [i for i in range(39)]
        # ytop = [koz_top for i in range(2)]
        # # y2 = [i for i in range(-50, 950)]
        # # xwall = [1500 for i in range(-50, 950)]
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.scatter([nodes1[0][0], optim_sol1[0], nodes1[0][2]],[nodes1[1][0],optim_sol1[1],nodes1[1][2]], label='Bezier Curve 1 Control Points')
    ax.plot(optimal_bez1[0],optimal_bez1[1], c='black',label='Quadratic Bezier curve')
    ax.plot(optimal_bez2[0],optimal_bez2[1], c='black')
    # ax.plot(optimal_bezPre[0], optimal_bezPre[1], c = 'black')
    valid1 = validity_checklatlon(koz_x, koz_bot, koz_top, wpts_all1, corridor, 'Curve 1', lr)
    valid2 = validity_checklatlon(koz_x, koz_bot, koz_top, wpts_all2, corridor, 'Curve 2', lr)
    print([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]])
    ax.scatter([nodes2[0][0], optim_sol2[0], nodes2[0][2]],[nodes2[1][0],optim_sol2[1],nodes2[1][2]], label='Bezier Curve 2 Control Points')
    # ax.scatter(optim_solPre[0], optim_solPre[1], marker = '*')
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
    # ax.plot([-83.2, -83.199], ytop, color = 'red', label = 'Emergency Vehicle Clearance Area', linestyle = '--')
    # ax.plot(xwall, y2, label = 'Flight Corridor Bound', linestyle = ':')
    ax.grid(True)
    ax.axis('equal')
    # ax.set_xlabel('X (ft)')
    # ax.set_ylabel('Y (ft)')
    ax.legend(loc = 'center left', fontsize = '8')
    # def
    plt.show()
    return wpts, x_wpts, y_wpts

