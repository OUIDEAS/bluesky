import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
 
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
 
 
def solve_optim(P0, P2, target_toa, turn_radius, guess, target_heading, velocity):
    def path_cost(P1):
        return (np.abs(path_length(P1, P0, P2, 1) - target_toa*velocity))
    
    
    cons = (
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2) - turn_radius},
            {'type': 'ineq', 'fun': lambda x: curvature(P0,x,P2)},
            {'type': 'ineq', 'fun': lambda x: x[0] - P0[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] + P2[0]},
            {'type': 'ineq', 'fun': lambda x: x[1] - P2[1]},
            {'type': 'ineq', 'fun': lambda x: np.abs(np.arctan2((P0[1]-x[1]), (P0[0]-x[0]))-target_heading)}
            )  
   
    val = minimize(path_cost,[guess[0],guess[1]], method='SLSQP', tol=1E-10, constraints=cons)
 
    print(val)
 
    return val.x , curvature(P0,val.x,P2)
 
 
 
if __name__ == "__main__":
 
    velocity  = 100
    turn_rate = 111111 # RAD/s
    turn_radius = velocity / turn_rate
   
    target_toa = 5
    uav_head = np.deg2rad(60)
 
 
    # target_length = target_toa/velocity
 
    print("VEHICLE VELOCITY:", velocity)
    print("VEHICLE TURN RADIUS:", turn_radius)
    # print("TARGET LENGTH: ", target_length)
 
    # GOAL: DESIGN P1 TO MEET CONSTRAINTS
    nodes1 = [np.array([0, 100, 400]).flatten(),np.array([0, 100, 100]).flatten()]
 
    # curve1 = bezier.Curve(nodes1, degree=len(nodes1[0])-1)
    bezier_moleski = manual_bez([nodes1[0][0],nodes1[1][0]], [nodes1[0][1],nodes1[1][1]],[nodes1[0][2],nodes1[1][2]], 200)
 
    #OPTIMIZE TEST:
    optim_sol, curv = solve_optim(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
                                  target_toa=target_toa,
                                  turn_radius=turn_radius,
                                  guess=[nodes1[0][1], nodes1[1][1]],
                                  target_heading=np.pi/4,
                                  velocity=velocity)
   
    optimal_bez = manual_bez([nodes1[0][0],nodes1[1][0]],
                             [optim_sol[0],optim_sol[1]],
                             [nodes1[0][2],nodes1[1][2]], 200)
 
    solved_heading = np.arctan2((optim_sol[1]-nodes1[1][0]),(optim_sol[1]-nodes1[0][0]))
 
    print("HEADING TO P1", solved_heading)
    print("REQUESTED HEADING: ", uav_head)
 
    print("SOLVED LENGTH: ", path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol[0],optim_sol[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1))
    print("UAV TURN RADIUS = ", turn_radius, "WHILE MIN RADIUS OF CURVATURE = ", curv)
 
    mx = np.mean([nodes1[0][0], nodes1[0][2]])
    my = np.mean([nodes1[1][0], nodes1[1][2]])
 
    center1x = np.mean([mx, nodes1[0][0]])
    center1y = np.mean([my, nodes1[1][0]])
    r1 = np.sqrt(
        (center1x-mx)**2 + (center1y-my)**2
    )
 
    center2x = np.mean([mx, nodes1[0][2]])
    center2y = np.mean([my, nodes1[1][2]])
    r2 = np.sqrt(
        (center2x-mx)**2 + (center2y-my)**2
    )
    circ1 = [center1x, center1y, r1]
    circ2 = [center2x, center2y, r2]
 
 
    print("ToA:", path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol[0],optim_sol[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1) / velocity)
 
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
 
    ax.plot(optimal_bez[0],optimal_bez[1], c='black',label='Quadratic Bezier curve', linestyle = '--')
 
    # ax.scatter([nodes1[0][0], optim_sol[0], nodes1[0][2]],[nodes1[1][0],optim_sol[1],nodes1[1][2]], label='Control Points')
    # nodes1 = [np.array([0, 100, 350]).flatten(),np.array([0, 100, 150]).flatten()]
 
    # # curve1 = bezier.Curve(nodes1, degree=len(nodes1[0])-1)
    # bezier_moleski = manual_bez([nodes1[0][0],nodes1[1][0]], [nodes1[0][1],nodes1[1][1]],[nodes1[0][2],nodes1[1][2]], 200)
 
    # #OPTIMIZE TEST:
    # optim_sol, curv = solve_optim(P0=[nodes1[0][0],nodes1[1][0]],P2=[nodes1[0][2], nodes1[1][2]],
    #                               target_toa=target_toa,
    #                               turn_radius=turn_radius,
    #                               guess=[nodes1[0][1], nodes1[1][1]],
    #                               target_heading=np.pi/4,
    #                               velocity=velocity)
   
    # optimal_bez = manual_bez([nodes1[0][0],nodes1[1][0]],
    #                          [optim_sol[0],optim_sol[1]],
    #                          [nodes1[0][2],nodes1[1][2]], 200)
    ax.plot([nodes1[0][0], nodes1[0][2]], [nodes1[1][0], nodes1[1][2]], label = 'Straight Line Path')
    ax.scatter([nodes1[0][0], optim_sol[0], nodes1[0][2]],[nodes1[1][0],optim_sol[1],nodes1[1][2]], label='Control Points')
    # ax.plot(optimal_bez[0],optimal_bez[1], c='black',label='Quadratic Bezier curve')
    # print("ToA:", path_length(P0=[nodes1[0][0],nodes1[1][0]], P1=[optim_sol[0],optim_sol[1]],P2=[nodes1[0][2], nodes1[1][2]], t=1) / velocity)

    ax.text(nodes1[0][0]+0.25,nodes1[1][0],  r'$\bf{p_0}$')
    ax.text(optim_sol[0]+0.25,optim_sol[1],  r'$\bf{p_1}$')
    ax.text(nodes1[0][2]+0.25,nodes1[1][2],  r'$\bf{p_2}$')
    ax.text(nodes1[0][2]-75,nodes1[1][2]+50,  r'Bezier TOA=5')

    # dist = np.hypot(nodes1[1][0]- nodes1[1][2], nodes1[0][0] - nodes1[0][2] )
    # t = dist/velocity
    # print(t)
    ax.text(nodes1[0][0],nodes1[1][0]-20,  r'Straight Line TOA=4.12')
    # ax.text(mx+0.25, my, r'$\bf{m}$')
    # ax.text(center1x+0.25, center1y, r'$\bf{C_1}$')
    # ax.text(center2x+0.25, center2y, r'$\bf{C_2}$')
 
    # ax.scatter(mx,my)
    # ax.scatter(center1x,center1y)
    # ax.scatter(center2x,center2y)
    # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
    # ax.add_patch(Circle((center2x, center2y), r2, color='black', fill=False))
 
    ax.grid(True)
    ax.axis('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.legend()
    # def
    plt.show()