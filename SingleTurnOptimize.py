import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
from scipy.optimize import fsolve
import math
import time
import scipy
import sympy as sp
from scipy.optimize import brentq
import Scenario1Optimize as MBO
from algebra_with_sympy import *

'''
TO DO: 



'''
def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

'''ENTRY TURN STUFF'''
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

def find_diff_entry(ba, nodes, velocity):
    pi = np.pi
    t_guess = 0.25
    mindiff = 100
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    tr = velocity**2 / (11.26*math.tan(ba[0]))
    h = tr + 750
    k = -1282
    circle_eq = lambda t: (Bx(t)-h)**2+(By(t)-k)**2 - tr**2

    S = fsolve(circle_eq, t_guess)
    t_final = 0

    for i in S:

        if i>=0 and i <= 1:
            index = int(np.round(float(i*200)))
            if index == 200:
                index = 199
            bez_angle = np.arctan2(By(i)-By(i-0.01), Bx(i)-Bx(i-0.01))

            x_l = [i for i in np.linspace(750, Bx(i), 200)] #Space between nominal path and bez
            y = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            # plt.plot(path[0], path[1])
            # plt.plot(x_l, y)
            # plt.axis('equal')
            # plt.show()

            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            diff = np.abs(bez_angle-int_angle)
            if diff < mindiff and np.abs(y[-1] - By(i))<=.00001:
                mindiff = diff
                # plt.plot(path[0], path[1])
                # plt.plot(x_l, y)
                # plt.scatter(path[0][index], path[1][index], color = 'yellow', marker = '*', s = 100)
                # plt.scatter(x_l[-1], y[-1], color = 'purple', marker = '*', s = 100)
                # plt.axis('equal')
                t_final = i
            else:
                t_final = 0

    return [mindiff, t_final]

def entryTurn(velocity, path, start, nodes):
    pi = np.pi
    t = sp.symbols('t', real = True)
    Bx = nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    sols = []
    
    # for t_index in np.linspace(0, 200, 1):
    for ba in np.linspace(25, 30, 100):
        '''Circle Parameters'''
        tr = velocity**2 / (11.26*np.tan(np.deg2rad(ba)))
        h = tr + start[0]
        k = start[1]
        circle_eq = (Bx-h)**2+(By-k)**2 - tr**2
        eqs = [circle_eq, t>0, t<1]
        S = sp.nonlinsolve(eqs,t)

        for i in S:
            if i[0].is_real:
                sols.append([ba, i[0]])
                
    return sols
        
'''EXIT TURN STUFF'''
def solve_optimExit(guess, max_bank, min_bank, min_t, nodes, velocity):
    #t is also being minimized, so part of constraints
    def path_cost(guess, nodes, velocity):
        diff = find_diff_exit(guess, nodes, velocity)
        print('DIFF FOR EXIT:', diff[0])
        return np.abs(diff[0]) 

    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank},
            {'type': 'ineq', 'fun': lambda x: x[1] - min_t}
    )
    val =  minimize(path_cost, guess, (nodes, velocity), method = 'SLSQP', tol = 1E-10, constraints=cons)
    return val.x

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

    tr = velocity**2 / (11.26*math.tan(ba))
    h = lambda t: Bx(t) - tr*math.cos(bezHead)
    k = lambda t: By(t) + tr*math.sin(bezHead)

    circle_eq = lambda y: (750-h(t_guess))**2+(y-k(t_guess))**2 - tr**2
    y_guess = 2000
    S = fsolve(circle_eq, y_guess) #gives y intersection
    t_final = .5

    y_l = [i for i in np.linspace(900, S, 200)]
    x = [750 for i in y_l]

    # print(S)
    for i in S:
        if i > 0: 
            y = [i for i in np.linspace(S, By(t_guess))]
            x_l = [h(t_guess) - np.sqrt(tr**2 - (y_y - k(t_guess))**2) for y_y in y]
            if x_l[0] ==750:# <= 0.01:
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
                    print('BLARGY', diff)
                    t_final = t_guess
                

    return [mindiff, t_final]
        
if __name__=='__main__':
    start = time.time()
    # wpts, x_wpts, y_wpts, bez1, bez2 = MBO.toCallOutside(188, np.deg2rad(20), 1.5*3.71, 1.5*3.71, np.deg2rad(90), 
    #                                     [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()],
    #                                     [np.array([1475, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()],
    #                                     1000, 0, 850, 1)
    
    # Sols = entryTurn(111.6, bez1, [750, -1282], [np.array([750, 1461.58471522, 1475]).flatten(),np.array([-50, -50, 450]).flatten()])
    # print(Sols)
    # end = time.time()
    # print(end-start)
    # print(np.deg2rad(29.22))
    velocity = 111.6
    nodes1 = [np.array([750, 1461.58471522, 1475]).flatten(),np.array([-50, -50, 450]).flatten()]
    optim_sol, t_val= solve_optimEntry(np.deg2rad(25), np.deg2rad(30), np.deg2rad(25), nodes1, velocity)

    nodes2 = [np.array([1475, 1500, 750]).flatten(),np.array([450, 900, 900]).flatten()]
    t_start2 = 0.35678391959798994
    optim_sol2 = solve_optimExit([np.deg2rad(30), 0.3], np.deg2rad(30), np.deg2rad(25), t_start2, nodes2, velocity)
    print(np.rad2deg(optim_sol2[0]), optim_sol2[1])
    print('DIFF', np.pi/2 - optim_sol2[0])
    end = time.time()
    print(end-start)
    # print(np.rad2deg(optim_sol[0]), t_val)
   