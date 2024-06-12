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
import MichaelBezierOptimize as MBO
from algebra_with_sympy import *

def solve_optimEntry(guess, max_bank,  min_bank):
    #controlling interception heading so path cost is heading?
    #Heading at t vs circle heading
    #Ideally diff is 0 -> target
    #So the diff closest to 0?
    def path_cost(ba):
        diff = find_diff_entry(ba)
        return np.abs(diff)
    cons = (
            {'type': 'ineq', 'fun': lambda x: max_bank - x[0]},
            {'type': 'ineq', 'fun': lambda x: x[0] - min_bank}
            
    )
    val = minimize(path_cost,guess, method='SLSQP', tol=1E-10, constraints=cons)
    return val.x

def funs(t, ba):
    tr = 111.6**2 / (11.26*math.tan(ba[0]))
    # print(tr)
    h = tr + 750
    k = -1282
    nodes = [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()]
    Bx = nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By =  nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    circle_eq = (Bx-h)**2+(By-k)**2 - tr**2
    return circle_eq

def find_diff_entry(ba):
    print(np.rad2deg(ba[0]))
    pi = np.pi
    # print(type(ba), ba[0])
    t = 0.25
    mindiff = 100
    nodes = [np.array([750, 1461.58471522, 1475]).flatten(),np.array([-50, -50, 450]).flatten()]
    path = manual_bez(P0 = [nodes[0][0], nodes[1][0]],
                      P1 = [nodes[0][1], nodes[1][1]],
                      P2 = [nodes[0][2], nodes[1][2]], 
                      points = 200)
    Bx = lambda t: nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = lambda t: nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    tr = 111.6**2 / (11.26*math.tan(ba[0]))
    print('TR', tr)
    h = tr + 750
    k = -1282
    circle_eq = lambda t: (Bx(t)-h)**2+(By(t)-k)**2 - tr**2
    # eqs = [Bx, By, circle_eq]
    S = fsolve(circle_eq, 0.25)
    # S = fsolve(funs, (0.25, ba))
    # t_values = np.linspace(0, 1, 500)
    # circle_eq_values = [circle_eq(t) for t in t_values]
    # plt.plot(path[0], path[1])
    # plt.plot(t_values, circle_eq_values)
    # plt.axhline(0, color='red', linestyle='--')
    # plt.xlabel('t')
    # plt.ylabel('circle_eq(t)')
    # plt.title('Circle Equation vs. t')
    # plt.show()
    print(Bx(0.601342376572023))
    print(S)
    print(Bx(S), By(S))
    # if isinstance(S, sp.Equality):
    # # Handle single equality case
    #     solutions = [S]
    #     print(type(solutions),solutions[0])
    # print(t)
    # S = np.array(S)
    for i in S:
        # print(i)
        if i>=0 and i <= 1:
            # print(i[0], i[0]*200, np.round(float(i[0]*200)))
            index = int(np.round(float(i*200)))
            print(index)
            if index == 200:
                index = 199
            # print(i[0]*200, index)
            # intersect_dist
            bez_angle = np.arctan2(By(i)-By(i-0.01), Bx(i)-Bx(i-0.01))
            print(np.rad2deg(bez_angle))
            x_l = [i for i in np.linspace(750, Bx(i), 200)] #Space between nominal path and bez
            y = [k+np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn
            # plt.plot(path[0], path[1])
            # plt.plot(x_l, y)
            # plt.axis('equal')
            # plt.show()
            print('TEXT', y[-1], path[1][index])
            print(x_l[-1], path[0][index])
            # if y[index]==path[1][index] and x_l[index]==path[0][index]:
            int_angle = np.arctan2(y[-1]-y[198], x_l[-1] - x_l[198])
            print(int_angle)
            diff = np.abs(bez_angle-int_angle)
            if diff < mindiff and np.abs(y[-1] - By(i))<=.00001:
                mindiff = diff
                plt.plot(path[0], path[1])
                plt.plot(x_l, y)
                plt.scatter(path[0][index], path[1][index], color = 'yellow', marker = '*', s = 100)
                plt.scatter(x_l[-1], y[-1], color = 'purple', marker = '*', s = 100)
                plt.axis('equal')
                ba2 = np.rad2deg(ba)
                # plt.title(f'{ba2}')
                # plt.show()
            
            # else:
            #     diff = 'NA'
            #     ind = i
                # print(mindiff)

    return mindiff

def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]


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
        # print(t)
        # S = np.array(S)
        for i in S:
            if i[0].is_real:
                # print(i[0])
                sols.append([ba, i[0]])
                
    # print(sols)
    return sols
        


        
if __name__=='__main__':
    start = time.time()
    wpts, x_wpts, y_wpts, bez1, bez2 = MBO.toCallOutside(188, np.deg2rad(20), 1.5*3.71, 1.5*3.71, np.deg2rad(90), 
                                        [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()],
                                        [np.array([1475, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()],
                                        1000, 0, 850, 1)
    
    # Sols = entryTurn(111.6, bez1, [750, -1282], [np.array([750, 1461.58471522, 1475]).flatten(),np.array([-50, -50, 450]).flatten()])
    # print(Sols)
    # end = time.time()
    # print(end-start)
    # print(np.deg2rad(29.22))
    optim_sol = solve_optimEntry(np.deg2rad(25), np.deg2rad(30), np.deg2rad(25))
    print(np.rad2deg(optim_sol[0]))