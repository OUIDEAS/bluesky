import numpy as np
from matplotlib import pyplot as plt    
from matplotlib.patches import Circle
from scipy.optimize import minimize
import math
import time
import scipy
import sympy as sp
from scipy.optimize import brentq
import MichaelBezierOptimize as MBO

# def solve_optimEntry(guess, max_bank):
#     def path_cost(ba):

#         return
#     cons = (
#             {'type': 'ineq', 'fun': lambda x: max_bank - x},
            
#     )

#     return

# def circlePath(ba, )





def entryTurn(velocity, path, start, nodes):
    pi = np.pi
    t = sp.symbols('t')
    Bx = nodes[0][1] + (nodes[0][0] - nodes[0][1]) * (1 - t)**2 + (nodes[0][2] - nodes[0][1]) * t**2
    By = nodes[1][1] + (nodes[1][0] - nodes[1][1]) * (1 - t)**2 + (nodes[1][2] - nodes[1][1]) * t**2
    sols = []
    
    # for t_index in np.linspace(0, 200, 1):
    for ba in np.linspace(25, 30, 10):
        '''Circle Parameters'''
        tr = velocity**2 / (11.26*np.tan(np.deg2rad(ba)))
        h = tr + start[0]
        k = start[1]
        circle_eq = (Bx-h)**2+(By-k)**2 - tr**2
        eqs = [circle_eq, t>0, t<1]
        S = sp.solve(eqs,t, dict = True, real = True)
        print(S[t])
        # S = np.array(S)
        for i in S:
            sols.append(S[t])
    print(sols)
    return S
        


        




    return
if __name__=='__main__':
    start = time.time()
    wpts, x_wpts, y_wpts, bez1, bez2 = MBO.toCallOutside(188, np.deg2rad(20), 1.5*3.71, 1.5*3.71, np.deg2rad(90), 
                                        [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()],
                                        [np.array([1475, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()],
                                        1000, 0, 850, 1)
    
    S = entryTurn(111.6, bez1, [750, -1282], [np.array([750, 1400, 1475]).flatten(),np.array([-50, 0, 450]).flatten()])
    end = time.time()
    print(end-start)