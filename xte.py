import sympy as sp
import numpy as np
from sympy import lambdify
from scipy.optimize import minimize
from scipy import integrate
import scipy
from scipy.optimize import Bounds
# import bezierSolver
 
class xte_class():
 
    def __init__(self):
        # Line PARAMS
        self.a = 0
        self.b = 1
        self.c = 0
       
        # CIRCLE PARAMS
        self.r = 1
        self.cx = 0
        self.cy = 0
 
        # SINE PARAMS
        self.amplitude = 0
        self.ordinary_frequency = 0.0001
        self.angular_frequency = 2 * np.pi * self.ordinary_frequency
        self.phase = 0
        self.y_shift = 0
 
        # SUPER ELLIPSE PARAMS
        self.sa = 250
        self.sb = 1750/2
        self.sn = 2
        self.sm = 5
 
        # EXP PARAMS
        self.exp_mod = 1/150
 
    def line_path(self,x):
        return self.a/self.b*x + self.c/self.b
 
    def exp_path(self, x):
        try:
            return np.exp(x * self.exp_mod)
        except:
            return sp.exp(x * self.exp_mod)
 
       
    def sine_path(self, x):
        try:
            return self.amplitude* sp.sin(self.angular_frequency * x) + self.y_shift
        except:
            return self.amplitude* np.sin(self.angular_frequency * x) + self.y_shift
 
    def xte_circle(self, x, y, ax):
        xte_list = []
        for point in range(len(x)):
            xtt = x[point]
            ytt = y[point]
 
            angle_between = np.arctan2((ytt - self.cy), (xtt - self.cx))
            dist =sp.sqrt((xtt - self.cx)**2 +(ytt - self.cy)**2 ) - self.r
            ax.plot([self.r*np.cos(angle_between)+self.cx, xtt],
                    [self.r*np.sin(angle_between)+self.cy, ytt],
                    color='black')
            xte_list.append(dist)
        return xte_list
 
    def super_ellipse(self, theta):
        try:
            return  [self.sa * np.abs((np.cos(theta)))**(2/self.sn) * np.sign(np.cos(theta)),  self.sb   * np.abs((np.sin(theta)))**(2/self.sm) * np.sign(np.sin(theta)) + self.sb]
        except:
             return [self.sa * np.abs((sp.cos(theta)))**(2/self.sn) * sp.sign(sp.cos(theta)), self.sb   * np.abs((sp.sin(theta)))**(2/self.sm) * sp.sign(sp.sin(theta))  + self.sb]
 
    def distance_superellipse(self, X, point):
        x_func, y_func = self.super_ellipse(X)
        return ( (point[0] - x_func)**2 + (point[1] - y_func)**2) **(1/2)
 
    def solve_nearest_superellipse(self, point):
        if point[0] < 0:
            init_guess = 3.14
        else:
            init_guess = 0
        results = scipy.optimize.minimize(self.distance_superellipse,[init_guess], method='Powell', args=(point),tol=0.00000000001)
        return results
 
    def distance(self, X,xtt,function,ytt):
        return sp.sqrt((X-xtt)**2 + (function-ytt)**2 )
 
    def XTE(self, x, y, symbolic_x, function, function_lambdified):
        # xte_sum = 0
        xte_list = []
        for point in range(len(x)):
            xtt = x[point]
            ytt = y[point]
            D = self.distance(symbolic_x,xtt,function,ytt)
            D_lambdify = lambdify([symbolic_x],D,modules='numpy')
            dp = minimize(D_lambdify,x0=np.array([xtt,ytt]),method='Powell')
            xdp = dp.x[0]
            ydp = function_lambdified(xdp)
            xte = self.distance(xdp, xtt, ydp, ytt)
            xte_list.append(xte)
        return(xte_list)
 
 
    def plot_dist(self, x, y, symbolic_x, function, function_lambdified, ax):
        xte_list = []
        for point in range(len(x)):
            xtt = x[point]
            ytt = y[point]
            D = self.distance(symbolic_x,xtt,function,ytt)
            D_lambdify = lambdify([symbolic_x],D,modules="numpy")
 
            dp = minimize(D_lambdify,x0=np.array([0]), method='Powell')
            xdp = dp.x[0]
           
            ydp = function_lambdified(xdp)
           
            # ax.plot([xdp,xtt],[ydp,ytt],"--",color='black')
 
            xte = self.distance(xdp, xtt, ydp, ytt)
            xte_list.append(xte)
        return xte_list
   
    def bezierXTE(self, x, y, P0, P1, P2):
        def cost(t):
            dist = (
            (P0[0] * (1-t)**2 + (1-t)*P1[0]*2*t + t**2 * (P2[0]) - x)**2 +
            (P0[1] * (1-t)**2 + (1-t)*P1[1]*2*t + t**2 * (P2[1]) - y)**2) ** (1/2)
            return dist
        bounds = Bounds([0], [1])
 
        val = minimize(cost, [1], bounds=bounds, method='Powell', tol=1E-11,)
        return val.x , cost(val.x)
 
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
if __name__ == "__main__":
    from matplotlib import pyplot as plt
    # from bezierSolver import quadraticBezier, quadraticBezierFromT
 
    p1x = -5
    p1y = -5
 
    p2x = 4
    p2y = 20
 
    p3x = 10
    p3y = 5
 
    vxL = np.linspace(-12,12,100)
    # vyL = np.linspace(-10,10,10)
    vyL = np.sin(vxL/3) * 6 + 5
 
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
 
    bez = manual_bez([p1x,p1y], [p2x, p2y],[p3x, p3y], 100)
 
    xteHandler = xte_class()
 
    ax.plot(bez[0],bez[1], c='black',label='Quadratic Bezier curve')
    ax.axis("equal")
 
    for vx in vxL:
        # for vy in vyL:
 
        vy = np.sin(vx/3) * 6 + 5
 
        plt.cla()
        ax.plot(bez[0],bez[1], c='black',label='Quadratic Bezier curve')
 
        # bez = quadraticBezier([p1x,p1y], [p2x, p2y],[p3x, p3y], 100)
 
        # xteHandler = xte_class()
       
        sol = xteHandler.bezierXTE(vx, vy, [p1x,p1y], [p2x, p2y],[p3x, p3y])
 
        nearestBez = find_bez_xy([p1x,p1y], [p2x, p2y],[p3x, p3y], sol[0])
 
        xCheck = np.linspace(vx, nearestBez[0],10)
 
        b = -(vy - nearestBez[1]) / (vx- nearestBez[0]) * vx + vy
        yCheck = (vy- nearestBez[1]) / (vx-nearestBez[0]) * xCheck + b
 
        # ax.plot(bez[0],bez[1], c='black',label='Quadratic Bezier curve')
        ax.scatter(nearestBez[0],nearestBez[1], color='black')
        ax.scatter(vx,vy, color='blue')
        ax.plot(xCheck, yCheck, '--')
        ax.plot(vxL, vyL, '--')
 
        plt.pause(0.1)
 