import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import math
import time

#Bezier Functions
show_animation = True


def calc_bezier_path(control_points, n_points=100):
    """
    Compute bezier path (trajectory) given control points.

    :param control_points: (numpy array)
    :param n_points: (int) number of points in the trajectory
    :return: (numpy array)
    """
    traj = []
    for t in np.linspace(0, 1, n_points):
        traj.append(bezier(t, control_points))

    return np.array(traj)


def bernstein_poly(n, i, t):
    """
    Bernstein polynom.

    :param n: (int) polynom degree
    :param i: (int)
    :param t: (float)
    :return: (float)
    """
    return scipy.special.comb(n, i) * t ** i * (1 - t) ** (n - i)


def bezier(t, control_points):
    """
    Return one point on the bezier curve.

    :param t: (float) number in [0, 1]
    :param control_points: (numpy array)
    :return: (numpy array) Coordinates of the point
    """
    n = len(control_points) - 1
    return np.sum([bernstein_poly(n, i, t) * control_points[i] for i in range(n + 1)], axis=0)


def bezier_derivatives_control_points(control_points, n_derivatives):
    """
    Compute control points of the successive derivatives of a given bezier curve.

    A derivative of a bezier curve is a bezier curve.
    See https://pomax.github.io/bezierinfo/#derivatives
    for detailed explanations

    :param control_points: (numpy array)
    :param n_derivatives: (int)
    e.g., n_derivatives=2 -> compute control points for first and second derivatives
    :return: ([numpy array])
    """
    w = {0: control_points}
    for i in range(n_derivatives):
        n = len(w[i])
        w[i + 1] = np.array([(n - 1) * (w[i][j + 1] - w[i][j])
                             for j in range(n - 1)])
    return w


def curvature(dx, dy, ddx, ddy):
    """
    Compute curvature at one point given first and second derivatives.

    :param dx: (float) First derivative along x axis
    :param dy: (float)
    :param ddx: (float) Second derivative along x axis
    :param ddy: (float)
    :return: (float)
    """
    return (dx * ddy - dy * ddx) / (dx ** 2 + dy ** 2) ** (3 / 2)


def plot_arrow(x, y, yaw, length=1.0, width=0.5, fc="r", ec="k"):  # pragma: no cover
    """Plot arrow."""
    if not isinstance(x, float):
        for (ix, iy, iyaw) in zip(x, y, yaw):
            plot_arrow(ix, iy, iyaw)
    else:
        plt.arrow(x, y, length * np.cos(yaw), length * np.sin(yaw),
                  fc=fc, ec=ec, head_width=width, head_length=width)
        plt.plot(x, y)


class Aircraft:
    def __init__(self, acid, v, h, posx, posy):
        self.acid = acid
        self.v = v
        self.h = h
        self.posx = posx
        self.posy = posy



def sim_step(ac, dt, targetwpt): #posx, posy, v, h, dt, acid, targetwpt
    posxnew = ac.posx + ac.v*np.cos(ac.h)*dt
    posynew = ac.posy + ac.v*np.sin(ac.h)*dt
    # if ac.acid == 'AC0':
    #     hnew = np.arctan2(targetwpt[1]-posynew, targetwpt[0]-posxnew)
    # else:
    hnew = ac.h
    return posxnew, posynew, hnew



def cacl_next_wpt(ac, wpts, dist):
    wpts_x = [x[0] for x in wpts]
    wpts_y = [y[1] for y in wpts]
    dx = [ac.posx - idx for idx in wpts_x]
    dy = [ac.posy - idy for idy in wpts_y]
    d = np.hypot(dx, dy)
    
    target_idx = np.argmin(d)
    
    if np.min(d) == 0:
        target_idx = np.argmin(d)+1
    if abs(wpts_x[target_idx] <= ac.posx and wpts_y[target_idx] <= ac.posy):
        target_idx+=1
    if d[target_idx] > dist:
        target_idx+=1
    
    return target_idx, d


def get_valid_points(boundx, boundy, innerbx, innerby, r):
    maxbx = np.max(boundx)
    maxby = np.max(boundy)
    minbx = np.min(boundx)
    minby = np.min(boundy)
    maxinx = np.max(innerbx)
    maxiny = np.max(innerby)
    mininx = np.min(innerbx)
    mininy = np.min(innerby)
    valid_points = [[i, j] for i in range(minbx, maxbx, r) for j in range(minby, maxby, r) if not (mininx <= i <= maxinx and mininy <= j <= maxiny)]
    x_coordinates = [point[0] for point in valid_points]
    y_coordinates = [point[1] for point in valid_points]
    return valid_points, x_coordinates, y_coordinates

def find_bezier_cp(vp, sx, sy, gx, gy, bx, by, ibx, iby):
    minCurv = 9999999900
    path_wpts = []
    for i in vp:
        if i!=[gx,gy] and i!=[sx,sy] and i[0]!=sx and i[0] > np.max(ibx):
            control_points = np.array([[sx, sy], [1500, -50], [i[0], i[1]], [1500, 3050], [gx,gy]])#[1500, -50],  [1500, 3050],
            path = calc_bezier_path(control_points, n_points=20)
            # if np.isin(path.T[0], i[0]).all() and np.isin(path.T[1], i[1]).all():
                # final_path = path
            t = 1  # Number in [0, 1]
            derivatives_cp = bezier_derivatives_control_points(control_points, 2)
            dt = bezier(t, derivatives_cp[1])
            ddt = bezier(t, derivatives_cp[2])
            # Radius of curvature
            c = curvature(dt[0], dt[1], ddt[0], ddt[1])
            if np.abs(c) < minCurv:
                minCurv = c
                print(i, c)
                final_path = path
                final_cp = control_points
    for i in range(len(final_path.T[0])):
        path_wpts.append([final_path.T[0][i], final_path.T[1][i]]) #building list of waypoints
    wpts_x = [x[0] for x in path_wpts]
    wpts_y = [y[1] for y in path_wpts]
    return final_path, final_cp, path_wpts, wpts_x, wpts_y

start_x = 750
start_y = -50
goal_x = 750
goal_y = 3050

boundx = []
boundy = []
innerbx = []
innerby = []
start_time = time.time()
for i in range(750, 1500):
    boundx.append(i)
    boundy.append(-100)
for i in range(-100, 3100):
    boundx.append(750)
    boundy.append(i)
for i in range(750, 1500):
    boundx.append(i)
    boundy.append(3100)
for i in range(-100, 3100):
    boundx.append(1500)
    boundy.append(i)
for i in range(750, 1000):
    boundx.append(i)
    innerbx.append(i)
    boundy.append(0)
    innerby.append(0)
for i in range(750, 1000):
    boundx.append(i)
    innerbx.append(i)
    boundy.append(3000)
    innerby.append(3000)
for i in range(0, 3000):
    boundx.append(1000)
    innerbx.append(1000)
    boundy.append(i)
    innerby.append(i)
for i in range(0, 3000):
    innerbx.append(750)
    innerby.append(i)

figure, ax = plt.subplots()
# ax.set_xlim(700,2000)
# ax.set_ylim(-200, 3200)
# ax.set_aspect('equal')
# plt.scatter(boundx, boundy)
valid_points, xvp, yvp = get_valid_points(boundx, boundy, innerbx, innerby, 50)
path, control_points, wpts, wpts_x, wpts_y = find_bezier_cp(valid_points, start_x, start_y, goal_x, goal_y, boundx, boundy, innerbx, innerby)
# print(wpts)
# plt.scatter(xvp, yvp, color = 'green')
# plt.plot(path.T[0], path.T[1], 'r')
# plt.scatter(wpts_x, wpts_y, marker = '*', color = 'yellow')
# plt.scatter(control_points.T[0], control_points.T[1], color = 'yellow')
#acid, v, h, posx, posy)
# ax.set_aspect('equal')
target_idx = -1
AC0 = Aircraft('AC0', 220, np.deg2rad(90), 750, 0)
EAC = Aircraft('EAC', 242, np.deg2rad(90), 750, -354) #330 = 15s
TOI = np.abs(AC0.posy - EAC.posy)/(EAC.v - AC0.v)
plt.scatter(750, 2200)
dist = np.hypot(wpts[target_idx][0]-AC0.posx, wpts[target_idx][1]-AC0.posy)
for i in range(0, 3100):
    dist_prev = dist
    EAC.posx, EAC.posy, EAC.h = sim_step(EAC, .01, [0,0])
    # target_idx, d = cacl_next_wpt(AC0, wpts, dist)
    # if dist >= dist_prev:
    #     target_idx+=1
    # if target_idx >= len(wpts):
    #     target_idx = -1
    # print(target_idx, wpts)
    AC0.posx, AC0.posy, AC0.h = sim_step(AC0, .01, wpts[target_idx])
    dist = np.hypot(EAC.posx-AC0.posx, EAC.posy-AC0.posy)
    # print(dist)
    # print(EAC.posy)
    # print(AC0.posx, AC0.posy, AC0.h, wpts[target_idx], target_idx)
    # plt.scatter(AC0.posx, AC0.posy, color = 'black')
    # plt.scatter(EAC.posx, EAC.posy, color = 'red')
    # circ = plt.Circle((EAC.posx, EAC.posy), 250, color = 'b', fill = False)
    # ax.add_artist(circ)
    if dist_prev < dist:
        print('Real TOI: ', i*0.01, 'Calced TOI: ', TOI)
        break
    if np.isclose(dist, 150, atol = 5):
        print(AC0.posy, EAC.posy, dist, i*0.01)
        print('EAC DIST TO TRAVEL: ', 3100 - EAC.posy)
        break
    # plt.pause(0.01)
# plt.grid()
plt.scatter(AC0.posx, AC0.posy, color = 'black')
plt.scatter(EAC.posx, EAC.posy, color = 'red')
stop_time = time.time()
print("TOTAL TIME: ", stop_time-start_time)



plt.show()
# print(xvp)
# calc_solution_space(boundx, boundy, innerbx, innerby, 0, 0, 0, 0, 0)

#