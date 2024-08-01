"""

Path planning with Bezier curve.

author: Atsushi Sakai(@Atsushi_twi)

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import sympy as sp
import math
from matplotlib.patches import Circle

show_animation = True


def calc_4points_bezier_path(sx, sy, syaw, ex, ey, eyaw, offset):
    """
    Compute control points and path given start and end position.

    :param sx: (float) x-coordinate of the starting point
    :param sy: (float) y-coordinate of the starting point
    :param syaw: (float) yaw angle at start
    :param ex: (float) x-coordinate of the ending point
    :param ey: (float) y-coordinate of the ending point
    :param eyaw: (float) yaw angle at the end
    :param offset: (float)
    :return: (numpy array, numpy array)
    """
    dist = np.hypot(sx - ex, sy - ey) / offset
    control_points = np.array(
        [[sx, sy],
         [sx + dist * np.cos(syaw), sy + dist * np.sin(syaw)],
         [ex - dist * np.cos(eyaw), ey - dist * np.sin(eyaw)],
         [ex, ey]])

    path = calc_bezier_path(control_points, n_points=100)

    return path, control_points


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


def main():
    """Plot an example bezier curve."""
    start_x = 1475.  # [m]
    start_y = 450.  # [m]
    start_yaw = np.radians(90.0)  # [rad]

    end_x = 750.  # [m]
    end_y = 900. # [m]
    end_yaw = np.radians(180.)  # [rad]
    # offset = 3.0

    # path, control_points = calc_4points_bezier_path(
    #     start_x, start_y, start_yaw, end_x, end_y, end_yaw, offset)

    # Note: alternatively, instead of specifying start and end position
    # you can directly define n control points and compute the path:
    control_points = np.array([[start_x, start_y],[1500., 900.],[end_x, end_y]])
    path = calc_bezier_path(control_points, n_points=100)
    # Display the tangent, normal and radius of cruvature at a given point
    # t = 0.86  # Number in [0, 1]
    int_point_y = []
    int_point_x = []
    angle = []
    t_vals = []
    entry_angle = []
    for t in np.linspace(0.35678391959798994, 1, 75):
        index = int(np.round(t*100))
        # print(index)
        # print(path[index])
        x_target, y_target = bezier(t, control_points)
        derivatives_cp = bezier_derivatives_control_points(control_points, 2)
        point = bezier(t, control_points)
        dt = bezier(t, derivatives_cp[1])
        ddt = bezier(t, derivatives_cp[2])
        # Radius of curvature
        radius = 1 / curvature(dt[0], dt[1], ddt[0], ddt[1])
        # Normalize derivative
        dt /= np.linalg.norm(dt, 2)
        tangent = np.array([point, point + dt])
        m = (tangent[1][1] - tangent[0][1])/(tangent[1][0]-tangent[0][0])
        y_int = tangent[0][1] - tangent[0][0]*m
        # print('y =',m,'x+',y_int)
        # print(y_int)
        # int_point = (750, m*750 + y_int, t)
        # print(int_point)
        

        for ba in np.linspace(25, 30, 5):
            print(ba, t)
            tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
            h, k = sp.symbols('h, k')
            A1 = -m
            B1 = 1
            C1 = -y_int
            h_tan = 750

            case1_eq1 = (A1 * h + B1 * k + C1) / sp.sqrt(A1**2 + B1**2) - tr

            case1_eq2 = (h - h_tan) - tr

            # Solve each system separately

            solution1 = sp.solve((case1_eq1, case1_eq2), (h, k))

            circ_int = np.sqrt(tr**2 - (750-float(solution1[h]))**2)+float(solution1[k])
            circ_int_prev = np.sqrt(tr**2 - (750.025-float(solution1[h]))**2)+float(solution1[k])
            if circ_int>0:
                # print(circ_int)
                int_point_y.append(circ_int)
                angle.append(ba)
                t_vals.append(t)
                # print(np.rad2deg(np.arctan2(float(circ_int_prev)-float(circ_int), 750.025-750)))
                entry_angle.append(np.rad2deg(np.arctan2(float(circ_int_prev)-float(circ_int), 750.025-750)))
                # print(entry_angle)
            # fig, ax = plt.subplots()
            # ax.plot(path.T[0], path.T[1], label="Bezier Path")
            # ax.plot(control_points.T[0], control_points.T[1],
            #         '--o', label="Control Points")
            # ax.plot(x_target, y_target)
            # ax.plot(tangent[:, 0], tangent[:, 1], label="Tangent")
            # # ax.plot(normal[:, 0], normal[:, 1], label="Normal")
            # ax.add_patch(Circle((solution1[h], solution1[k]), tr, color='green', fill=False))
            # plt.ylim(450, 2000)
            # plt.xlim(450, 2000)
            # circ_int = np.sqrt(tr**2 - (750-float(solution1[h]))**2)+float(solution1[k])
            # print(circ_int, tr**2 - (750-float(solution1[h]))**2)
            # ax.scatter(750, circ_int, s = 500, marker = '*', color = 'yellow')
            # plot_arrow(start_x, start_y, start_yaw)
            # plot_arrow(end_x, end_y, end_yaw)
            # ax.legend()
            # ax.axis("equal")
            # ax.grid(True)
            # plt.show()
        # normal = np.array([point, point + [- dt[1], dt[0]]])
        # curvature_center = point + np.array([- dt[1], dt[0]]) * radius
        # circle = plt.Circle(tuple(curvature_center), radius,
        #                     color=(0, 0.8, 0.8), fill=False, linewidth=1)

    # assert path.T[0][0] == start_x, "path is invalid"
    # assert path.T[1][0] == start_y, "path is invalid"
    # assert path.T[0][-1] == end_x, "path is invalid"
    # assert path.T[1][-1] == end_y, "path is invalid"
    # print(path.T[0])
    mindex = int_point_y.index(min(int_point_y))
    print(min(int_point_y), int_point_y.index(min(int_point_y)))
    print(angle[mindex])
    print(t_vals[mindex])
    print(entry_angle[mindex])
    if show_animation:  # pragma: no cover
            fig, ax = plt.subplots()
            ax.plot(path.T[0], path.T[1], label="Bezier Path")
            ax.plot(control_points.T[0], control_points.T[1],
                    '--o', label="Control Points")
            ax.plot(x_target, y_target)
            ax.plot(tangent[:, 0], tangent[:, 1], label="Tangent")
            # ax.plot(normal[:, 0], normal[:, 1], label="Normal")
            # ax.add_patch(Circle((center1x, center1y), r1, color='black', fill=False))
            # ax.add_artist(circle2)
            plot_arrow(start_x, start_y, start_yaw)
            plot_arrow(end_x, end_y, end_yaw)
            ax.legend()
            ax.axis("equal")
            ax.grid(True)
            plt.show()


# def main2():
#     """Show the effect of the offset."""
#     start_x = 10.0  # [m]
#     start_y = 1.0  # [m]
#     start_yaw = np.radians(180.0)  # [rad]

#     end_x = -0.0  # [m]
#     end_y = -3.0  # [m]
#     end_yaw = np.radians(-45.0)  # [rad]

#     for offset in np.arange(1.0, 5.0, 1.0):
#         path, control_points = calc_4points_bezier_path(
#             start_x, start_y, start_yaw, end_x, end_y, end_yaw, offset)
#         assert path.T[0][0] == start_x, "path is invalid"
#         assert path.T[1][0] == start_y, "path is invalid"
#         assert path.T[0][-1] == end_x, "path is invalid"
#         assert path.T[1][-1] == end_y, "path is invalid"

#         if show_animation:  # pragma: no cover
#             plt.plot(path.T[0], path.T[1], label="Offset=" + str(offset))

#     if show_animation:  # pragma: no cover
#         plot_arrow(start_x, start_y, start_yaw)
#         plot_arrow(end_x, end_y, end_yaw)
#         plt.legend()
#         plt.axis("equal")
#         plt.grid(True)
#         plt.show()


if __name__ == '__main__':
    main()
    # main2()