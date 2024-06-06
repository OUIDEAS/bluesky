import numpy as np
import sympy as sp
import scipy
import math
import matplotlib.pyplot as plt

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

def manual_bez(P0, P1, P2, points):
    t = np.linspace(0,1,points)
    return [(P1[0] + (P0[0]-P1[0])*(1 - t)**2 +(P2[0]-P1[0])*t**2),(P1[1] + (P0[1]-P1[1])*(1 - t)**2 +(P2[1]-P1[1])*t**2)]

if __name__ == '__main__':
    velocity = 111.6
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
    y_path = [i for i in range(900, 3000)]
    x = [750 for i in y_path]
    int_angles = []

    nodes2 = [np.array([1475, 1000, 750]).flatten(),np.array([450, 1000, 900]).flatten()]
    t_start = 0.35678391959798994
    # t_start = 0.5
    cp = np.array([[nodes2[0][0],nodes2[1][0]], [1500, 900], [nodes2[0][2],nodes2[1][2]]])

    optimal_bez2 = manual_bez([nodes2[0][0],nodes2[1][0]],
                                [1500, 900],
                                [nodes2[0][2],nodes2[1][2]], 200)
    
    for t in np.linspace(t_start, .9, 75):
        for ba in np.linspace(25, 45, 100):
            index = int(np.round(t*200))
            bezPos = [optimal_bez2[0][index], optimal_bez2[1][index]]
            # print(bezPos)
            bezPosPrev = [optimal_bez2[0][index-1], optimal_bez2[1][index-1]]
            bezHead = np.arctan2(bezPos[1] - bezPosPrev[1], bezPos[0] - bezPosPrev[0])
            # print(np.rad2deg(bezHead))
            tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
            h = bezPos[0] - tr*math.cos(bezHead)
            k = bezPos[1] + tr*math.sin(bezHead)
            # print(h, k)

            # x, y = sp.symbols('x, y')
            circ_int = np.sqrt(tr**2 - (750-h)**2)+k #Will be used to check if it intercepts
            # if circ_int>0:
            x_l = [i for i in np.linspace(750, bezPos[0])] #Space between nominal path and bez
            y = [k-np.sqrt(tr**2 - (x-h)**2) for x in x_l] #arc created by turn

            # y_path = [i for i in range(900, 2500)] #nominal path
            # path = [750 for i in y_path] #nominal path
            # print(x_l)
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
                int_angles.append(int_angle)
    
    print(len(exitTOA_l))
    # print(exitTOA_l)
    mindex = exitTOA_l.index(min(exitTOA_l))
    tFinal = int(np.round(t_vals[mindex]*200))
    realT = t_vals[mindex]
    centerFinal = center[mindex]
    bezPosFinal = [optimal_bez2[0][tFinal], optimal_bez2[1][tFinal]]
    trFinal = 111.6**2/(11.26*np.tan(np.deg2rad(angles[mindex])))
    x_l = [i for i in np.linspace(750, bezPosFinal[0])] #Space between nominal path and bez
    y = [centerFinal[1]-np.sqrt(trFinal**2 - (x-centerFinal[0])**2) for x in x_l] #arc created by turn

    print('T VALUE:', realT)
    print('TOA:', exitTOA_l[mindex])
    print('BANK ANGLE:', angles[mindex])
    print('INTERCEPTION ANGLE:', np.rad2deg(int_angles[mindex]))
    print('INTERCEPTION POINT:', (750,int_point[mindex]))


    plt.figure()
    plt.plot(x_l, y, linestyle = 'dashed')
    plt.plot(x, y_path, color = 'orange')
    plt.axis('equal')
    plt.plot(optimal_bez2[0], optimal_bez2[1], color = 'black')
    plt.show()
            # plt.pause(.1)

            # print(ba, t)
            # '''Find tangent'''
            # derivatives_cp = bezier_derivatives_control_points(cp, 2)
            # point = bezier(t, cp)
            # dt = bezier(t, derivatives_cp[1])
            # dt /= np.linalg.norm(dt, 2)
            # tangent = np.array([point, point + dt])
            # m = (tangent[1][1] - tangent[0][1])/(tangent[1][0]-tangent[0][0])
            # y_int = tangent[0][1] - tangent[0][0]*m

            # '''Solve Circles'''
            # tr = 111.6**2/(11.26*math.tan(np.deg2rad(ba)))
            # h, k = sp.symbols('h, k')
            # A1 = -m
            # B1 = 1
            # C1 = -y_int
            # h_tan = 750

            # case1_eq1 = (A1 * h + B1 * k + C1) / sp.sqrt(A1**2 + B1**2) - tr

            # case1_eq2 = (h - h_tan) - tr
            # case1_eq3 = (optimal_bez2[0][int(np.round(t*200))] - h)**2 + ((optimal_bez2[1][int(np.round(t*200))] - k))**2 -tr**2

            # # Solve each system separately
            # solution1 = sp.solve((case1_eq1, case1_eq2), (h, k))
            # h_var = float(solution1[h])
            # k_var = float(solution1[k])
            # sols.append([h_var, k_var])
            # circ_int = np.sqrt(tr**2 - (750-h_var)**2)+k_var
            # circ_int_prev = np.sqrt(tr**2 - (750.025-h_var)**2)+k_var

            # if circ_int>0:
            #     # print(ba, t)
            #     if int(np.round(t*200)) < 200:
            #         x_index = int(np.round(t*200))
            #     else:
            #         x_index = 199

            #     x_exit = [i for i in np.linspace(750, optimal_bez2[0][x_index], 100)]
            #     x_exit_l.append(x_exit)

            #     y_exit = [k_var-np.sqrt(tr**2 - (x-h_var)**2) for x in x_exit]
            #     y_exit_l.append(y_exit)

            #     central_angle = pi + np.arctan2(y_exit[-1]-k_var, x_exit[-1] - h_var)
            #     central_angle_l.append(central_angle)

            #     exitLength = 2*pi*tr * (central_angle/(2*pi))
            #     exitLength_l.append(exitLength)

            #     exitTOA = exitLength/velocity
            #     exitTOA_l.append(exitTOA)

            #     int_point.append(circ_int)
            #     t_vals.append(t)
            #     angles.append(ba)
            #     print('t:', t, 'bank angle:',  ba, 'exit TOA:', exitTOA)
            #     plt.figure(1)
            #     plt.plot(optimal_bez2[0], optimal_bez2[1], color = 'black')
            #     plt.plot(x, y, color = 'orange')
            #     plt.plot(x_exit, y_exit, linestyle = 'dashed')
            #     plt.axis('equal')
            #     plt.show()
                
                
    
    
    # mindex = exitTOA_l.index(min(exitTOA_l))


    # finalExitTOA = exitTOA_l[mindex]
    # finalExitLength = exitLength_l[mindex]
    # finalCentral_angle = central_angle_l[mindex]
    # final_t = t_vals[mindex]

    # if int(np.round(final_t*200)) < 200:
    #     x_index = int(np.round(final_t*200))
    # else:
    #     x_index = 199
    
    # # x_exit = [i for i in np.linspace(750, path[0][x_index], 100)]

    # # y_exit = [np.sqrt(tr**2 - (x-sols[mindex][0])**2)+sols[mindex][1] for x in x_exit]
    # x_exitf = x_exit_l[mindex]
    # y_exitf = y_exit_l[mindex]
    # print(angles[mindex], final_t, finalExitTOA)



# x_exit, y_exit, central_angle2, exitLength, exitTOA, h2, k2, bank, exit_p = exitPath(velocity, t_start2, cp, optimal_bez2)