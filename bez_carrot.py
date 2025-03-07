import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define symbols
t, x, y = sp.symbols('t x y')

# Define control points for the quadratic Bezier curve
p1x, p1y =  350446.77356838, 2732167.00004947
p2x, p2y =  378669.86220166, 2735683.48242048
p3x, p3y =  380966.05460668, 2751252.55571667
xc = p1x

r = 2000
yc = p1y + r * 1 / 2

dt = 1
v = 136
# Generate Bezier curve points for plotting
t_vals = np.linspace(0, 1, 500)
bx_vals = p1x * (1 - t_vals)**2 + 2 * (1 - t_vals) * t_vals * p2x + t_vals**2 * p3x
by_vals = p1y * (1 - t_vals)**2 + 2 * (1 - t_vals) * t_vals * p2y + t_vals**2 * p3y

# Plot the results
fig, ax = plt.subplots(figsize=(8, 8))
# Formatting the plot
ax.set_aspect('equal', adjustable='datalim')
ax.legend()
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.grid()
historyx = []
historyy = []
tS = 0
heading = np.pi/4
turnrate = np.deg2rad(3)
while yc <= p3y:
    plt.cla()
    ax.plot(bx_vals, by_vals, label="Bezier Curve", linewidth=2, color='b')
    
    # Define the Bezier curve equations
    bx = p1x * (1 - t)**2 + 2 * (1 - t) * t * p2x + t**2 * p3x
    by = p1y * (1 - t)**2 + 2 * (1 - t) * t * p2y + t**2 * p3y

    # Define the circle equation
    circle_eq = (x - xc) **2 + (y - yc) **2 - r**2

    # Define the parametric equations for the Bezier curve
    bx_eqn = bx - x
    by_eqn = by - y

    sol = fsolve((circle_eq, bx_eqn, by_eqn), (t,x,y), (tS, xc+r, yc+r))
    print(sol)
    tS = sol[0]
    sol = np.array(sol).flatten()

    hcmd_new = np.arctan2(np.float64(sol[2]-yc), np.float64(sol[1]-xc))

    rt_cmd =  np.clip((hcmd_new - heading)/dt, -turnrate, turnrate)

    heading = heading + rt_cmd * dt

    xc = dt * v * np.cos(heading) + xc
    yc = dt * v * np.sin(heading) + yc

    historyx.append(xc)
    historyy.append(yc)


    ax.scatter(sol[1],sol[2])
    ax.scatter(xc, yc)
    # plt.show()

    # Plot the circle
    circle = plt.Circle((xc, yc), r, color='r', fill=False, label="Circle")
    ax.add_artist(circle)

    ax.plot(historyx, historyy)
    plt.pause(0.1)
