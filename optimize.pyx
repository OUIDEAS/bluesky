# cython: language_level=3
import numpy as np
cimport numpy as np
from scipy.optimize import minimize
import math
ctypedef np.float64_t DTYPE_t
from libc.math cimport fabs



# Helper Functions
cpdef get_line_o(np.ndarray p1, np.ndarray p2, double t):
    if t < 0:
        t = 0
    elif t >= 0.9999:
        t = 1
    x = p1[0] + (p2[0] - p1[0]) * t
    y = p1[1] + (p2[1] - p1[1]) * t
    return [x, y]


cpdef get_travel_xy(np.ndarray wps, np.ndarray point, int gi, list travel, double v):
    cdef list new_time = []
    cdef double total = 0
    for i in range(len(travel)):
        if gi > i:
            new_time.append(0)
            total += 0
        elif gi == i and i + 1 < len(wps):
            d = np.hypot(wps[i + 1][0] - point[0], wps[i + 1][1] - point[1])
            new_time.append(d / v)
            total += d / v
        elif gi == i and i + 1 >= len(wps):
            d = np.hypot(point[0] - wps[i][0], point[1] - wps[i][1])
            new_time.append(d / v)
            total += d / v
        elif gi < i:
            new_time.append(travel[i])
            total += travel[i]
    return new_time, total

cpdef total_travel_o(np.ndarray waypts, double v):
    cdef list dists = []
    cdef list travel = []
    cdef double total = 0
    for i in range(len(waypts) - 1):
        d = np.hypot(waypts[i + 1][0] - waypts[i][0], waypts[i + 1][1] - waypts[i][1])
        dists.append(d)
        travel.append(d / v)
        total += d / v
    print(f'Total STAR Route Travel Time: {total}')
    return dists, travel, total




cpdef curvature_o(np.ndarray P0, np.ndarray P1, np.ndarray P2):
    cdef double mx = np.mean([P0[0], P2[0]])
    cdef double my = np.mean([P0[1], P2[1]])
    cdef double m[2]

    # Manually assign the values to the array elements
    m[0] = mx
    m[1] = my
    cdef double center1x = np.mean([mx, P0[0]])
    cdef double center1y = np.mean([my, P0[1]])
    cdef double r1 = np.sqrt((center1x - mx) ** 2 + (center1y - my) ** 2)
    cdef double center2x = np.mean([mx, P2[0]])
    cdef double center2y = np.mean([my, P2[1]])
    cdef double r2 = np.sqrt((center2x - mx) ** 2 + (center2y - my) ** 2)
    cdef list circle1 = [center1x, center1y, r1]
    cdef list circle2 = [center2x, center2y, r2]
    cdef double p1_from_c1 = np.sqrt((circle1[0] - P1[0]) ** 2 + (circle1[1] - P1[1]) ** 2)
    cdef double p1_from_c2 = np.sqrt((circle2[0] - P1[0]) ** 2 + (circle2[1] - P1[1]) ** 2)
    cdef double area = np.abs(P0[0] * P1[1] + P1[0] * P2[1] + P2[0] * P0[1] - P0[1] * P1[0] - P1[1] * P2[0] - P2[0] * P0[0]) / 2
    cdef double kmax
    if p1_from_c1 <= circle1[2]:
        kmax = area / (np.sqrt((P0[0] - P1[0]) ** 2 + (P0[1] - P1[1]) ** 2)) ** 3
    elif p1_from_c2 <= circle2[2]:
        kmax = area / (np.sqrt((P1[0] - P2[0]) ** 2 + (P1[1] - P2[1]) ** 2)) ** 3
    else:
        kmax = (np.sqrt((m[0] - P1[0]) ** 2 + (m[1] - P1[1]) ** 2)) ** 3 / (area ** 2)
    cdef double roc = 1 / kmax
    return roc

cpdef path_length_o(np.ndarray P1, np.ndarray P0, np.ndarray P2, double t):
    cdef double ax = P0[0] - 2 * P1[0] + P2[0]
    cdef double ay = P0[1] - 2 * P1[1] + P2[1]
    cdef double bx = 2 * P1[0] - 2 * P0[0]
    cdef double by = 2 * P1[1] - 2 * P0[1]
    cdef double A = 4 * (ax ** 2 + ay ** 2)
    cdef double B = 4 * (ax * bx + ay * by)
    cdef double C = bx ** 2 + by ** 2
    cdef double b = B / (2.0 * A)
    cdef double c = C / A
    cdef double u = t + b
    cdef double k = c - (b * b)
    cdef double L = 0.5 * np.sqrt(A) * ((u * np.sqrt((u * u) + k)) -
                                        (b * np.sqrt((b * b) + k)) +
                                        (k * np.log(np.abs((u + np.sqrt((u * u) + k)) / (b + np.sqrt((b * b) + k))))))
    return L

cpdef toa_diff_o(np.ndarray P1, np.ndarray P0, np.ndarray P2, double t, double toa, double v):
    cdef double ax = P0[0] - 2 * P1[0] + P2[0]
    cdef double ay = P0[1] - 2 * P1[1] + P2[1]
    cdef double bx = 2 * P1[0] - 2 * P0[0]
    cdef double by = 2 * P1[1] - 2 * P0[1]
    cdef double A = 4 * (ax ** 2 + ay ** 2)
    cdef double B = 4 * (ax * bx + ay * by)
    cdef double C = bx ** 2 + by ** 2
    cdef double b = B / (2.0 * A)
    cdef double c = C / A
    cdef double u = t + b
    cdef double k = c - (b * b)
    cdef double L = 0.5 * np.sqrt(A) * ((u * np.sqrt((u * u) + k)) -
                                        (b * np.sqrt((b * b) + k)) +
                                        (k * np.log(np.abs((u + np.sqrt((u * u) + k)) / (b + np.sqrt((b * b) + k))))))
    cdef double total = L / v
    return np.abs(total - toa)
# Function to find Bézier curve coordinates at a given time t
def find_bez_xy_o(np.ndarray  P0, np.ndarray  P1, np.ndarray  P2, double t) -> np.ndarray :
    Bx = lambda z: P1[0] + (P0[0] - P1[0]) * (1 - z)**2 + (P2[0] - P1[0]) * z**2
    By = lambda z: P1[1] + (P0[1] - P1[1]) * (1 - z)**2 + (P2[1] - P1[1]) * z**2
    x = Bx(t)
    y = By(t)
    return np.array([x, y])

# Function to manually calculate Bézier curve segments
def manual_bez_partial_o(np.ndarray P0, np.ndarray P1, np.ndarray P2, int points, double t_start) -> np.ndarray :
    t = np.linspace(t_start, 1, points)
    x_vals = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    y_vals = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    return np.array([x_vals, y_vals])

def manual_bez_o(np.ndarray  P0, np.ndarray  P1, np.ndarray  P2, int points) -> np.ndarray :
    t = np.linspace(0, 1, points)
    x_vals = P1[0] + (P0[0] - P1[0]) * (1 - t)**2 + (P2[0] - P1[0]) * t**2
    y_vals = P1[1] + (P0[1] - P1[1]) * (1 - t)**2 + (P2[1] - P1[1]) * t**2
    return np.array([x_vals, y_vals])




def central_angle_o(np.ndarray  center, np.ndarray  point1, np.ndarray  point2) -> tuple:
    """
    Computes the central angle between two points on a circle.

    Parameters:
        center (np.ndarray): (x, y) coordinates of circle center.
        point1 (np.ndarray): First point (x, y) on the circle.
        point2 (np.ndarray): Second point (x, y) on the circle.

    Returns:
        angle_radians (float): Angle in radians.
        angle_degrees (float): Angle in degrees.
    """
    # Compute vectors from center to points
    v1 = point1 - center
    v2 = point2 - center

    # Compute central angle using the dot product formula
    dot_product = np.dot(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle_radians = np.arccos(np.clip(dot_product / norms, -1.0, 1.0))

    # Convert to degrees
    angle_degrees = np.rad2deg(angle_radians)

    return angle_radians, angle_degrees


def entryPath_o(double velocity, np.ndarray[double] ba, np.ndarray intersect, np.ndarray pos, int lr, double head) -> tuple:
    pi = np.pi
    '''Entry Into Bezier Curve'''
    
    # Ensure ba is a numpy ndarray and process it correctly
    ba = np.array(ba, dtype=np.float64)

    # Calculate turn radius (tr)
    if ba[0] > 0:
        tr = (velocity * 1.94384)**2 / (11.26 * math.tan(np.deg2rad(ba[0])))
    else:
        tr = (velocity * 1.94384)**2 / (11.26 * math.tan(np.deg2rad(ba[0])))
    tr *= 0.3048  # Convert units to meters
    
    # Calculate the entry position coordinates
    h = pos[0] + lr * tr
    k = pos[1]

    # Calculate entry path coordinates
    if lr == 1:
        y_entry = np.linspace(pos[1], intersect[1], 50)
        x_entry = h - lr * np.sqrt(tr**2 - (y_entry - k) ** 2)
    else:
        x_entry = np.linspace(pos[0], intersect[0], 50)
        y_entry = k + np.sqrt(tr**2 - (x_entry - h) ** 2)

    # Handle NaN values in the entry path
    nan_mask = np.isnan(y_entry)
    if nan_mask[0]:  # If first value is NaN, reset manually
        y_entry[0] = pos[1]
        x_entry[0] = pos[0]
    if nan_mask[-1]:  # If last value is NaN, try alternative calculation
        y_entry = np.linspace(pos[1], intersect[1], 50)
        x_entry = h - lr * np.sqrt(tr**2 - (y_entry - k) ** 2)

    # Remove NaN values (list comprehension)
    x_entry, y_entry = list(x_entry), list(y_entry)
    cleaned_data = [(x, y) for x, y in zip(x_entry, y_entry) if not math.isnan(y)]
    x_entry, y_entry = zip(*cleaned_data) if cleaned_data else ([], [])

    # Calculate the central angle and entry length
    ar, ad = central_angle_o(np.array([h, k]), np.array([x_entry[0], y_entry[0]]), np.array([x_entry[-1], y_entry[-1]]))
    entryLength = tr * ar  # Entry length
    entryTOA = entryLength / velocity  # Time of arrival at entry

    return list(x_entry), list(y_entry), ar, entryLength, entryTOA, h, k

