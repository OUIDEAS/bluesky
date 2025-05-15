import numpy as np
from pyproj import Proj
import matplotlib.pyplot as plt
def __WSG84_To_Meters_Single(waypoint, home, projobject):
    p = projobject
    homeX, homeY = p(home[1], home[0])  # Note the order: lon, lat
    lon = waypoint[1]
    lat = waypoint[0]
    x, y = p(lon, lat)  # Note the order: lon, lat
    x = (x - homeX)
    y = (y - homeY)
    return [x, y]

def Meters_To_WSG84(waypoints, home):
        # convert position back to LAT/LONg
        p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)  
        homeX, homeY = p(home[1], home[0])
        waypoints = np.array(waypoints)
        print(len(waypoints))
        asize = waypoints.shape
        lats = []
        lons = []
        # if (len(asize) >= 1):
        waypoints_LongLat = []
        for pt in range(0, len(waypoints)):
            x = (waypoints[pt][0] + homeX)
            y = (waypoints[pt][1] + homeY)
            # print(x,y)
            lon, lat = p(x,y,inverse=True)
            altitude = 150
            if(len(waypoints[0])>2):
                altitude = waypoints[pt][2]
            waypoints_LongLat.append([lat, lon, altitude])
            lats.append(lat)
            lons.append(lon)
        return waypoints_LongLat, lats, lons

        # else:
        #     x = (waypoints[0] + homeX)
        #     y = (waypoints[1] + homeY)
        #     lon, lat = p(x,y,inverse=True)
        #     altitude = 0
        #     if(len(waypoints)>2):
        #         altitude = waypoints[2]
        #     return [lat, lon, altitude]

def LongLat_To_WSG84_Meters(waypoints, home):
    p = Proj(proj='utm',zone=17,ellps='WGS84', preserve_units=False)
    waypoints = np.array(waypoints)
    asize = waypoints.shape
    
    waypoints_meters = []
    if(len(asize) > 1): #m x 3 matrix       
        for pt in waypoints:
            if len(pt) > 2:
                newpt = __WSG84_To_Meters_Single(pt[0], home,projobject=p)
                if(asize[1]==3):#altitude
                    newpt=[newpt[0],newpt[1],pt[1], pt[2]]
            else:
                newpt = __WSG84_To_Meters_Single(pt, home, projobject=p)
                if(asize[1]==3):#altitude
                    newpt=[newpt[0],newpt[1],pt[1]]

            waypoints_meters.append(newpt)
    else:
        waypoints_meters=(__WSG84_To_Meters_Single([waypoints[0],waypoints[1]], home,projobject=p))
        altitude=0
        if(len(waypoints)==3):
            altitude = waypoints[2]
        waypoints_meters=[waypoints_meters[0],waypoints_meters[1], altitude]

    return waypoints_meters


def interp_wps(waypts, n, v):
    interps = np.zeros((len(waypts)-1, 2, n))
    dists = []
    travel = []
    # interps[0][0] = (8, 10)
    # print(interps[0][0])
    # interps[0][0][0] = 9
    # interps[0][1][0] = 10
    # print(interps[0])
    # print(interps)
    total = 0
    for i in range(0, len(waypts)-1):
        dx = waypts[i+1][0] - waypts[i][0]
        dy = waypts[i+1][1] - waypts[i][1]
        m = (waypts[i+1][1] - waypts[i][1])/(waypts[i+1][0] - waypts[i][0])
        d = np.hypot(dx, dy)
        dists.append(d)
        travel.append(d/v)
        total+= d/v
        x = np.linspace(waypts[i][0], waypts[i+1][0], n, endpoint=False)
        interps[i][0][0], interps[i][1][0] = x[0], waypts[i][1]
        for j in range(1, len(x)):
            y = waypts[i][1] + m * (x[j]-waypts[i][0])
            interps[i][0][j] = x[j]
            interps[i][1][j] = y
    print(f'Total STAR Route Travel Time: {total}')
    return interps, dists, travel, total


kcmh = [39.997472, -82.891194]
favus = [40.058917, -82.639417]
teeze = [40.0865, -82.848111]
elupy = [39.831083, -83.074778]
xavyr = [39.906167, -82.647333]


# # TEEZE WAYPOINTS
melzz = [40.313417, -83.247194]
dubln = [40.202944, -83.132306]
trlgy = [40.167944, -83.061139]
polrs = [40.111194, -82.953444]
taces = [40.090111, -82.916333]

teeze_path = [melzz, dubln, trlgy, polrs, taces, teeze]
teeze_wps = LongLat_To_WSG84_Meters(teeze_path, kcmh)
teeze_pts, teeze_dists, teeze_toa, teeze_total_t = interp_wps(teeze_wps, 100, 57.412)

# FAVUS WAYPOINTS
# bugzz = [40.565, -82.454056]
# cbuss = [40.325306, -82.5405]
# molls = [40.132139, -82.609194]
# ordiw = [39.989306, -82.666056]
# bazel = [39.981917, -82.704556]

# favus_path =[cbuss, molls, favus, ordiw, bazel]
# favus_wps = LongLat_To_WSG84_Meters(favus_path, kcmh)

# ELUPY WAYPOINTS
# jaktz = [39.591028, -83.419583]
# rscot = [39.722389,  -83.286306]
# obetz = [39.787667, -83.158389]
# edwib = [39.877472, -82.984861]
# gagbe = [39.907167, -82.927278]
# jesce = [39.903556, -82.858889]

# elupy_path = [rscot, obetz, elupy, edwib, gagbe, jesce]
# elupy_wps = LongLat_To_WSG84_Meters(elupy_path, kcmh)

# XAVYR WAYPOINTS
# scrlt = [39.502917, -82.350833]
# brtus = [39.730944, -82.473083]
# guber = [39.963222, -82.670889]
# bkeye = [39.982056, -82.706694]

# xavyr_path = [brtus, xavyr, guber, bkeye]
# xavyr_wps = LongLat_To_WSG84_Meters(xavyr_path, kcmh)


teeze_pts, teeze_dists, teeze_toa, teeze_total_t = interp_wps(teeze_wps, 100, 54.84)
print(teeze_wps)
print(teeze_dists)
teeze_dists = np.flip(teeze_dists)
new_teeze = []
for i in range(len(teeze_wps)-1, -1, -1):
    new_teeze.append(teeze_wps[i])
# print(new_teeze)
new_wp_meters = []  # Start with the first waypoint
new_wp_meters.append(new_teeze[0])

angles = []

for i in range(0, len(new_teeze)-1):
    current = new_teeze[i]
    next = new_teeze[i+1]
    xd = next[0] - current[0]
    yd = next[1] - current[1]
    a = np.arctan2(yd, xd)
    # print(a, np.pi)
    if 0>a:
        a=a+(2*np.pi)
        # print(a)
    elif 2*np.pi < a:
        a-=2*np.pi
    angles.append(a)
# print(angles)
for i in range(0, len(teeze_dists)):
    xd = 0.3* teeze_dists[i] * np.cos(angles[i])
    yd = 0.3* teeze_dists[i] * np.sin(angles[i])

    new_wp_meters.append([new_wp_meters[i][0] + xd, new_wp_meters[i][1]+ yd])
    # new_x = current[0] + 0.1*xd
    # new_y = current[1] + 0.1*yd
flip_wps = []
for i in range(len(new_wp_meters)-1, -1, -1):
    flip_wps.append(new_wp_meters[i])
    # print(f'INDICES {i, i+1}, XD: {xd}, YD: {yd}')
# print(new_wp_meters)
# for i in range(len(teeze_wps)-1, 1, -1):
#     print(i)
#     current = teeze_wps[i]
print(flip_wps)
flip_ll, l, lo = Meters_To_WSG84(flip_wps, kcmh)
for i, wp in enumerate(flip_ll):
    print(f"{wp}")
# print(new_list)
t_x, t_y = zip(*teeze_wps)
nx, ny = zip(*flip_wps)

plt.plot(t_x, t_y, color = 'blue')
plt.plot(nx, ny, color = 'orange')
plt.grid()
plt.axis('equal')
plt.show()

# # FAVUS WAYPOINTS
# bugzz = [40.565, -82.454056]
# cbuss = [40.325306, -82.5405]
# molls = [40.132139, -82.609194]
# ordiw = [39.989306, -82.666056]
# bazel = [39.981917, -82.704556]

# favus_path =[cbuss, molls, favus, ordiw, bazel]
# favus_wps = LongLat_To_WSG84_Meters(favus_path, kcmh)
# # favus_pts, favus_dists, favus_toa, favus_total_t = interp_wps(favus_wps, 100, 57.412)
# favus_pts, favus_dists, favus_toa, favus_total_t = interp_wps(favus_wps, 100, 54.84)

# # ELUPY WAYPOINTS
# jaktz = [39.591028, -83.419583]
# rscot = [39.722389,  -83.286306]
# obetz = [39.787667, -83.158389]
# edwib = [39.877472, -82.984861]
# gagbe = [39.907167, -82.927278]
# jesce = [39.903556, -82.858889]

# elupy_path = [rscot, obetz, elupy, edwib, gagbe, jesce]
# elupy_wps = LongLat_To_WSG84_Meters(elupy_path, kcmh)
# # elupy_pts, elupy_dists, elupy_toa, elupy_total_t = interp_wps(elupy_wps, 100, 57.412)
# elupy_pts, elupy_dists, elupy_toa, elupy_total_t = interp_wps(elupy_wps, 100, 54.84)

# # XAVYR WAYPOINTS
# scrlt = [39.502917, -82.350833]
# brtus = [39.730944, -82.473083]
# guber = [39.963222, -82.670889]
# bkeye = [39.982056, -82.706694]

# xavyr_path = [brtus, xavyr, guber, bkeye]
# xavyr_wps = LongLat_To_WSG84_Meters(xavyr_path, kcmh)
# # xavyr_pts, xavyr_dists, xavyr_toa, xavyr_total_t = interp_wps(xavyr_wps, 100, 57.412)
# xavyr_pts, xavyr_dists, xavyr_toa, xavyr_total_t = interp_wps(xavyr_wps, 100, 54.84)