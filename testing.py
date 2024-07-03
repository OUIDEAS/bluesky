import numpy as np
from pyproj import Proj
import matplotlib.pyplot as plt

def Meters_To_WSG84(waypoints, home):
        # convert position back to LAT/LONg
        p = Proj(proj='utm',zone=17, ellps='WGS84', preserve_units=False)
        homeX, homeY = p(home[1], home[0])
        waypoints = np.array(waypoints)
        asize = waypoints.shape
        if (len(asize) > 1):
            waypoints_LongLat = []
            for pt in waypoints:
                x = (pt[0] + homeX)
                y = (pt[1] + homeY)
                lon, lat = p(x,y,inverse=True)
                altitude = 150
                if(len(pt)>2):
                    altitude = pt[2]
                waypoints_LongLat.append([lat, lon, altitude])
            return waypoints_LongLat

        else:
            x = (waypoints[0] + homeX)
            y = (waypoints[1] + homeY)
            lon, lat = p(x,y,inverse=True)
            altitude = 0
            if(len(waypoints)>2):
                altitude = waypoints[2]
            return [lat, lon, altitude]


# lat-lon 39.42, -83.2 lon = x

# home = [39.42, -83.2]
# border = [0, 2200]
# border2 = Meters_To_WSG84(border, home)
# print(border2)
# border2[1] = home[1]
# plt.scatter(home[1], home[0])
# plt.scatter(border2[1], border2[0])
# plt.show()

nodes = np.zeros((4, 2), dtype = object)
print(nodes)
nodes[0][0]= np.array([1, 2, 23, 4, 53])
print(nodes[0][0], nodes[0][1])
nodes = np.array([[0, 0], [0,0]], dtype = object)
print(nodes)
nodes[0][0]= np.array([1, 2, 23, 4, 53])
print(nodes[0][0], nodes[0][1])
