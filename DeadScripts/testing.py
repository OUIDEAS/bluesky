import numpy as np
from pyproj import Proj, transform
import matplotlib.pyplot as plt
import math
import utm
from utm import from_latlon as llutm
from utm import to_latlon as utmll

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

def __WSG84_To_Meters_Single(waypoint, home, projobject):
    p = projobject
    homeX, homeY = p(home[1], home[0])

    lon = waypoint[0]
    lat = waypoint[1]
    x,y = p(lat, lon)
    x = (x - homeX)
    y = (y - homeY)
        
    return [x, y]

class gnss_converter():
    def __init__(self):
        self.R = 6378137.0
 
    def merc_y(self, lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2)) * self.R
 
    def merc_x(self, lon):
        return math.radians(lon) * self.R

def get_utm_zone(longitude):
    return int((longitude + 180) / 6) + 1
# lat-lon 39.42, -83.2 lon = x

# home = [39.42, -83.2]
# border = [0,0]
# border2 = Meters_To_WSG84(border, home)
# print(border2)
# border2[1] = home[1]
# plt.scatter(home[1], home[0])
# plt.scatter(border2[1], border2[0])
# plt.show()

# nodes = np.zeros((4, 2), dtype = object)
# print(nodes)
# nodes[0][0]= np.array([1, 2, 23, 4, 53])
# print(nodes[0][0], nodes[0][1])
# nodes = np.array([[0, 0], [0,0]], dtype = object)
# print(nodes)
# nodes[0][0]= np.array([1, 2, 23, 4, 53])
# print(nodes[0][0], nodes[0][1])

# rad=gnss_converter()

# node111y = gnss_converter.merc_y(self=rad, lat = 39.4189156870712)
# node101x = gnss_converter.merc_x(self = rad, lon = -82.19730091056)
# node211y = gnss_converter.merc_y(self=rad, lat = 39.417680050070295)
# node201x = gnss_converter.merc_x(self = rad, lon = -82.19738940535501)
# print(node101x, node111y, node201x, node211y)
# m = (node111y - node211y)/(node101x-node201x)
# print(m)
# b = node211y - m*node201x
# print(b)

# node111y = utm.from_latlon(39.4189156870712, -82.19730091056)
# print(node111y[0])
# hsadf = llutm(39.4189156870712, -82.19730091056)
# print(hsadf)
# gasp = utmll(hsadf[0], hsadf[1], hsadf[2], hsadf[3])
# print(gasp[0])

# nodeslat = np.array([39, 39.1, 39.2, 39.3])
# nodeslon = np.array([-82, -82, -82, -82])
# ff = llutm(nodeslat, nodeslon)
# print(ff)
# ff = [ff[0], ff[1]]
# print(ff[0])
# p = Proj(proj='utm',zone=17, ellps='WGS84', preserve_units=False)

# homexy= [39.1977676440383, -82.2]
# pos = [39.199348426751705, -82.2]
# new = __WSG84_To_Meters_Single([pos[1],pos[0]], [homexy[1], homexy[0]], p)
# print(new)
# m = (new[1]-0)/(new[0]-0)
# print(m)
# plt.scatter(homexy[0], homexy[1])
# plt.scatter(pos[1], pos[0])
# plt.plot([pos[1], homexy[0]], [pos[0], homexy[1]])
# plt.axis('equal')
# plt.show()
# plt.scatter(new[0], new[1])
# plt.scatter(0, 0)
# plt.show()

# home = [39.4149676440383, -82.2]
# pos = [39.4165484267517, -82.2]

# # Get the UTM zone based on the longitude
# utm_zone = get_utm_zone(home[1])
# print(utm_zone)

# # Define the projection
# p = Proj(proj='utm', zone=17, ellps='WGS84')



# # Call the function
# new = __WSG84_To_Meters_Single(pos, home, p)
# print(new)

velocity  = 111.6 #m/s
turn_rate = np.deg2rad(3.5) # RAD/s
print(velocity/turn_rate)
tr = velocity/turn_rate
ba = np.arctan(111.6**2/(tr*11.6))
print(np.rad2deg(ba))
# turn_rate = np.deg2rad(4.50) # RAD/s
turn_radius = 111.6**2/(11.26*math.tan(np.deg2rad(60)))
# turn_radius*=0.3048
# print(np.rad2deg(velocity/turn_radius))
# new = Meters_To_WSG84([0, 46], [39.4, -82.2])
# print(new[0]-39.4)