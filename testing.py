import numpy as np
from pyproj import Proj, transform
import matplotlib.pyplot as plt
import math
import utm
from utm import from_latlon as llutm
from utm import to_latlon as utmll
import pandas as pd
from numpy import array
import pyclothoids as pc
from pyclothoids import Clothoid
from pyclothoids import SolveG2

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
            # x = (waypoints[0] + homeX)
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

def get_line(p1, p2, t):
    # m = (p2[1]-p1[1])/(p2[0]-p1[0])
    x = p1[0] + (p2[0]-p1[0])*t
    y = p1[1] + (p2[1]-p1[1])*t
    return [x,y]



kcmh = [39.997472, -82.891194]
teeze = [40.0865, -82.848111]
# TEEZE WAYPOINTS
melzz = [40.313417, -83.247194]
dubln = [40.202944, -83.132306]
trlgy = [40.167944, -83.061139]
polrs = [40.111194, -82.953444]
taces = [40.090111, -82.916333]

teeze_path = [melzz, dubln, trlgy, polrs, taces, teeze]
teeze_wps = LongLat_To_WSG84_Meters(teeze_path, kcmh)

# print(teeze_wps)
# plt.plot([teeze_wps[0][0],teeze_wps[1][0],teeze_wps[2][0],teeze_wps[3][0],teeze_wps[4][0],teeze_wps[5][0]],
#          [teeze_wps[0][1],teeze_wps[1][1],teeze_wps[2][1],teeze_wps[3][1],teeze_wps[4][1],teeze_wps[5][1]], linewidth = 2, color = 'blue')
# for i in range(len(teeze_wps)):
#     plt.scatter(teeze_wps[i][0], teeze_wps[i][1], color = 'orange', s=50)
# for i in range(len(teeze_wps)-1):
#     p = get_line(teeze_wps[i], teeze_wps[i+1], 0.1)
#     plt.scatter(p[0], p[1], marker = 's', s = 50)
#     p = get_line(teeze_wps[i], teeze_wps[i+1], 0.11257983)
#     plt.scatter(p[0], p[1], marker = 's', s = 50)
# plt.show()

# entry = np.zeros((1000, 20, 20), dtype = int)
# print(entry)
# entry[500][5][6] = 10
# print(entry[500][5])

# start_end = [[0,0,0] for _ in range(20)]
# start_end[5][1] = 20
# start_end[6][1] = 30
# print(np.max(start_end))
turn_radius = (57.412*1.94384)**2/(11.26*math.tan(np.deg2rad(73)))
turn_radius*=0.3048
c0 = SolveG2(229, -391, np.pi/2, 0, 450, 225/2, np.pi/2, 0)#, Dmax=1000, dmax=100)
# c0 = Clothoid.G1Hermite(229, -391, np.pi/2, 450, 225/2, np.pi/2)

# print(c0.length/57.412)
t1 = (c0[0].length + c0[1].length + c0[2].length)/57.412
print(t1)
# print(c0[1])
# p = c0.SampleXY(200)
# c1 = Clothoid.G1Hermite(450, 225/2, np.pi/2, 229, 637.5, np.pi/2)
c1 = SolveG2(450, 225/2, np.pi/2, 0, 229, 637.5, np.pi/2, 0)#, Dmax=100, dmax=100)
t2 = (c1[0].length + c1[1].length + c1[2].length)/57.412
print(t2)
print(t1+t2)
# p1 = c1.SampleXY(200)
# print(c1.length/57.412)
# print(c0.dk)
plt.figure()
plt.axis('equal')
for i in c0:
        plt.plot( *i.SampleXY(500), color = 'black' )
for i in c1:
        plt.plot( *i.SampleXY(500), color = 'black' )
# plt.plot(p[0], p[1])
# plt.plot(p1[0], p1[1])
plt.show()

# home = [39.4172, -82.2]
# wp = [0, -1078.6]
# newll = Meters_To_WSG84(wp, home)
# print(newll)
# fleet_dist = np.zeros(100)
# print(fleet_dist)
# gate = np.zeros((20, 1), dtype = str)
# # print(gate, gate[0][0])
# t_t = np.zeros((10, 2))
# print(t_t)
# t_t[0] = [8, 15]
# print(t_t[0])

# etas = [[918.9021328172605, 'AC0', 9, 'II'], [867.9018810362703, 'AC1', 0, 'II'], [883.228067292088, 'AC2', 4, 'II'], [886.523202144327, 'AC3', 6, 'II'], [889.9883177182583, 'AC4', 7, 'IV'], [876.5375112693564, 'AC5', 2, 'IV'], [893.8323160667973, 'AC6', 8, 'II'], [885.3667127197147, 'AC7', 5, 'I'], [873.777925461747, 'AC8', 1, 'IV'], [877.612124103327, 'AC9', 3, 'II']]
# gates = {'I': [], 'II': [], 'III': [], 'IV': []}

# # Group aircraft by gate dynamically
# for i in etas:
#     gates[i[3]].append(i)
# # print(gates)
# # Print aircraft by gate
# for gate, aircraft in gates.items():
#     print(f'GATE {gate} AIRCRAFT: {aircraft}')

# # Sort lists dynamically and store in dictionaries
# fcfs_sorted = {gate: sorted(aircraft, key=lambda x: x[2]) for gate, aircraft in gates.items()}
# etas_sorted = {gate: sorted(aircraft, key=lambda x: x[0]) for gate, aircraft in gates.items()}
# # print(etas_sorted)

# # for gate, aircraft in fcfs_sorted.items():
# #     c=1
# #     b=0
# #     for i in aircraft:
# #         plt.scatter(b,c, marker = 's')
# #         plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
# #         c+=1
# #     plt.text(b, c+0.5, 'FCFS')
# #     b+=1
# #     c=1
# #     for i in etas_sorted[gate]:
# #         plt.scatter(b,c, marker = 's')
# #         plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
# #         c+=1
# #     plt.text(b, c+0.5, 'Position Shift')
# #     plt.ylim((0, c+0.6))
# #     plt.xlim((-.1, 1.6))
# #     plt.show()
# # print(fcfs_sorted)
# # if fcfs_sorted == etas_sorted:
#     # print('True')
# # staff_it = {'I': [], 'II': [], 'III': [], 'IV': []}
# staff_it = {key: [vals[0]] if vals else [] for key, vals in fcfs_sorted.items()}
# # print(staff_it)
# # for key, vals in fcfs_sorted.items():
# #     # print(vals[0])
# #     if vals: 
# #         staff_it[key].append(vals[0])
# #     else:
# #         pass
# # print(fcfs_sorted['I'][0])
# # print(staff_it)
# t_it = 15
# # print(fcfs_sorted)
# for key, vals in fcfs_sorted.items():
#         for i in range(1, len(vals)):
#             if staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
#                 # print(vals[i])
#                 staff_it[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3]])
#                 # print(staff_it[key])
#             else:
#                 # print(vals[i])
#                 # print(vals[i][1], vals[i][2])
#                 staff_it[key].append([staff_it[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3]])
#                 # print(staff_it[key])
# # print(gates)
# # print(staff_it)
# t_tran = [[15, 30], [30, 15], [15, 30], [30, 15]]
# # rta = {
# #     key: [
# #         [vals[i][0]+t_tran, vals[i][2], vals[i][1], vals[i][3]]
# #         for i in range(len(vals))
# #         ]   for key, vals in fcfs_sorted.items()}
# # print(rta)
# # rta = {'I': [], 'II': [], 'III': [], 'IV': []}
# rta = []
# rta2 = []
# ind = 0
# for key, vals in staff_it.items():
#     for i in vals:
#         rs = [i[0]+t_tran[ind][0], i[0]+t_tran[ind][1]]
#         # print(t_tran[ind])
#         if ind <=1:
#             rta.append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#             rta2.append([np.max(rs), i[2],i[1], i[3], f'R{np.argmax(rs)+1}'])
#         else:
#             rta.append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+3}'])
#             rta2.append([np.max(rs), i[2],i[1], i[3], f'R{np.argmax(rs)+3}'])

#     ind+=1
#     # print(ind)

#         # if key == 'I' or key == 'II':
#         #     rs = [i[0]+t_tran[0], i[0]+t_tran[1]]
#         #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#         #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#         # if key == 'III' or key == 'IV':
#         #     rs = [i[0]+t_tran[2], i[0]+t_tran[3]]
#         #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#         #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
# cp = []
# for i in rta:
#     cp.append(i)
# cp = sorted(cp, key = lambda x: x[0])
# cp = [[v[0], i, v[2], v[3], v[3]] for i, v in enumerate(cp)]
# print(cp)
# c=0
# print(rta)
# for key, vals in staff_it.items():
    
#     for i in range(len(vals)):
#         if c<=4:
#             plt.scatter(c, fcfs_sorted[key][i][0], marker = 's')
#             plt.scatter(c+1, vals[i][0], marker = 's')
            
#             plt.text(c-0.1, fcfs_sorted[key][i][0], f'{fcfs_sorted[key][i][1]}')
#             plt.text(c+1-0.1, vals[i][0], f'{vals[i][1]}')

#             plt.plot([c, c+1], [fcfs_sorted[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
#         # plt.plot([c+1, c+2], [vals[i][0], rta[key][i][0]], linestyle = '--', color = 'k')
#         else:
#             plt.scatter(c, fcfs_sorted[key][i][0], marker = 's')
#             plt.scatter(c-1, vals[i][0], marker = 's')
            
#             plt.text(c+0.1, fcfs_sorted[key][i][0], f'{fcfs_sorted[key][i][1]}')
#             plt.text(c-1-0.1, vals[i][0], f'{vals[i][1]}')

#             plt.plot([c, c-1], [fcfs_sorted[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
#     c+=2
#     if c == 4:
#         c+=4

# c=0
# for i in cp:
#     print(f'{i[2]} ', end='')
#     plt.scatter(5, i[0], marker = 's')
#     # if 'A' in i[2]:
#         # plt.text(0.9, i[0], f'{i[2]}', weight = 'bold')
#     # else:
#     plt.text(5, i[0], f'{i[2]}', weight = 'bold')
#     c+=1
# # flat = np.array(rta).flatten()
# print('\n')
# plt.show()
# rta_d = {'R1': [], 'R2': [], 'R3': [], 'R4': []}
# for i in rta:
#     rta_d[i[4]].append(i)
# print(rta_d)

# non_rta = {'R1': [], 'R2': [], 'R3': [], 'R4': []}
# for i in rta2:
#     non_rta[i[4]].append(i)

# sta_p = {key: [vals[0]] if vals else [] for key, vals in rta_d.items()}
# non_sta = {key: [vals[0]] if vals else [] for key, vals in non_rta.items()}
# for key, vals in rta_d.items():
#     for i in range(len(vals)):
#         if sta_p[key][i-1] and vals[i][0] > sta_p[key][i][0]+t_it:
#             sta_p[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4]])
#         else:
#             sta_p[key].append([sta_p[key][i][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4]])

# print(sta_p)
# print(non_sta)
# c=0
# # for key, vals in rta_d.items():
# #         for i in vals:
# #                 print(i)
# #                 plt.scatter(c, i[0], marker = 's')
# #                 plt.text(c-0.1, i[0], f'{i[2]}')
# #         plt.text(c, i[0], f'RTA\nRunway {key}')
# #         c+=2
    
# # plt.show()
# print(rta_d)
# rta_d['R2'].pop(0)
# print(rta_d)
# # print(flat)

# run = [[888.8175377558848, 'AC12', 4, 'II', 'R2', 878.5770332638525, array([1173.43539141,  354.67517028]), 888.8175377558848, 356.91567477052934], [923.9173798556541, 'AC2', 0, 'I', 'R1', 922.0469734375755, array([124.1144099 , 952.33195789]), 923.9173798556541, 117.98481631807965], [1102.6895197165802, 'AC13', 1, 'I', 'R1', 882.9768717744743, array([ 227.71264794, 1075.70255123]), 882.9768717744743, 219.71264794210595], [1170.5675840934268, 'AC15', 2, 'I', 'R1', 893.4962438216261, array([285.07134027, 801.49813605]), 893.4962438216261, 277.07134027180075], [1227.811481910163, 'AC5', 3, 'I', 'R1', 914.4492969761197, array([312.42601917, 806.72993939]), 918.9173798556541, 308.8941020545087], [1249.8534161060657, 'AC18', 5, 'II', 'R2', 870.9120590485333, array([922.54891878, 381.13039964]), 873.8175377558848, 376.0358783501807], [1265.7434577297333, 'AC14', 6, 'I', 'R1', 912.2329546984605, array([ 358.14165272, 1162.65771837]), 913.9173798556541, 351.82607787407926], [1273.5754591091872, 'AC17', 7, 'II', 'R2', 872.538339515919, array([1277.3853397 ,  396.47872311]), 878.8175377558848, 394.75792135330244], [1303.081215809987, 'AC7', 8, 'II', 'R2', 868.8175377558848, array([1340.9433664 ,  442.26367805]), 868.8175377558848, 434.26367805410223], [1374.908633062144, 'AC16', 9, 'II', 'R2', 899.0090568066238, array([790.10313089, 483.89957626]), 899.0090568066238, 475.8995762555202], [1403.8369014534178, 'AC8', 10, 'II', 'R2', 872.7517313083106, array([1418.57988914,  516.95355725]), 883.8175377558848, 520.019363697533], [1452.0872023023658, 'AC6', 11, 'I', 'R1', 908.9173798556541, array([ 551.16982245, 1249.2197597 ]), 908.9173798556541, 543.1698224467117], [1661.092333058534, 'AC10', 12, 'IV', 'R1', 898.3288742524088, array([ 770.73926789, 1325.50522509]), 898.3409697085295, 762.7513633500045], [1780.8477972316637, 'AC9', 13, 'III', 'R2', 879.9046428861502, array([1661.7048789 ,  908.94315435]), 879.9046428861502, 900.9431543455134], [1796.2447915825546, 'AC1', 14, 'IV', 'R1', 910.5518706286291, array([ 893.69292095, 1373.87534795]), 910.5518706286291, 885.6929209539255], [1827.7977480223892, 'AC4', 15, 'III', 'R2', 896.1004886208009, array([1683.3577176,  939.6972594]), 896.1004886208009, 931.6972594015883], [1963.475258627524, 'AC3', 16, 'IV', 'R1', 888.3409697085295, array([1083.13428892, 1400.16345689]), 888.3409697085295, 1075.1342889189946], [1971.9777561496694, 'AC0', 17, 'III', 'R2', 904.3280019933156, array([1718.98188316, 1075.64975416]), 904.3280019933156, 1067.6497541563538], [2212.40863444712, 'AC19', 18, 'IV', 'R1', 923.3011692757384, array([1297.10746517, 1455.47381182]), 923.3011692757384, 1289.1074651713818], [2265.52337582773, 'AC11', 19, 'IV', 'R2', 893.1223533804006, array([1544.43979018, 1379.96378979]), 893.3409697085295, 1372.1824061192006]]

# index = [i for i, entry in enumerate(run) if entry[1] == 'AC2']
# seq_dat = run[index[0]]
# print(seq_dat)

# letters = ['A', 'B', 'C' ,'D']
# for i in letters:
#     globals()[f'fcfs{i}'] = []
# print(fcfsA)
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

velocity  = 106.68 #m/s
turn_rate = np.deg2rad(5) # RAD/s
print(velocity/turn_rate)
tr = velocity/turn_rate
ba = np.arctan(velocity**2/(tr*11.6))
# print(np.rad2deg(ba))
# # turn_rate = np.deg2rad(4.50) # RAD/s
turn_radius = 106.68**2/(11.26*math.tan(np.deg2rad(5)))
turn_radius*=0.3048
print(np.rad2deg(54.884/turn_radius))
# turn_radius*=0.3048
# print(np.rad2deg(velocity/turn_radius))
# new = Meters_To_WSG84([0, 46], [39.4, -82.2])
# print(new[0]-39.4)
# events = []
# event_l = {
#     'event': 'Simulation Start',
#     'timeStamp': 0,
#     'ACID': 'AX0, AX1, EX0',
#     'ExpNum': 0,
#     'Category': 'All Aircraft',
#     'ExpType': 'yes'
# }
# events.append(event_l)
# note_events = pd.DataFrame(events)
# print(note_events)
# note_event = {
#                             'event': 'Bezier Generation',
#                             'timeStamp': 100,
#                             'ACID': 'AX0',
#                             'ExpNum': 0,
#                             'Category': 'Fleet Aircraft',
#                             'ExpType': 'yes'
#                         }
# events.append(note_event)
# note_events = pd.DataFrame(events)
# print(note_events)
# exptype = 'SingleBez'
# path = f'{exptype}JSONsnote_events'
# print(path)


