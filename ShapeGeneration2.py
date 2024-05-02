'''
Michael Variny,  mv570819@ohio.edu
Travis Moleski,  tm769613@ohio.edu
Nathaniel Hawes, nh979416@ohio.edu

Script that is a process for generating a polygon given 2 pairs of (x,y) coordinates
and generates an A* grid and list of KOZ points to allow for more freedom in the
shape used in PRACAS collision detection

'''


# from distutils.archive_util import make_archive
from doctest import master
import math
from multiprocessing.spawn import prepare
import time
import numpy as np
from matplotlib import pyplot as plt, use
from pyparsing import line
from scipy.interpolate import CubicSpline
from shapely import geometry
from shapely.prepared import prep
# from torch import equal
import cProfile, re, pstats, io, csv
from pstats import SortKey

start_time = time.time()
res_list = []
koz_time_list = []
time_list = []

class Generation():
    def __init__(self, n, res, paint_middle):
        self.numberOfPoints = n
        self.resolution = res
        self.paint_middle = paint_middle 

    def sin_func(self, x, a, b):
        return a * np.sin(b * x)

    def lin_func(self, x, m, y1, x1):
        return (m*(x) + (y1-(m*x1)))

    def cubic_func(self, a, b, c, d, x):
        return ((a * x ** 3) +(b*x**2) +(c*x) + d)

    
    '''
    Shape2 Function: Generates points for the polygon to be generated in "plotting"
    '''
    def shape2(self):
       # master_lines = []
        x2 = -82.891
        y2 = 40.02
        x1 = -82.94
        y1 = 39.97
        mu = 0
        sigma = 0
        
        x = np.linspace(x1, x2, n)
        m = (y2 - y1) / (x2 - x1)
        CurveType = 'None'

        '''
        Creates the top curve
        '''
        if CurveType == 'linear':
            y = [self.lin_func(x[i], m , y1, x1) for i in range(0,n)]
        elif CurveType == 'cubic':
            y = [self.cubic_func(0.01, 0.1, 1, y1, x[i]) for i in range(0,n)]# + np.random.normal(mu, sigma, size=1) # Plots a cubic function
            for i in range(1, 100):
                y[0], y[-1] = y1, y2
                p = np.polyfit(x, y, 3)
                y = np.polyval(p, x)
                # print(y[0], y[-1])
            # print(y[0], y[-1])
            y[0], y[-1] = y1, y2 # endpoints end up being within 0.00000000000001 so I just set them to the values
        else:
            y = [y2 for i in range(0,n)]
            # for i in range(0,n):
            #     y.append(y2) 
            
        
        cs = CubicSpline(x, y, bc_type='natural')
        ynew = y#cs(x)

        master_lines = [[x[point], ynew[point]] for point in range(0,n)] # Creates a list of combined x and y points
        
        '''
        Creates the line using the second given set of bounds (x2,y2)
        '''
        xp2 = x2
        yp2 = np.linspace(y2, y1, n)
        for point in range(0,n):
            master_lines.append([xp2, yp2[point]])

        '''
        Creates the bottom curve/baseline
        '''        
        xp1 = np.linspace(x1, x2, n)
        yp1 = [y1 for i in range(0,n)]

        # yp1[0] = y1
        # yp1[-1] = y2

        # yp1, xp1 = np.flip(yp1), np.flip(xp1)
        for point in range(0,n):
            master_lines.append([xp1[point], yp1[point]])
        
        '''
        Creates the line using the first given set of bounds (x1,y1)
        '''
        # yp3 = np.linspace(0, y1, n)
        xp3 = x1
        yp3 = np.linspace(y2, y1, n)
        for point in range(0,n):
            master_lines.append([xp3, yp3[point]])

        return master_lines

    def varTurnRate(self,TurnRate_right,TurnRate_left, LookAheadTime,V):
        self.V = V
        self.tau = LookAheadTime
        x1 = 0                          # Initial X and Y Coordinates
        y1 = 0
        A = V * LookAheadTime           # Arc Length

        if TurnRate_left != 0:
            rL     = (V*(360/TurnRate_left))/(2*np.pi)   # Left Turning Radius
            AnormL = A / rL                              # Normalized Arc with the left radius
            x2L    = rL * np.sin(AnormL)                 # Find the x coordinate after max left turn
            y2L    = rL * (1 - np.cos(AnormL))           # Find the y coordinate after max left turn
        else:
            x2L = V * LookAheadTime
            y2L = 0

        if TurnRate_right != 0:
            rR = (V*(360/TurnRate_right))/(2*np.pi)      # Right turning radius
            AnormR = A / rR                              # Normalized Arc with the right radius
            x2R = rR * np.sin(AnormR)                    # Find the x coordinate after max right turn
            y2R = -rR * (1 - np.cos(AnormR))             # Find the y coordinate after max right turn (Y coordinate flipped because initial heading = 0)
        else:
            x2R = V * LookAheadTime
            y2R = 0

        xL = np.linspace(x1,x2L,self.numberOfPoints)     # Create x and y points along line from initial x to final x (after max turns)
        yL = np.linspace(y1,y2L,self.numberOfPoints)     # ----> All these x and y points are used primarily for the polygon creation

        xR = np.linspace(x2R,x1,self.numberOfPoints)
        yR = np.linspace(y2R,y1,self.numberOfPoints)

        # Begin to compile all points/calculations
        master_line = []

        # Creates the top line (left turning boundary)--------------------------------------------------------------------------------
        master_line = [[xL[point], yL[point]] for point in range(0,self.numberOfPoints)] # Creates a list of combined x and y points
        
        # Creates the line at the right end of the graph------------------------------------------------------------------------------
        # Create the top line on the right (0 < theta < ThetaMaxLeft)
        if TurnRate_left != 0:
            xSL = np.zeros(n)
            ySL = np.zeros(n)

            ThRangeL = np.linspace(TurnRate_left,0,n, endpoint = False)
            for i in range(n):
                r = (V*(360/ThRangeL[i]))/(2*np.pi)
                Anorm = A/r
                xSL[i] = r * np.sin(Anorm)
                ySL[i] = r * (1- np.cos(Anorm))
            
            for point in range(n):
                master_line.append([xSL[point],ySL[point]])

        # Create the bottom line on the right (ThetaMaxRight < theta < 0)
        if TurnRate_right != 0:
            xSR = np.zeros(n)
            ySR = np.zeros(n)

            xSR[0] = V * LookAheadTime
            ySR[0] = 0

            ThRangeR = np.linspace(0,TurnRate_right,n)
            for i in range(1,n):
                r = (V*(360/ThRangeR[i]))/(2*np.pi)
                Anorm = A/r
                xSR[i] = r * np.sin(Anorm)
                ySR[i] = -r * (1- np.cos(Anorm))

            for point in range(n):
                master_line.append([xSR[point],ySR[point]])

        # Creates the bottom line (right turning boundary)----------------------------------------------------------------------------
        for point in range(n):
            master_line.append([xR[point], yR[point]])
       
        # Creates the line at the left end of the graph-------------------------------------------------------------------------------
        xp3 = x1
        yp3 = np.linspace(-y1, y1, self.numberOfPoints)
        for point in range(self.numberOfPoints):
            master_line.append([xp3, yp3[point]])

        return master_line

    '''
    toPolygon:
    Uses shapely to convert the (x,y) points to shapely 'Points' and generate a polygon and a prepared polygon
    '''
    def toPolygon(self, master_lines):

        
        master_list = [geometry.Point(p[0], p[1]) for p in master_lines]
        master_list.append(master_list[0])
        master_list.append(master_list[-1]) 

        polygon = geometry.Polygon([[p.x, p.y] for p in master_list])
        

        # create prepared polygon
        prep_polygon = prep(polygon)

        return prep_polygon, polygon, master_list

    '''
    makeGrid:
    Generates all the points in the A* Grid.
    Size is dictated by the variable 'res.' This is a res x res grid.
    res = 20 -> 20x20
    '''

    def makeGrid(self, poly):
        
        latmin, lonmin, latmax, lonmax = poly.bounds
        
        points = [geometry.Point(p1, p2) for p1 in np.linspace(latmin, latmax, res) 
        for p2 in np.linspace(lonmin, lonmax, res)]

        return points

    '''
    makeKOZ: 
    Generates the keep out zone and returns a list of shapely points
    '''
    def makeKOZ(self, ppoly, points):
        # '''
        # Uses shapely to convert the (x,y) points to shapely 'Points' and generate a polygon and a prepared polygon
        # '''
        # master_list = [geometry.Point(p[0], p[1]) for p in xytot]
        # master_list.append(master_list[0])
        # master_list.append(master_list[-1]) 

        # polygon = geometry.Polygon([[p.x, p.y] for p in master_list])
        
        # latmin, lonmin, latmax, lonmax = poly.bounds
        #print(latmin, lonmin, latmax, lonmax )

        # create prepared polygon
        # prep_polygon = prep(polygon)
        
        # '''
        # Generates all the points in the A* Grid
        # '''
        # points = [geometry.Point(p1, p2) for p1 in np.linspace(latmin-2, latmax+2, res) 
        # for p2 in np.linspace(lonmin-2, lonmax+2, res)]
        
        '''
        Identifies points inside the polygon for the KOZ
        Also makes a list of points that can be compared to another
        '''
        if paint_middle == True:
            keep_out_points = [i for i in filter(ppoly.contains, points)]

        #list_ko_points = [(i.x, i.y) for i in keep_out_points]

        '''
        Creates a square using shapely that moves along the grid 
        The square moves up a column then starts at the bottom of the next column
        The square is then compared to the polygon, if the two overlap at all then the 4 points
        are added to the KOZ.
        Variable 'skip' is used to make sure that there is no weird polygon whenever the counter gets
        to the top of the grid. It increments by 'res'.
        '''
        skip = res-1
        for i in range(0, len(points)-res-1):
            if i == skip:
                i = i + 1
                skip = skip + res
            p1 = points[i]
            p2 = points[i+1]
            p3 = points[i+res+1]
            p4 = points[i+res]
            square = geometry.Polygon([p1, p2, p3, p4]) #Creates a square of gridpoints that moves along the grid
            if ppoly.intersects(square) :#and set(square.exterior.coords) & set(list_ko_points):
                keep_out_points.append(p1)
                keep_out_points.append(p2)
                keep_out_points.append(p3)
                keep_out_points.append(p4)
        
        return keep_out_points
    def removeDupes(self, keep_out_points, points):
        # po2 = []
        # [po2.append(x) for x in points if x not in po2 and keep_out_points]
        return points
    def plotting(self, points, pKOZ, poly, title):
        KOZ_Fill      = False                 # Plotting Options
        A_Star        = True
        PlotGrid      = True
        EdgeMarkers   = True
        HeadingVector = False

        fig, ax = plt.subplots()           # Construct Figure
        fig.set_size_inches((8, 8))        # Set Figure Size

        if PlotGrid:
            ax.grid()

        if EdgeMarkers:
            for item in xytot:
                ax.scatter(item[0], item[1], None, 'g') 

        if KOZ_Fill:
            xs,ys = poly.exterior.xy
            ax.fill(xs,ys,fc='r',ec='black',label='Keep Out Zone')
        else:
            x,y = poly.exterior.xy
            ax.plot(x,y,'r', label='Keep Out Zone') 

        if HeadingVector:
            x_heading = [0, self.V*self.tau/2]
            y_heading = [0,0]
            ax.plot(x_heading,y_heading,'-k',label='Heading Vector')
            ax.plot(x_heading[0],y_heading[0],color='yellow',marker='*',markeredgecolor = 'k',markersize = 20,label='Own Uav')
            ax.plot(x_heading[1],y_heading[1],color='k',marker='>',markersize = 12)

        if A_Star:
            for p in points:
                ax.scatter(p.x, p.y, color = 'blue') 

            for p in pKOZ: #paints invalid points red
                ax.scatter(p.x,p.y,color='red')
        
        # plt.axis('equal')
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.legend(loc = 'lower left', frameon=True,edgecolor='k',framealpha=1)
        plt.title(title)


# for i in range(0, 19):
if __name__ == "__main__":
    t_start = time.time()
    res = 10
    n = res
    paint_middle = True
    kozObj = Generation(n, res, paint_middle)

    Profilers        = False    # This variable turns the profilers for all the functions on or off
    VarTurnRatePlots = False     # Produces proof of concept plots for VTR

    if Profilers == True:
        '''
        Profilers for all of the functions
        Last line is plt.show()
        '''
        profiler = cProfile.Profile()
        profiler.enable()

        t11 = time.time()
        xytot = kozObj.varTurnRate(20,10,10,7)
        # xytot = kozObj.shape2()
        t12 = time.time()

        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        result = io.StringIO()
        pstats.Stats(profiler,stream=result).print_stats()
        result=result.getvalue()
        # chop the string into a csv-like buffer
        result='ncalls'+result.split('ncalls')[-1]
        result='\n'.join([','.join(line.rstrip().split(None,5)) for line in result.split('\n')])
        # save it to disk
        with open('Shape2.csv', 'w+') as f:
            # f=open(result.rsplit('.')[0]+'.csv','w')
            f.write(result)
            f.close()

        profiler = cProfile.Profile()
        profiler.enable()
        
        t21 = time.time()
        ppoly, poly, master_list = kozObj.toPolygon(xytot)
        t22 = time.time()

        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        result = io.StringIO()
        pstats.Stats(profiler,stream=result).print_stats()
        result=result.getvalue()
        # chop the string into a csv-like buffer
        result='ncalls'+result.split('ncalls')[-1]
        result='\n'.join([','.join(line.rstrip().split(None,5)) for line in result.split('\n')])
        # save it to disk
        with open('toPolygon.csv', 'w+') as f:
            # f=open(result.rsplit('.')[0]+'.csv','w')
            f.write(result)
            f.close()
        
        profiler = cProfile.Profile()
        profiler.enable()

        t31 = time.time()
        points = kozObj.makeGrid(poly)
        t32 = time.time()

        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        result = io.StringIO()
        pstats.Stats(profiler,stream=result).print_stats()
        result=result.getvalue()
        # chop the string into a csv-like buffer
        result='ncalls'+result.split('ncalls')[-1]
        result='\n'.join([','.join(line.rstrip().split(None,5)) for line in result.split('\n')])
        # save it to disk
        with open('makeGrid.csv', 'w+') as f:
            # f=open(result.rsplit('.')[0]+'.csv','w')
            f.write(result)
            f.close()
        
        profiler = cProfile.Profile()
        profiler.enable()

        t41 = time.time()
        pKOZ = kozObj.makeKOZ(ppoly, points)
        t42 = time.time()

        profiler.disable()

        t51 = time.time()
        po2 = kozObj.removeDupes(pKOZ, points)
        t52 = time.time()

        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        result = io.StringIO()
        pstats.Stats(profiler,stream=result).print_stats()
        result=result.getvalue()
        # chop the string into a csv-like buffer
        result='ncalls'+result.split('ncalls')[-1]
        result='\n'.join([','.join(line.rstrip().split(None,5)) for line in result.split('\n')])
        # save it to disk
        with open('pKOZ.csv', 'w+') as f:
            # f=open(result.rsplit('.')[0]+'.csv','w')
            f.write(result)
            f.close()

        ttot = ((t42-t41) + (t32 - t31) + (t22-t21) + (t12-t11) + (t52-t51))
        tkoz = (t42-t41)
        tdupes = (t52-t51)

        print('TOTAL TIME OF FUNCTIONS WITHOUT PLOTTING... ', ttot)
        print('TOTAL KOZ GENERATION TIME... ',tkoz )
        print(tdupes)
        # res_list.append(res)
        # time_list.append(ttot)
        # koz_time_list.append(tkoz)
        # res = res + 10
        # plot = kozObj.plotting(points, pKOZ, poly); 
        # plt.plot(res_list, time_list, "-b", label = "Total Time")
        # plt.plot(res_list, koz_time_list, "-g", label = "KOZ only time")
        # plt.title("Resolution vs Times")
        # plt.xlabel("Resolution(number)")
        # plt.ylabel("Time(s)")
        # plt.yticks(np.linspace(koz_time_list[0], koz_time_list[18], 20))
        # plt.legend(loc = "upper left")
        plt.show()
        
    else:
        # xytot = kozObj.varTurnRate(15,5,10,7)
        xytot = kozObj.shape2()
        ppoly, poly, master_list = kozObj.toPolygon(xytot)
        points = kozObj.makeGrid(poly)
        
        pKOZ = kozObj.makeKOZ(ppoly, points)
        
        plot = kozObj.plotting(points, pKOZ, poly,'Test')
        plt.show()

    if VarTurnRatePlots:
        VTR1 = kozObj.varTurnRate(10,10,10,7)   # Constant TR
        VTR2 = kozObj.varTurnRate(20,10,10,7)   # Variable TR I
        VTR3 = kozObj.varTurnRate(5,10,10,7)    # Variable TR II
        VTR4 = kozObj.varTurnRate(10, 0,10,7)   # Lost Control of Right Turn

        ppoly1, poly1, master_list1 = kozObj.toPolygon(VTR1)
        ppoly2, poly2, master_list2 = kozObj.toPolygon(VTR2)
        ppoly3, poly3, master_list3 = kozObj.toPolygon(VTR3)
        ppoly4, poly4, master_list4 = kozObj.toPolygon(VTR4)

        points1 = kozObj.makeGrid(poly1)
        points2 = kozObj.makeGrid(poly2)
        points3 = kozObj.makeGrid(poly3)
        points4 = kozObj.makeGrid(poly4)

        pKOZ1 = kozObj.makeKOZ(ppoly1, points1)
        pKOZ2 = kozObj.makeKOZ(ppoly2, points2)
        pKOZ3 = kozObj.makeKOZ(ppoly3, points3)
        pKOZ4 = kozObj.makeKOZ(ppoly4, points4)

        plot1 = kozObj.plotting(points1, pKOZ1, poly1,'Constant Turn Rate KOZ\n(10 \u00B0/s both sides)')
        plot2 = kozObj.plotting(points2, pKOZ2, poly2,'Variable Turn Rate KOZ\n(20 \u00B0/s right, 10 \u00B0/s left)')
        plot3 = kozObj.plotting(points3, pKOZ3, poly3,'Variable Turn Rate KOZ\n(5 \u00B0/s right, 10 \u00B0/s left)')
        plot4 = kozObj.plotting(points4, pKOZ3, poly4,'Loss of Turn Capabilities\n(10 \u00B0/s right, 0 \u00B0/s left)')

        plt.show()