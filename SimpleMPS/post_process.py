import numpy as np
import matplotlib.pyplot as plt 

#Poly_file contains the coordiantes of the vertices of the polygon
#as well as the edge connectivity data showing which vertex is connected to which.
poly_file = open('./run/polygon.txt','r')
edge_poly_list = list()
for line in poly_file:
    vals = line.split()
    if len(vals)==1:
        #num_rows specifies how many rows are in the edge_list and 
        #how many rows are in the vertices list
        num_rows = int(vals[0]);
    else:
        edge_poly_list.append(vals)
poly_file.close()

edge_poly_array = np.asarray(edge_poly_list)
#poly contains the coordinates for the vertices
poly = edge_poly_array[:num_rows]
poly.astype(np.float)
#edge_list contains the edge connectivity list
#Cast as int because they are indices into poly
edge_list = edge_poly_array[num_rows:]
edge_list.astype(np.int)

#Plot the polygon first.
for row in edge_list:
    if len(row)!=1:
        point1 = poly[row[0]]
        point2 = poly[row[1]]
        xlist = np.array([point1[0],point2[0]])
        ylist = np.array([point1[1],point2[1]])
        plt.plot(xlist,ylist,'k')

#Placement_point_list contains the coordiantes fo each placement point
#One row of placement_point_list contains a placement_point
placement_point_list = list()
place_file = open('./run/points_for_placement.txt','r')
for line in place_file:
    point_vals = line.split()
    if len(point_vals)!=1:
        plt.plot(point_vals[0],point_vals[1],'rx',mfc='none',markersize=8)
        placement_point_list.append(point_vals)
place_file.close()
#Cast as float since these are floating point coordinates
placement_points = np.asarray(placement_point_list,dtype=float)

#Repeat for coverage points
cover_file = open('./run/points_for_coverage.txt','r')
for line in cover_file:
    point_vals = line.split()
    if len(point_vals)!=1:
        #'bx' indicates plot blue x's for each point
        plt.plot(point_vals[0],point_vals[1],'b+', mfc='none', markersize=12)
cover_file.close()

#We now want to plot the solution footprints. For now, we assume circles.
#The sol_points is a list containing indices into placement_points that specifies
#which placement_point is the center of a footprint.
sol_points = list()
sol_file = open('./run/solution.txt','r')
for line in sol_file:
    sol_vals = line.split()
    if len(sol_vals) !=1:
        sol_points.append(sol_vals[0])

#These are indices into placement points, so cast as int.
sol_points = np.asarray(sol_points,dtype=int)
#Index mis-match by 1, so subtract 1.
sol_values = placement_points[sol_points-1]


#Now we draw circles. We parameterize a circle with the function (xc+r*cos(t),yc+r*sin(t))
#where (xc,yc) is the center of the circle, r is the radius, and t is a parameter tha runs from
#0 to 2*pi.
t = np.linspace(0,2*np.pi,50)
radius = .25
for row in sol_values:
    xc = row[0]
    yc = row[1]
    xvals = radius*np.cos(t)
    xvals += xc
    yvals = radius*np.sin(t)
    yvals += yc
    plt.plot(xvals,yvals,'g')
    plt.plot(xc,yc,'o',mfc='none',markersize=12)
    

plt.show()     
