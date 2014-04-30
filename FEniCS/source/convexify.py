
from dolfin import *
import os
import numpy as np

#====================================
#convexify
#====================================

def middlePoint(a,b):
  return (a+b)/2

def mult(trafo,coeffs):
  c = []
  for i in range(len(coeffs)):
    y = sum([float(trafo[i,j]*coeffs[j]) for j in range(len(coeffs))])
    c.append(y)
  return c

def convexify(mesh, V, u_array):
	trafo = as_matrix([ [1   ,   0,   0,0,0,0], \
						[   0,   1,   0,0,0,0], \
						[   0,   0,   1,0,0,0], \
						[-0.5,   0,-0.5,0,2,0], \
						[   0,-0.5,-0.5,2,0,0], \
						[-0.5,-0.5,   0,0,0,2]])

	trafo_inv = as_matrix([ [1   ,   0,   0,  0,  0,  0], \
							[0   ,   1,   0,  0,  0,  0], \
							[0   ,   0,   1,  0,  0,  0], \
							[0   ,0.25,0.25,  0,0.5,  0], \
							[0.25,   0,0.25,0.5,  0,  0], \
							[0.25,0.25,   0,  0,  0,0.5]])

	coord = mesh.coordinates()


	#file to write points into
	tmpfile = open('points_to_convexify.dat', "w")

	dofmap = V.dofmap()

	suspicious_begin = 0
	suspicious_end = 0

	# Iterate over cells and collect dofs
	for cell in cells(mesh):
	  cell_ind = cell.index()
	  vert_inds = cell.entities(0)

	  pointsLagrange = [coord[vert_inds[0]], coord[vert_inds[1]], coord[vert_inds[2]], \
				middlePoint(coord[vert_inds[1]], coord[vert_inds[2]]), \
				middlePoint(coord[vert_inds[2]], coord[vert_inds[0]]), \
				middlePoint(coord[vert_inds[0]], coord[vert_inds[1]])]
	 
	 #bezier control points
	  points = [coord[vert_inds[0]], coord[vert_inds[1]], coord[vert_inds[2]], \
				middlePoint(coord[vert_inds[2]], coord[vert_inds[0]]), \
				middlePoint(coord[vert_inds[1]], coord[vert_inds[2]]), \
				middlePoint(coord[vert_inds[0]], coord[vert_inds[1]])]
	 
	  #calculate lagrange coefficients
	  coeffs = [u_array[i] for i in dofmap.cell_dofs(cell_ind)]
	  #coeffsLagrange = map(u_, pointsLagrange)
	  
	  #transform into bezier basis  
	  coeffBezier = mult(trafo,coeffs)
	  
	  coeffsdoppelt = mult(trafo_inv, coeffBezier)

	  edge_ends = [(0,3), (2,3), (2,4), (1,4), (1,5), (0,5), (3,5), (3,4), (4,5)] 

	  middle_pts = [middlePoint(points[i],points[j]) for i,j in edge_ends]
	  middle_val_pts = [middlePoint(coeffBezier[i],coeffBezier[j]) for i,j in edge_ends]
	  
	  #points_with_middle_pts[len(points_with_middle_pts):] = points
	  
	  if cell_ind == 84 or cell_ind == 67:
		print "cell index ", cell_ind
#		print "coeffsLagrange ", coeffsLagrange
		print "coeffsBezier ", coeffBezier
		print "coeffs ^2 ", coeffsdoppelt
		print "coeffs ", coeffs
		print "points ", points
		print "middle_pts ", middle_pts
		print "middle_val_pts ", middle_val_pts
		
	  #write into file
	  for i in range(len(points)):
		  tmpfile.write("%.10e %.10e %.10e \n" % (points[i][0],points[i][1],coeffBezier[i]))

	  #for i in range(len(middle_pts)):
		#tmpfile.write("%.10e %.10e %.10e \n" % (middle_pts[i][0],middle_pts[i][1],middle_val_pts[i]))
	
	#close file after writing points  
	tmpfile.close()

	#convexify with c++ code
	os.system("../bin/Convexify points_to_convexify.dat convexified_points.dat")

	#open file with convexified coefficients
	tmpfile = open('convexified_points.dat', 'r')

	coeff_cplusplus = np.empty(len(u_array))
	u_convex_array = np.empty(len(u_array))

	print "len u", len(u_array)
	print "cplus plus len ", len(coeff_cplusplus)

	#update coefficients
	i = 0
	new_coord = 0
	ignore = False

	"""
	for line in tmpfile:
	  print i, " ", new_coord
	  if ignore:
		i = i+1
		if i == 9:
		  ignore = False
		  i = 0
		continue
	  coordinates = line.split(" ")
	  print "coord ", coordinates
	  coeff_cplusplus[new_coord] = coordinates[2]
	  new_coord = new_coord +1
	  i = i+1
	  if i==6:
		ignore = True
		i = 0
	  
	print "coeffs c++ ", coeff_cplusplus
"""
	for line in tmpfile:
	  coordinates = line.split(" ")
	  coeff_cplusplus[i] = coordinates[2]
	  i = i+1
	  
	#print "coeffs c++ ", coeff_cplusplus


	i=0
	#transform from bezier to lagrange basis
	for cell in cells(mesh):
	  cell_ind = cell.index()
	  coeffs = coeff_cplusplus[i*6: i*6+6]
#	  print "22coeffs", coeffs
	  coeffs = mult(trafo_inv,coeffs)
	  for j in range(6):
		u_convex_array[dofmap.cell_dofs(cell_ind)[j]] = coeffs[j]
	  if cell_ind == 84 or cell_ind == 67:
		print "cell index ", cell_ind
		print "coeffsBezier ", coeff_cplusplus[i*6: i*6+6]
		print "coeffsLagrange ", coeffs
	  i = i+1

	print "max difference ", max(u_array-u_convex_array), " at ", (u_array-u_convex_array).argmax()
	u_array = u_convex_array
	return u_convex_array

