"""
"""

from dolfin import *
import scipy.io
import numpy as np

import os

def mult(trafo,coeffs):
  c = []
  for i in range(len(points)):
    y = sum([float(trafo[i,j]*coeffs[j]) for j in range(len(points))])
    c.append(y)
  return c

#parameters['reorder_dofs_serial'] = False


# Create mesh and define function space
n = 5

mesh = UnitSquareMesh(n, n)
#plot(mesh)
V = FunctionSpace(mesh, 'DG', 2)
bigMesh = refine(refine(mesh))
bigV = FunctionSpace(bigMesh, 'DG', 2)

k = 1
print float(k)/n
edgeMesh = RectangleMesh(1-3*float(3*k)/n, 1-float(2*k)/n, 1, 1, 3*k, k)
edgeV = FunctionSpace(edgeMesh, 'DG', 2)
edgePlotV = FunctionSpace(refine(refine(edgeMesh)), 'DG', 2)

# Define boundary conditions
#u0 = Constant(0.0) #const rhs
#u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
#u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1

#rhs
#f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
f = Constant(7.0)#simpleMongeAmpere
#f = Constant(1.0) #simpleMongeAmpere2
#f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

#define exact solution
u_e = interpolate(u0, V)


#====================================
#define laplace's iteration
#====================================


#start solution
u_ = Function(V)      # the most recently computed solution
#u_ = u_e

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)

#start solution via laplace
#penalty
sigma = 200*7

w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
#coeff = as_matrix([[1,0],[0,1]])
coeff = as_matrix([[4,-3],[-3,4]])

#define bilinear form
a = inner(coeff*nabla_grad(v), nabla_grad(u))*dx \
  - inner(jump(v,n),avg(coeff*nabla_grad(u)))*dS\
  - inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS\
  + Constant(sigma)('+')/h('+')* jump(u)*jump(v)*dS \
  - v*inner(n,coeff*nabla_grad(u))*ds \
  - u*inner(n,coeff*nabla_grad(v))*ds \
  + Constant(sigma)/h*v*u*ds
#  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds
#  + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \

#define rhs functional
L = inner(-Constant(2.0)*f,v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds

solve(a == L, u_)


def h(x,y):
  return max([.5*x, y-0.5, 2*x+y-1, -5*x+y-4, 2-10*x*x+10*y*y])


#interactive()
#exit()

#create array
u_array = u_.vector().array()
#print "u_array before ", u_array


#====================================
#convexify
#====================================


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


#print g.vector().array()

def middlePoint(a,b):
  return (a+b)/2

coord = mesh.coordinates()


#file to write points into
tmpfile = open('points_to_convexify.dat', "w")

dofmap = V.dofmap()

suspicious_begin = 0
suspicious_end = 3

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
  coeffsLagrange = map(u_, pointsLagrange)
  
  
  #transform into bezier basis  
  coeffBezier = mult(trafo,coeffs)
  
  coeffsdoppelt = mult(trafo_inv, coeffBezier)

  edge_ends = [(0,3), (2,3), (2,4), (1,4), (1,5), (0,5), (3,5), (3,4), (4,5)] 

  middle_pts = [middlePoint(points[i],points[j]) for i,j in edge_ends]
  middle_val_pts = [middlePoint(coeffBezier[i],coeffBezier[j]) for i,j in edge_ends]
  
  #points_with_middle_pts[len(points_with_middle_pts):] = points

  for j in dofmap.cell_dofs(cell_ind):
    if 569 == j:
      print "entry 569"
      print "cell index ", cell_ind
      print "coeffsLagrange ", coeffsLagrange
      print "coeffsBezier ", coeffBezier
      print "coeffs ^2 ", coeffsdoppelt
      print "coeffs ", coeffs

  
  if cell_ind >= suspicious_begin and cell_ind <= suspicious_end:
    print "cell index ", cell_ind
    print "coeffsLagrange ", coeffsLagrange
    print "coeffsBezier ", coeffBezier
    print "coeffs ^2 ", coeffsdoppelt
    print "coeffs ", coeffs
    print "points ", points
    print "middle_pts ", middle_pts
    print "middle_val_pts ", middle_val_pts
    
  #write into file
  for i in range(len(points)):
      tmpfile.write("%.10e %.10e %.10e \n" % (points[i][0],points[i][1],coeffBezier[i]))

  for i in range(len(middle_pts)):
    tmpfile.write("%.10e %.10e %.10e \n" % (middle_pts[i][0],middle_pts[i][1],middle_val_pts[i]))
    

#    tmpfile.write(str(points[i][0]) + " " + str(points[i][1]) + " " + str(coeffBezier[i]) + '\n')

#close file after writing points  
tmpfile.close()

#write xml output
File('mesh.xml') << mesh
File('func_to_convexify.xml') << u_
#exit()
#convexify with c++ code
os.system("../bin/Convexify points_to_convexify.dat convexified_points.dat")

#open file with convexified coefficients
tmpfile = open('convexified_points.dat', 'r')

#create new convex function
u_convex = Function(V)

coeff_cplusplus = np.empty(len(u_array))
#coeff_cplusplus[len(coeff_cplusplus):] = np.array(u_array, copy=True)
#coeff_cplusplus= []
u_convex_array = u_convex.vector().array()

print "len u", len(u_array)
print "cplus plus len ", len(coeff_cplusplus)


#update coefficients
i = 0
new_coord = 0
ignore = False

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

i=0
#transform from bezier to lagrange basis
for cell in cells(mesh):
  cell_ind = cell.index()
  coeffs = coeff_cplusplus[i*6: i*6+6]
  print "22coeffs", coeffs
  coeffs = mult(trafo_inv,coeffs)
  for j in range(6):
    u_convex_array[dofmap.cell_dofs(cell_ind)[j]] = coeffs[j]
  if cell_ind >= suspicious_begin and cell_ind <= suspicious_end:
    print "cell index ", cell_ind
    print "coeffsBezier ", coeff_cplusplus[i*6: i*6+6]
    print "coeffsLagrange ", coeffs
  i = i+1

 
#print "u_convex after ", u_convex_array
#print "u_array after ", u_array

print "max difference ", max(u_array-u_convex_array), " at ", (u_array-u_convex_array).argmax()

u_convex.vector()[:] = u_convex_array
#u_.vector()[:] = u_array

#plot results
plot(project(u_, bigV), title = 'solution')
plot(project(u_convex,bigV), title = 'convexified solution')
#plot(project(u_-u_convex, bigV), title = 'difference convex and not convex')
plot(u_-u_convex, title = 'difference convex and not convex')

#plot(project(u_-u_convex, edgePlotV), title = 'difference convex and not convex on the edge')
#plot(project(u_, edgePlotV), title = 'solution on the edge')
#plot(project(u_convex, edgePlotV), title = 'convex. sol on the edge')

#plot(project(abs(u_-u_e),bigV), title = 'error in non convexified solution')
#plot(project(abs(u_convex-u_e), bigV), title = 'error in conv. sol.')

#plot Determinant of Hessian
test_convex = det(grad(grad(u_)))
test_convex2 = det(grad(grad(u_convex)))

bigTest_convex = project(test_convex, bigV)
bigTest_convex2 = project(test_convex2, bigV)
plot(test_convex, title = 'determinant of hessian')
plot(test_convex2, title = 'determinant of convexified hessian')

plot(project(f, bigV), title = 'right-hand side')
plot(project(abs(f-test_convex), bigV), title = 'rhs - determinant of hessian')
plot(project(abs(f-test_convex2), bigV), title = 'rhs - determinant of convexified hessian')


#hold plot
interactive()

#update start solution
u_.assign(u_convex)

exit()

#====================================
#define brenner's iteration
#====================================


#penalty
sigma = 10


#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)


u  = TrialFunction(V)
v  = TestFunction(V)

coeff = cofac (grad(grad(u)))

F  = (f-det(coeff) )*v*dx

F = F + (dot(avg(coeff)*(nabla_grad(u)('+')) , n('+')) +  dot(avg(coeff)*(nabla_grad(u)('-')), n('-')) )*avg(v)*dS

F = F + (dot(coeff*(nabla_grad(v)) , n) )*(u-u0)*ds

F = F + Constant(sigma)/h*(u-u0)*v*ds


F = action(F, u_)

J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, boundary)

problem = NonlinearVariationalProblem(F, u_, bc, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
#info(prm, True)

prm['newton_solver']['absolute_tolerance'] = 1E-6
prm['newton_solver']['relative_tolerance'] = 1E-8
prm['newton_solver']['maximum_iterations'] = 10
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['report'] = True
#prm['linear_solver'] = 'gmres'
#prm['preconditioner'] = 'ilu'
#prm['krylov_solver']['absolute_tolerance'] = 1E-9
#prm['krylov_solver']['relative_tolerance'] = 1E-7
#prm['krylov_solver']['maximum_iterations'] = 1000
#prm['krylov_solver']['gmres']['restart'] = 40
#prm['krylov_solver']['preconditioner']['ilu']['fill_level'] = 0
set_log_level(PROGRESS)

solver.solve()

#examine error
u_e_array = u_e.vector().array()
u_array = u_.vector().array()
print 'Errornorm:', errornorm(u_,u_e)

 # Plot solution and mesh
plot(u_)
plot(abs(u_-u_e))
#plot(mesh)

#Hold plot
interactive()


# Dump solution to file in VTK format
s = 'MongeAmpere.pvd'
file = File(s)
file << u_
  
