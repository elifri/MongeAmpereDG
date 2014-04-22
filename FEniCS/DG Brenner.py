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
mesh = UnitSquareMesh(20, 20)
#mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)

# Define boundary conditions
#u0 = Constant(0.0) #const rhs
#u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
#u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1

#rhs
#f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
#f = Constant(4.0)#simpleMongeAmpere
f = Constant(1.0) #simpleMongeAmpere2
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
sigma = 20

w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
coeff = as_matrix([[1,0],[0,1]])

#define bilinear form
a = inner(as_matrix(coeff)*nabla_grad(u), nabla_grad(v))*dx \
  - inner(jump(v,n),avg(coeff*nabla_grad(u)))*dS\
  - inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS\
  + Constant(sigma)('+')/h('+')* jump(u)*jump(v)*dS \
  + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \
  - v*inner(n,as_matrix(grad(grad(w)))*nabla_grad(u))*ds \
  - u*inner(n,as_matrix(grad(grad(w)))*nabla_grad(v))*ds \
  + Constant(sigma)/h*v*u*ds
#  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds

#define rhs functional
L = inner(-Constant(2.0)*f,v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds

solve(a == L, u_)


def h(x,y):
  return max([.5*x, y-0.5, 2*x+y-1, -5*x+y-4, 2-10*x*x+10*y*y])

#u_ = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0+0.25*sin(pi*x[0]*x[1]*150)'),V)
#u_ = interpolate(Expression('sin(pi*(x[0]-x[1]))'),V)#/exp(2*(x[0]*x[0]+x[1]*x[1]))+x[0]*x[0]+x[1]*x[1]'),V)

#plot(u_, title = 'solution')
#plot(u_e, title= 'exact solution')
#plot(abs(u_-u_e), title = 'error')

#interactive()
#exit()

#create array
u_array = u_.vector().array()
print "u_array before ", u_array


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


#g = Function(V)
#g = interpolate(Constant(0.0),V)
#g.vector()[10]=1

#print g.vector().array()

def middlePoint(a,b):
  return (a+b)/2

coord = mesh.coordinates()

hullPoints = []

#file to write points into
tmpfile = open('points_to_convexify.dat', "w")

dofmap = V.dofmap()

# Iterate over cells and collect dofs
for cell in cells(mesh):
  cell_ind = cell.index()
  vert_inds = cell.entities(0)

  pointsLagrange = [coord[vert_inds[0]], coord[vert_inds[1]], coord[vert_inds[2]], \
            middlePoint(coord[vert_inds[1]], coord[vert_inds[2]]), \
            middlePoint(coord[vert_inds[2]], coord[vert_inds[0]]), \
            middlePoint(coord[vert_inds[0]], coord[vert_inds[1]])]
 
  points = [coord[vert_inds[0]], coord[vert_inds[1]], coord[vert_inds[2]], \
            middlePoint(coord[vert_inds[2]], coord[vert_inds[0]]), \
            middlePoint(coord[vert_inds[1]], coord[vert_inds[2]]), \
            middlePoint(coord[vert_inds[0]], coord[vert_inds[1]])]
 
  #calculate lagrange coefficients
  coeffs = [u_array[i] for i in dofmap.cell_dofs(cell_ind)]
  coeffsLagrange = map(u_, points)
  
#  print "coeffs ", coeffs
#  print "coeffsLagrange ", coeffsLagrange
  
  #transform into bezier basis  
  coeffBezier = mult(trafo,coeffs)

#  print "coeffsBezier ", coeffBezier
  
  coeffs = mult(trafo_inv, coeffBezier)
  
  #write into file
  for i in range(len(points)):
    tmpfile.write(str(points[i][0]) + " " + str(points[i][1]) + " " + str(coeffBezier[i]) + '\n')

#close file after writing points  
tmpfile.close()

#write xml output
File('mesh.xml') << mesh
File('func_to_convexify.xml') << u_

#convexify with c++ code
os.system("../bin/Convexify points_to_convexify.dat convexified_points.dat")

#open file with convexified coefficients
tmpfile = open('convexified_points.dat', 'r')

#create new convex function
u_convex = Function(V)

coeff_cplusplus = np.array(u_array, copy=True)
u_convex_array = u_convex.vector().array()
#update coefficients
i = 0
for line in tmpfile:
  coordinates = line.split(" ")
  coeff_cplusplus[i] = coordinates[2]
  i = i+1

i=0
#transform from bezier to lagrange basis
for cell in cells(mesh):
  cell_ind = cell.index()
  coeffs = coeff_cplusplus[i*6: i*6+6]
  coeffs = mult(trafo_inv,coeffs)
  for j in range(6):
    u_convex_array[dofmap.cell_dofs(cell_ind)[j]] = coeffs[j]
  i = i+1
 
#print "u_convex after ", u_convex_array
#print "u_array after ", u_array

#print "difference ", u_array-u_convex_array

#for i in range(len(u_array)/6): #loop over all cells
#  coeffBezier = u_array[i*6: i*6+6]
#  #print i, coeffBezier
#  u_array[i*6: i*6+6] = mult(trafo_inv,coeffBezier)


u_convex.vector()[:] = u_convex_array
#u_.vector()[:] = u_array

#plot results
plot(u_convex, title = 'convexified solution')
plot(u_, title = 'solution')
plot(u_-u_convex, title = 'difference convex and not convex')
plot(u_convex-u_e, title = 'error in conv. sol.')

#hold plot
interactive()

#u_.assign(u_convex)

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

F  = (f-det(grad(grad(u))) )*v*dx

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

prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
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
  
