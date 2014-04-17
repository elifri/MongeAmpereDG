"""
"""

from dolfin import *
import scipy.io
import numpy as np

import os

def mult(trafo,coeffs):
  c = []
  for i in range(len(points)):
#    print "here comes ", trafo[i,2]*coeffs[2]
#    print "here is the array ", ([ float(trafo[i,j]*coeffs[j]) for j in range(len(points))])
    y = sum([float(trafo[i,j]*coeffs[j]) for j in range(len(points))])
    c.append(y)
  return c

parameters['reorder_dofs_serial'] = False


# Create mesh and define function space
mesh = UnitSquareMesh(30, 30)
#mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)

# Define boundary conditions
#u0 = Constant(0.0) #const rhs
#u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1

#rhs
#f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
#f = Constant(4.0)#simpleMongeAmpere
#f = Constant(1.0) #simpleMongeAmpere2
f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

#define exact solution
u_e = interpolate(u0, V)


#====================================
#define laplace's iteration
#====================================


#start solution
u_ = Function(V)      # the most recently computed solution
u_ = u_e

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

plot(u_, title = 'solution')
plot(abs(u_-u_e), title = 'error')

#interactive()

#====================================
#convexify
#====================================


trafo = as_matrix([ [1   ,   0,   0,0,0,0], \
                    [   0,   1,   0,0,0,0], \
                    [   0,   0,   1,0,0,0], \
                    [-0.5,   0,-0.5,2,0,0], \
                    [   0,-0.5,-0.5,0,2,0], \
                    [-0.5,-0.5,   0,0,0,2]])

trafo_inv = as_matrix([ [1   ,   0,   0,  0,  0,  0], \
                        [0   ,   1,   0,  0,  0,  0], \
                        [0   ,   0,   1,  0,  0,  0], \
                        [0.25,   0,0.25,0.5,  0,  0], \
                        [   0,0.25,0.25,  0,0.5,  0], \
                        [0.25,0.25,   0,  0,  0,0.5]])


#g = Function(V)
#g = interpolate(Constant(0.0),V)
#g.vector()[10]=1

#print g.vector().array()

def middlePoint(a,b):
  return (a+b)/2

#vertex_to_dof_map = V.dofmap().vertex_to_dof_map(mesh)
#dof_to_vertex_map = V.dofmap().dof_to_vertex_map(mesh)

coord = mesh.coordinates()

hullPoints = []

#file to write points into
tmpfile = open('points_to_convexify.dat', "w")

# Iterate over cells and collect dofs
for cell in cells(mesh):
  cell_ind = cell.index()
  vert_inds = cell.entities(0)
 
  print "cell", cell.index()
  points = [coord[vert_inds[0]], coord[vert_inds[1]], coord[vert_inds[2]], \
            middlePoint(coord[vert_inds[0]], coord[vert_inds[1]]), \
            middlePoint(coord[vert_inds[1]], coord[vert_inds[2]]), \
            middlePoint(coord[vert_inds[2]], coord[vert_inds[0]])]
 
  #calculate lagrange coefficients
  coeffs = map(u_,points)
  #print(coeffs)

#for i in range(len(u_array)):
#  dof_to_vertex_map[0]
#  #points = [ dof_to_vertex_map[j] for j in range(3)]
#  points.append([middlePoint(points[0], points[1]), \
#                 middlePoint(points[1], points[2]), \
#                 middlePoint(points[2], points[0])])
#  coeffs = u_array[6*i:6*i+6]
  
#  print coeffs
 
  #transform into bezier basis  
  coeffBezier = mult(trafo,coeffs)
  
  #write into file
  for i in range(len(points)):
    tmpfile.write(str(points[i][0]) + " " + str(points[i][1]) + " " + str(coeffBezier[i]) + '\n')

#close file after writing points  
tmpfile.close()


#wait for convexification
#input('Press bar after the Convexificaton process!')
os.system("../bin/Convexify points_to_convexify.dat convexified_points.dat")

#open file with convexified coefficients
tmpfile = open('convexified_points.dat', 'r')

#create array
u_array = u_.vector().array()
print len(u_array)

#update coefficients
i = 0
for line in tmpfile:
  coordinates = line.split(" ")
  u_array[i] = coordinates[2]
  i = i+1

print u_array

#transform from bezier to lagrange basis
for i in range(len(u_array)/6): #loop over all cells
  coeffBezier = u_array[i*6: i*6+6]
  #print i, coeffBezier
  u_array[i*6: i*6+6] = mult(trafo_inv,coeffBezier)

#create new convex function
u_convex = Function(V)

u_convex.vector()[:] = u_array

#plot results
plot(u_convex, title = 'convexified solution')
plot(u_, title = 'solution')
plot(u_convex-u_e, title = 'error in conv. sol.')

#hold plot
interactive()

#exit()

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
prm['newton_solver']['maximum_iterations'] = 100
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
file << u
  
