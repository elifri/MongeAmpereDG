"""
"""

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from convexify import convexify
import math

# Create mesh and define function space
n = 20
deg = 3

mesh = UnitSquareMesh(n, n)
#plot(mesh)
V = FunctionSpace(mesh, 'CG', deg, 'crossed')
bigMesh = refine(refine(mesh))
bigV = FunctionSpace(bigMesh, 'CG', deg, 'crossed')

# Define boundary conditions
#u0 = Constant(0.0) #const rhs
u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
#u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1

class Error(Expression):
  def eval(self, v, x):
    s = (x[0]-1./2)**2+(x[1]-1./2)**2
    if s < 1./10:
      v[0] = math.exp(-1./(1-100*s**2))
    else:
      v[0]=0


#error = Expression('1/2*x[0]*x[0]+1/2*x[1]*x[1]')#add noisy data
#error = Expression('0.1*pow(x[0],2)+0.01*pow(x[1],2)')#add noisy data
#error = Expression('0.01*sin(3*x[0]*x[1]+pi)')#add noisy data
error = Error()

print "u0 ", u0

#rhs
f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
#f = Constant(7.0)#simpleMongeAmpere
#f = Constant(1.0) #simpleMongeAmpere2
#f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

#define exact solution

u_e = interpolate(u0, V)

#define startsolution

u_ = Function(V)


#====================================
#define laplace's iteration
#====================================
#u_= MA_iteration(mesh, V, u0, f, 1)
#print "u_array before ", u_array

#==================
#convexify
#===================
"""
#create new convex function
u_convex = Function(V)

u_array = u_.vector().array()
u_convex.vector()[:] = convexify(mesh, V, u_array)
print "u_array ", u_array
print "u_convex_array ", u_convex.vector().array()

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
"""
#hold plot
#interactive()

#update start solution
#u_.assign(u_convex)


#choose between "identity" and disturbed exact solution
u_.assign(u_e-0.01*interpolate(error, V))
#u_.assign(u_e)

#u_ = interpolate(Expression('2*x[0]*x[0]+2*x[1]*x[1]+x[0]*x[1]'),V)

#====================================
#define brenner's iteration
#====================================


#penalty
sigma = 100


#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)

u  = TrialFunction(V)
v  = TestFunction(V)

coeff = cofac (grad(grad(u)))

F  = (f-det(coeff) )*v*dx

F = F + (dot(avg(coeff)*(nabla_grad(u)('+')) , n('+')) +  dot(avg(coeff)*(nabla_grad(u)('-')), n('-')) )*avg(v)*dS

F = F - (dot(coeff*(nabla_grad(v)) , n) )*(u-u0)*ds

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
prm['newton_solver']['relative_tolerance'] = 1E-10
prm['newton_solver']['maximum_iterations'] = 30
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
plot(project(u_,bigV), title = 'solution')
plot(project(abs(u_-u_e),bigV), title = 'error')

plot(det(grad(grad(u_))), title = 'determinant of hessian')

plot(project(abs(f-det(grad(grad(u_)))), bigV), title = 'rhs - determinant of hessian')

#plot(mesh)

#Hold plot
interactive()


# Dump solution to file in VTK format
s = 'MongeAmpere.pvd'
file = File(s)
file << u_
  
