"""
"""
#source ~/Work/FEniCS/share/dolfin/dolfin.conf
#, damit dolfin geladen werden kann

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from convexify import convexify
import math

def frobenius_product(a,b):
  return a[0,0]*b[0,0] + a[1,0]*b[1,0] + a[0,1]*b[0,1] + a[1,1]*b[1,1]


# Create mesh and define function space
n = 20
deg = 1
deg_hessian = 0

mesh = UnitSquareMesh(n, n)
#plot(mesh)
V = FunctionSpace(mesh, 'CG', deg)
Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))

W = V*Sigma

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


#u_ = interpolate(Expression('2*x[0]*x[0]+2*x[1]*x[1]+x[0]*x[1]'),V)

#====================================
#define brenner's iteration
#====================================

#define startsolution

u_ = Function(W)

#choose between "identity" and disturbed exact solution
u1_ = Function(V)
u1_.assign(u_e-1*interpolate(error, V))
assign(u_.sub(0), u1_)
assign(u_.sub(1),project(as_matrix([[1,0],[0,1]]), Sigma))
#assign(u_.sub(1), interpolate(Expression((('1.0','0.0'),('0.0','1.0'))),Sigma))
# [Constant(1.), Constant(0.), Constant(0.), Constant(1.)])



#penalty
sigmaC = 50
sigmaG = 50
sigmaB = 50


#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)

(u, w)= TrialFunctions(W)

print type(w)

(v, mu)  = TestFunctions(W)
#u = Function(V)
#w00 = Function(Sigma)
#w10 = Function(Sigma)
#w01 = Function(Sigma)
#w11 = Function(Sigma)

#v = TrialFunction(V)
#mu00 = TrialFunction(Sigma)
#mu10 = TrialFunction(Sigma)
#mu01 = TrialFunction(Sigma)
#mu11 = TrialFunction(Sigma)


#mu_matrix = as_matrix([[mu[0,0],mu[0,1],[mu[1,0],mu[1,1]]]])

F = 0 
#(f-det(DH^2 u))*v
#F = u*v*dx
F  = (f- det(w))*v*dx

#jump in gradient
#F = F + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS

#jump in function
F = F + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS

#test with hessian
#F = F + frobenius_product(w, mu)*dx

#correction term
#F = F - dot(avg(mu)*nabla_grad(u)('+'),n('+'))*dS \
#      - dot(avg(mu)*nabla_grad(u)('-'),n('-'))*dS

#boundary conditions
F = F + Constant(sigmaB)/h*(u-u0)*v*ds



F = action(F, u_)

J  = derivative(F, u_)   # Gateaux derivative in dir. of u

u0_boundary = Function(W)
assign(u0_boundary.sub(0), interpolate(u0,V))
#assign(u0_boundary.sub(1), grad(grad(u0)))
assign(u0_boundary.sub(1),project(as_matrix([[1,0],[0,1]]), Sigma))

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(W, u0_boundary, boundary)


problem = NonlinearVariationalProblem(F, u_, None, J)
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
print 'Errornorm:', errornorm(u_.sub(0),u_e)

 # Plot solution and mesh
plot(project(u_.sub(0),bigV), title = 'solution')
plot(project(abs(u_.sub(0)-u_e),bigV), title = 'error')

plot(det(grad(grad(u_.sub(0)))), title = 'determinant of hessian')

plot(project(abs(f-det(grad(grad(u_.sub(0))))), bigV), title = 'rhs - determinant of hessian')

#plot(mesh)

#Hold plot
interactive()


# Dump solution to file in VTK format
s = 'MongeAmpere.pvd'
file = File(s)
file << u_
  
