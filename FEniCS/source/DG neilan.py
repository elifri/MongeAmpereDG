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

def frobenius_product(a00, a10, a01, a11, b00, b10, b01, b11):
  return a00*b00 + a10*b10 + a01*b01 + a11*b11 


# Create mesh and define function space
n = 20
deg = 1
deg_hessian = 0

mesh = UnitSquareMesh(n, n)
#plot(mesh)
V = FunctionSpace(mesh, 'CG', deg, 'crossed')
Sigma = FunctionSpace(mesh, 'DG', deg_hessian, 'crossed')

W = MixedFunctionSpace([V,Sigma, Sigma, Sigma, Sigma])

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

u_ = Function(W)

#choose between "identity" and disturbed exact solution
#u_.assign(u_e-0.01*interpolate(error, V))
u_.assign(u_e)

#u_ = interpolate(Expression('2*x[0]*x[0]+2*x[1]*x[1]+x[0]*x[1]'),V)

#====================================
#define brenner's iteration
#====================================

#penalty
sigmaC = 50
sigmaG = 50
sigmaB = 50


#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)

uBig= TrialFunction(W)
u, w00, w10, w01, w11 = split(uBig)
(v, mu00, mu10, mu01, mu11)  = TestFunctions(W)
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


mu_matrix = as_matrix([[mu00,mu01],[mu10,mu11]])

F = 0 
#(f-det(DH^2 u))*v
F = u*v*dx
#F  = (f- (w00*w11-w01*w10))*v*dx

#jump in gradient
#F =  Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS

#jump in function
#F = F + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS

#test with hessian
#F = F + frobenius_product(w00, w10, w01, w11, mu00, mu10, mu01, mu11)*dx

#correction term
#F = F - dot(avg(mu_matrix)*nabla_grad(u)('+'),n('+'))*dS \
#      - dot(avg(mu_matrix)*nabla_grad(u)('-'),n('-'))*dS

#boundary conditions
#F = F + Constant(sigmaB)/h*(u-u0)*v*ds


F = action(F, u_)

J  = derivative(F, u_, uBig)   # Gateaux derivative in dir. of u

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
  
