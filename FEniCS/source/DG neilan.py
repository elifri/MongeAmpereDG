"""
"""
#source ~/Work/FEniCS/share/dolfin/dolfin.conf
#, damit dolfin geladen werden kann

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from MA_iterated_August import start_iteration
from convexify import convexify
import math

def frobenius_product(a,b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

def frobenius_product2(a,b):
  return a[0,0]*b[0] + a[0,1]*b[1] + a[1,0]*b[2] + a[1,1]*b[3]


def determinant(a):
  return ((a[0]*a[3]) - (a[1]*a[2]))

def matrix_mult(A,b):
  return (as_matrix([[A[0],A[1]],[A[2],A[3]]])*b)
# return [A[0]*b[0] + A[1]*b[1], A[2]*b[0] + A[3]*b[1]]

# Create mesh and define function space
n = 4
deg = 2
deg_hessian = 2

mesh = UnitSquareMesh(n, n, "crossed")
#plot(mesh)
V = FunctionSpace(mesh, 'CG', deg)
Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)

W = V*Sigma

bigMesh = refine(refine(mesh))
bigV = FunctionSpace(bigMesh, 'CG', deg, 'crossed')

#define rhs
f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
#f = Constant(7.0)#simpleMongeAmpere
#f = Constant(1.0) #simpleMongeAmpere2
#f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

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

#define exact solution

u_e = interpolate(u0, V)


#u_ = interpolate(Expression('2*x[0]*x[0]+2*x[1]*x[1]+x[0]*x[1]'),V)

u_ = Function(W)

#define startsolution

#choose between "identity" and disturbed exact solution
u1_ = start_iteration(mesh, V, u0, f)
#u1_.assign(u_e)
#u1_.assign(u_e-1*interpolate(error, V))
assign(u_.sub(0), u1_)
#assign(u_.sub(1),project(as_matrix([[1,0],[0,1]]), Sigma))
#assign(u_.sub(1), interpolate(Expression((('1.0','0.0','0.0','1.0'))),Sigma))

#plot(u1_, bigV, title ='start solution')
#plot(project(abs(u_.sub(0)-u_e),bigV), title = 'error')

assign(u_.sub(1), [project((u1_.dx(0)).dx(0),Sigma_single), \
                   project((u1_.dx(0)).dx(1),Sigma_single), \
                   project((u1_.dx(1)).dx(0),Sigma_single), \
                   project((u1_.dx(1)).dx(1),Sigma_single)])



plot(project(u_.sub(0),bigV), title = 'startsolution')
plot(project(abs(u_.sub(0)-u_e),bigV), title = 'start error')
#plot(determinant(u_.sub(1)), title = 'determinant of hessian')

#plot(u_.sub(1)[0], title = 'first entry of hessian')
#plot(u_.sub(1)[1], title = 'second entry of hessian')
#plot(u_.sub(1)[2], title = 'third entry of hessian')
#plot(u_.sub(1)[3], title = 'fourth entry of hessian')

interactive()  


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
F  = (f- determinant(w))*v*dx

#jump in gradient
#F = F + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS

#jump in function
F = F + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS

#test with hessian
F = F + frobenius_product(w, mu)*dx

#piecewise hessian

F = F - frobenius_product2(grad(grad(u)),mu)*dx

#correction term
F = F + dot(matrix_mult(avg(mu),nabla_grad(u)('+')) ,n('+'))*dS \
      + dot(matrix_mult(avg(mu),nabla_grad(u)('-')) ,n('-'))*dS

#boundary conditions
F = F + Constant(sigmaB)/h*(u-u0)*v*ds

F = action(F, u_)


J  = derivative(F, u_)   # Gateaux derivative in dir. of u

u0_boundary = Function(W)
assign(u0_boundary.sub(0), interpolate(u0,V))
#assign(u0_boundary.sub(1), grad(grad(u0)))
assign(u0_boundary.sub(1),interpolate(Expression((('1.0','0.0','0.0','1.0'))),Sigma))

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(W, u0_boundary, boundary)


problem = NonlinearVariationalProblem(F, u_, None, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
#info(prm, True)

prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-10
prm['newton_solver']['maximum_iterations'] = 50
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
plot(mesh, title='mesh')
plot(project(u_.sub(0),bigV), title = 'solution')
plot(project(abs(u_.sub(0)-u_e),bigV), title = 'error')

interactive()  

plot(determinant(u_.sub(1)), title = 'determinant of hessian')

#plot(project(abs(f-determinant(u_.sub(0))), bigV), title = 'rhs - determinant of hessian')

#plot(mesh)

#Hold plot
interactive()  
