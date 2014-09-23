"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *

def frobenius_product(a,b):
  return a[0,0]*b[0,0] + a[0,1]*b[0,1] + a[1,0]*b[1,0] + a[1,1]*b[1,1]


# Create mesh and define function space
mesh = UnitSquareMesh(6, 4)
V = FunctionSpace(mesh, 'Lagrange', 1)
W = TensorFunctionSpace(mesh, 'DG', 0, shape=(2,2))
P = V*W

# Define boundary conditions
u0 = Constant('0.0')
# Define rhs
class rhs(Expression):    
  def eval(self, v, x):
    v[0] = 1.0

f = rhs(element = FiniteElement("Quadrature", triangle, 3))

# Define variational problem
u = TrialFunction(P)
u1, u2 = split(u)
v1, v2 = TestFunctions(P)

h = CellSize(mesh)
n = FacetNormal(mesh)

F = 0
F = f*v1*dx
F = F - det(u2)*v1*dx
F = F + Constant(50.0)('+')*avg(h)* jump(nabla_grad(u1),n)*jump(nabla_grad(v1),n)*dS
F = F + frobenius_product(u2,v2)*dx
F = F + (dot(avg(v2)*nabla_grad(u1)('+'),n('+'))+dot(avg(v2)*nabla_grad(u1)('-'),n('-')))*dS
F = F + Constant(50.0)/h*(u1-u0)*v1*ds

u_ = Function(P)
F = action(F, u_)

J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u

problem = NonlinearVariationalProblem(F, u_, None, J)
solver  = NonlinearVariationalSolver(problem)

set_log_level(PROGRESS)
solver.solve()


# Plot solution and mesh
#plot(u_.sub(0))
# Hold plot

mesh = refine(mesh)
V = FunctionSpace(mesh, 'Lagrange', 1)
W = TensorFunctionSpace(mesh, 'DG', 0, shape=(2,2))
P = V*W


# Define variational problem
U = TrialFunction(P)
u, w = split(U)
v, mu = TestFunctions(P)

h = CellSize(mesh)
n = FacetNormal(mesh)


F = 0
F = f*v*dx
F = F - det(w)*v*dx
F = F + Constant(50.0)('+')*avg(h)* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS
F = F + frobenius_product(w,mu)*dx
F = F + (dot(avg(mu)*nabla_grad(u)('+'),n('+'))+dot(avg(mu)*nabla_grad(u)('-'),n('-')))*dS
F = F + Constant(50.0)/h*(u-u0)*v*ds
#F = (f-det(w))*v*dx - Constant(50.0)*(u-u0)*v*ds

u_ = Function(P)
F = action(F, u_)

J  = derivative(F, u_, U)   # Gateaux derivative in dir. of u

problem = NonlinearVariationalProblem(F, u_, None, J)
solver  = NonlinearVariationalSolver(problem)

set_log_level(PROGRESS)
solver.solve()


# Plot solution and mesh
plot(u_.sub(0), title='tada')
interactive()