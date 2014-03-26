"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *
#from numpy import *

# Create mesh and define function space
mesh = UnitSquare(10, 10)
mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)

# Define boundary conditions
u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]')
#u0 = Constant(0.0)
u_e = interpolate(u0, V)

def u0_boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

#define variables
#coeff = as_matrix([[ 4, -3],[ -3, 4]])
coeff = as_matrix([[1,0],[0,1]])
sigma = Constant(7.0*20.0)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Constant(8.0)

h = CellSize(mesh)
n = FacetNormal(mesh)

a = -inner(coeff*nabla_grad(u), nabla_grad(v))*dx \
  + inner(jump(v,n),avg(coeff*nabla_grad(u)))*dS\
  + inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS \
  + v*inner(n,coeff*nabla_grad(u))*ds \
  + u*inner(n,coeff*nabla_grad(v))*ds \
  + sigma('+')/h('+')* jump(u)*jump(v)*dS \
  + sigma/h*v*u*ds
  
A = assemble(a)
values = A.array()

print values

L = inner(f,v)*dx + u0*dot(n,coeff*nabla_grad(v))*ds +sigma/h *u0*v*ds
b= assemble(L)
#print b.array()

# Compute solution
u = Function(V)
#solve(a == L, u)
solve(A,u.vector(),b)

#examine error

u_e_array = u_e.vector().array()
u_array = u.vector().array()
print 'Max error:', abs(u_e_array - u_array).max()

# Plot solution and mesh
plot(u)
#plot(u-u_e)
#plot(mesh)

# Dump solution to file in VTK format
file = File('poisson.pvd')
file << u

# Hold plot
interactive()
