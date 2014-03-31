"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *
import scipy.io

# Create mesh and define function space
mesh = UnitSquareMesh(1, 1)
mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)

# Define boundary conditions
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]')
u0 = Constant(0.0)
u_e = interpolate(u0, V)

def u0_boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

#define variables
#coeff = as_matrix([[ 4, -3],[ -3, 4]])
coeff = as_matrix([[1,0],[0,1]])
sigma = Constant(7.0*20.0*10)
#sigma = sigma*1*10
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Constant(8.0)

h = Constant(1.0)#CellSize(mesh)
n = FacetNormal(mesh)

print(h)


a = -inner(coeff*nabla_grad(u), nabla_grad(v))*dx \
  + jump(v)*avg(inner(n,coeff*nabla_grad(u)))*dS\
  + jump(u)*avg(inner(n,coeff*nabla_grad(v)))*dS \
  + sigma('+')/h('+')* jump(u)*jump(v)*dS \
  + v*inner(n,coeff*nabla_grad(u))*ds \
  + u*inner(n,coeff*nabla_grad(v))*ds \
  + sigma/h*v*u*ds

#old version:  + inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS \ 
  
#A = assemble(a)
#values = A.array()

L = inner(f,v)*dx + u0*dot(n,coeff*nabla_grad(v))*ds +sigma/h *u0*v*ds
#b= assemble(L)

A, b = assemble_system(a,L)

scipy.io.savemat('../compare_poisson/A.mat', {'A': A.array(), 'b': b.array()})


#print b.array()

# Compute solution
u = Function(V)
solve(a == L, u)#, solver_parameters={"linear_solver":"bicgstab", "preconditioner":"jacobi"})
#solve(A,u.vector(), b)

#examine error

u_e_array = u_e.vector().array()
u_array = u.vector().array()
print 'Errornorm:', errornorm(u,u_e)

# Plot solution and mesh
plot(u)
#plot(u-u_e)
plot(mesh)

#g = TrialFunction(V)
#g=interpolate(Constant(0.0), V)
#g.vector()[25]=1
#g.vector()[0]=1
#g.vector() = g_array
#print(g.vector())
#plot(g)

#print g(0,0), " ", g(0.125,0), " ", g(0, 0.125), " ", g(0.25,0.25)

# Dump solution to file in VTK format
file = File('poisson.pvd')
file << u

# Hold plot
#interactive()
