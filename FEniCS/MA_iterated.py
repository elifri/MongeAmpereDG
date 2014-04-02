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

#define poblem

# Define boundary conditions
u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]')
#u0 = Constant(0.0)

#rhs
f = Constant(8.0) 

#define exact solution
u_e = interpolate(u0, V)

#define variables

#start solution
w = Function(V)
w = interpolate(u0, V)
#and the cofactor matrix of its hessian
#coeff = as_matrix(grad(grad(w)))
coeff = cofac (grad(grad(w)))

#create output in right space
wxx = project(coeff[0,0], FunctionSpace(mesh, 'DG',0))
wyx = project(coeff[1,0], FunctionSpace(mesh, 'DG',0))
wxy = project(coeff[0,1], FunctionSpace(mesh, 'DG',0))
wyy = project(coeff[1,1], FunctionSpace(mesh, 'DG',0))

#plot(f_x)

print wxx.vector().array()
print wyx.vector().array()
print wxy.vector().array()
print wyy.vector().array()


#coeff = as_matrix([[ 4, -3],[ -3, 4]])
#coeff = as_matrix([[1,0],[0,1]])

#penalty
sigma = Constant(7.0*20.0)

#maximum number of iterations
max_it = 1

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)

#define bilinear form
a = -inner(as_matrix(coeff)*nabla_grad(u), nabla_grad(v))*dx \
  + inner(jump(v,n),avg(coeff*nabla_grad(u)))*dS\
  + inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS\
  + sigma('+')/h('+')* jump(u)*jump(v)*dS \
  + v*inner(n,as_matrix(grad(grad(w)))*nabla_grad(u))*ds \
  + u*inner(n,as_matrix(grad(grad(w)))*nabla_grad(v))*ds \
  + sigma/h*v*u*ds

#define rhs functional
L = inner(f,v)*dx + u0*dot(n,coeff*nabla_grad(v))*ds +sigma/h *u0*v*ds


#iterate
u = Function(V)

for iteration in range(0,max_it):

  #update penalty
  sigma = sigma*(iteration+1)*10;

  #dump matrices
#  A, b = assemble_system(a,L)

#  scipy.io.savemat('A'+str(iteration)+'.mat', {'A': A.array(), 'b': b.array()})

  # Compute solution
  solve(a == L, u)#, solver_parameters={"linear_solver":"bicgstab", "preconditioner":"jacobi"})
  #solve(A,u.vector(), b)

  #examine error
  u_e_array = u_e.vector().array()
  u_array = u.vector().array()
  print 'Errornorm:', errornorm(u,u_e)

  # Plot solution and mesh
  plot(u)
  #plot(u-u_e)
  #plot(mesh)

  # Dump solution to file in VTK format
  s = 'poisson'+str(iteration)+'.pvd'
  file = File(s)
  file << u
  
  #Hold plot
  interactive()

  #temp_array = u.vector().array() #copy u's coefficients
  #w.vector().set_local(temp_array)
  
  #w=u.copy(deepcopy = True)
  
  #update w
  w.assign(u)
  
  print coeff

