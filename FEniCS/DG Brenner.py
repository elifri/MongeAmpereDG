"""
"""

from dolfin import *
import scipy.io

# Create mesh and define function space
mesh = UnitSquareMesh(30, 30)
mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)


# Define boundary conditions
#u0 = Constant(0.0) #const rhs
#u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2

#rhs
#f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
f = Constant(4.0)#simpleMongeAmpere
#f = Constant(1.0) #simpleMongeAmpere2

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

plot(u_)
plot(abs(u_-u_e))

interactive()


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

#question necessary? (fenics implies boundary conditions)
F = F + Constant(sigma)/h*(u-u0)*v*ds


F = action(F, u_)

J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, boundary)

problem = NonlinearVariationalProblem(F, u_, bc, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
info(prm, True)

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
  
