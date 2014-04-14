"""
"""

from dolfin import *
import scipy.io

# Create mesh and define function space
mesh = UnitSquareMesh(10, 10)
mesh = refine(mesh)
V = FunctionSpace(mesh, 'CG', 2)

#define poblem

# Define boundary conditions
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Constant(0.0) #const rhs
u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2

#rhs
#f = Constant(-8.0) #simpleMongeAmpere
#f = Constant(-2.0) #const rhs
f = Constant(-2.0) #simpleMongeAmpere2

#define exact solution
u_e = interpolate(u0, V)

#define variables

#start solution
w = Function(V)
#w = interpolate(u0, V)
w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
#and the cofactor matrix of its hessian

#penalty
sigma = 7.0*50.0*10

def tensor_jump(u, n):
    return outer(u('+'), n('+')) + outer(u('-'), n('-'))
def my_cofac(u):
#    return as_matrix[ [u.dx(1).dx(1) , -u.dx(1).dx(0)], [-u.dx(0).dx(1) , u.dx(0).dx(0)]]
    return cofac (grad(grad(u)))
def avg_cofac(u):
#    return (my_cofac(u)('+') + my_cofac(u)('-') )
    return avg(my_cofac(u))


#define geometry for penalty and normals
h = CellSize(mesh)
n = FacetNormal(mesh)


u  = TrialFunction(V)
v  = TestFunction(V)
u_ = Function(V)      # the most recently computed solution

coeff = avg(cofac (grad(grad(u))))

print type(dot(coeff*nabla_grad(u) , as_vector((n('+')[0], n('+')[1]))))
print rank(jump(grad(u) , n))
print rank(inner(f-det(grad(grad(u))),v))

F  = (f-det(grad(grad(u))) )*v*dx

#F = 0

F = F + (dot(coeff*(nabla_grad(u)('+')) , n('+')) +  dot(avg_cofac(u)*(nabla_grad(u)('-')), n('-')) )*avg(v)*dS
#  + avg( cofac(grad(grad(u))) ) * jump(nabla_grad(u),n)*avg(v)*dS\
#  + (dot(avg_cofac(u)*(nabla_grad(u)('+')) , n('+')) +  dot(avg_cofac(u)*(nabla_grad(u)('-')), n('-')) )*dS\

F = F + (dot(coeff*(nabla_grad(v)('+')) , n('+')) +  dot(avg_cofac(u)*(nabla_grad(v)('-')), n('-')) )*avg(u-u0)*dS

#  - (dot(nabla_grad(v)('+') , n('+')) +  dot((nabla_grad(v)('-')), n('-')) )*avg(u-u0)*dS\
#  - (dot(avg(coeff)*nabla_grad(v)('+') , n('+')) +  dot(avg_cofac(u)*(nabla_grad(v)('-')), n('-')) )*avg(u-u0)*dS\  
#  - jump( cofac(grad(grad(u)))*nabla_grad(v),n )*avg(u-u0)*ds\

F = F + Constant(sigma)/h*(u-u0)*v*ds


F = action(F, u_)

J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, boundary)

problem = NonlinearVariationalProblem(F, u_, bc, J)
solver  = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0
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
u_array = u.vector().array()
print 'Errornorm:', errornorm(u,u_e)

 # Plot solution and mesh
plot(u)
plot(u-u_e)
#plot(mesh)

# Dump solution to file in VTK format
s = 'poisson'+str(iteration)+'.pvd'
file = File(s)
file << u
  
#Hold plot
interactive()
