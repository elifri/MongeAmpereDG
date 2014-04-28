"""
"""

from dolfin import *
import scipy.io

def MA_iteration(mesh, V, u0, f, max_it):
  #define variables

  #define exact solution
  u_e = interpolate(u0, V)

  #start solution
  w = Function(V)
  #w = interpolate(u0, V)
  w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
  #and the cofactor matrix of its hessian
  coeff = as_matrix([[1,0],[0,1]])

  #penalty
  sigma = 7.0*50.0*10

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

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
  L = inner(f,v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds

  #iterate
  u = Function(V)

  for iteration in range(0,max_it):

    #update penalty
    sigma = sigma*(iteration+1)*10;

    print sigma
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


    # Dump solution to file in VTK format
    s = 'poisson'+str(iteration)+'.pvd'
    file = File(s)
    file << u
    
    #update w
    w.assign(u)
    
    coeff = cofac(grad(grad(w)))
    
    #create output in right space
    wxx = project(coeff[0,0], FunctionSpace(mesh, 'DG',0))
    wyx = project(coeff[1,0], FunctionSpace(mesh, 'DG',0))
    wxy = project(coeff[0,1], FunctionSpace(mesh, 'DG',0))
    wyy = project(coeff[1,1], FunctionSpace(mesh, 'DG',0))

    #plot(f_x)

    print "iteration"+str(iteration)+ " cofactor entries"
    print wxx.vector().array()
    print wyx.vector().array()
    print wxy.vector().array()
    print wyy.vector().array()

    scipy.io.savemat('wxx'+str(iteration)+'.mat', {'wxx': wxx.vector().array()})
    scipy.io.savemat('wxy'+str(iteration)+'.mat', {'wxy': wxy.vector().array()})
    scipy.io.savemat('wyx'+str(iteration)+'.mat', {'wyx': wyx.vector().array()})
    scipy.io.savemat('wyy'+str(iteration)+'.mat', {'wyy': wyy.vector().array()})
    
    return u;
  


# Create mesh and define function space
mesh = UnitSquareMesh(40, 40)
mesh = refine(mesh)
V = FunctionSpace(mesh, 'DG', 2)


#define poblem

# Define boundary conditions
#u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
#u0 = Constant(0.0) #const rhs
u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2

#rhs
#f = Constant(-8.0) #simpleMongeAmpere
#f = Constant(-2.0) #const rhs
f = Constant(-2.0) #simpleMongeAmpere2

#maximum number of iterations
max_it = 1

u = MA_iteration(mesh, V, u0, f, max_it)

# Plot solution and mesh
plot(u)
#plot(u-u_e)
#plot(mesh)

#Hold plot
interactive()
