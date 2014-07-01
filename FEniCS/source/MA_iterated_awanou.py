"""
"""#source ~/Work/FEniCS/share/dolfin/dolfin.conf, damit dolfin geladen werden kann

from dolfin import *
import scipy.io

def MA_iteration_Awanou(mesh, V, u0, f, max_it,w):
  #define variables

  #define exact solution
  u_e = interpolate(u0, V)
  
  #cofactor matrix of startsolution's hessian
  coeff = cofac(grad(grad(w)))
  #coeff = as_matrix([[1,0],[0,1]])
  
  print(coeff)


  #penalty
  sigma = 50.0 

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  #define bilinear form
  a = inner(nabla_grad(v), coeff*nabla_grad(u))*dx \
    - inner(jump(v,n),avg(coeff*nabla_grad(u)))*dS\
    - inner(jump(u,n),avg(coeff*nabla_grad(v)))*dS\
    + Constant(sigma)('+')/h('+')* jump(u)*jump(v)*dS \
    + Constant(sigma/2)('+')/h('+')*jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
    - v*inner(n,coeff*nabla_grad(u))*ds \
    - u*inner(n,coeff*nabla_grad(v))*ds \
    + Constant(sigma)/h*v*u*ds

  #  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds
  # + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \
#
#    

  print(div(grad(w)))

  #define rhs functional
  L = inner(-f-div(nabla_grad(w))+det(grad(grad(w))),v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds


  #iterate
  u = Function(V)

  for iteration in range(0,max_it):
    #dump matrices
    A, b = assemble_system(a,L)
    print ("b", b.array())

  #  scipy.io.savemat('A'+str(iteration)+'.mat', {'A': A.array(), 'b': b.array()})

    # Compute solution
    solve(a == L, u)#, solver_parameters={"linear_solver":"bicgstab", "preconditioner":"jacobi"})
    #solve(A,u.vector(), b)

    #examine error
    u_e_array = u_e.vector().array()
    u_array = u.vector().array()
    print 'Errornorm:', errornorm(u,u_e)

    #==================
    #convexify
    #===================
    #convexify function

#    u.vector()[:] = convexify(mesh, V, u_array)
#    print "u_array ", u_array
#    print "u_convex_array ", u_convex.vector().array()

    #plot results

 #   plot(project(u_convex,bigV), title = 'convexified solution')
    #plot(project(u_-u_convex, bigV), title = 'difference convex and not convex')
    #plot(u_-u_convex, title = 'difference convex and not convex')

  #  plot(project(abs(u_-u_e),bigV), title = 'error in non convexified solution')
    
    #hold plot
    #interactive()

    # Dump solution to file in VTK format
    s = 'MongeAmpereAwanou'+str(iteration)+'.pvd'
    file = File(s)
    file << u
    
    #update w
    w.assign(u)
    
    #get information abuot second derivatives
    
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

    #scipy.io.savemat('wxx'+str(iteration)+'.mat', {'wxx': wxx.vector().array()})
    #scipy.io.savemat('wxy'+str(iteration)+'.mat', {'wxy': wxy.vector().array()})
    #scipy.io.savemat('wyx'+str(iteration)+'.mat', {'wyx': wyx.vector().array()})
    #scipy.io.savemat('wyy'+str(iteration)+'.mat', {'wyy': wyy.vector().array()})
  
    plot(project(u,bigV), title = 'solution'+str(iteration)+' Awanou')
    plot(project(abs(u-u_e),bigV), title = 'error'+str(iteration)+' Awanou')

    plot(det(grad(grad(u))), title = 'determinant of hessian'+str(iteration)+' Awanou')
    #plot(abs(wxx-(u_e.dx(0)).dx(0)), title = 'error in dxx'+str(iteration))

    #plot(project(abs(f-det(grad(grad(u)))), bigV), title = 'rhs - determinant of hessian'+str(iteration))
    
    #interactive()
    
  return u;
  

if __name__ == "__main__":
  # Create mesh and define function space
  deg = 2
  
  mesh = UnitSquareMesh(20, 20, 'crossed')
  V = FunctionSpace(mesh, 'DG', deg)
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'DG', deg)


  #define poblem

  # Define boundary conditions
  #u0 = Constant(0.0) #const rhs
  #u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  #u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  #u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1


  #rhs
  #f = Constant(1.0) #const rhs
  #f = Constant(7.0) #simpleMongeAmpere
  f = Constant(1.0) #simpleMongeAmpere2
  #f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
  #f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

  #exact solution
  u_e = interpolate(u0, V)

  #start solution
  w = Function(V)
  #error = Expression('x[0]*x[1]*(1-x[0]) + x[0]*x[1]*(1-x[1])')
  error = Expression('0.0')
  
  #choose between "identity" and disturbed exact solution
  #w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
  w.assign(u_e)
  
  
  #maximum number of iterations
  max_it = 10

  u = MA_iteration_Awanou(mesh, V, u0, f, max_it,w)

 # Plot solution and mesh
  #plot(project(u,bigV), title = 'solution')
  plot(project(abs(u-u_e),bigV), title = 'error')

  #plot(det(grad(grad(u))), title = 'determinant of hessian')

  #plot(project(abs(f-det(grad(grad(u)))), bigV), title = 'rhs - determinant of hessian')

  #Hold plot
  #interactive()
