"""
"""#source ~/Work/FEniCS/share/dolfin/dolfin.conf, damit dolfin geladen werden kann

from dolfin import *
import scipy.io
import math

def MA_iteration_Benamou(mesh, V, u0, f, sigma, max_it,w):
  
  #define variables

  #define exact solution
  u_e = interpolate(u0, V)
  
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  #define bilinear form
  a = inner(nabla_grad(v), nabla_grad(u))*dx \
    - inner(jump(v,n),avg(nabla_grad(u)))*dS\
    - inner(jump(u,n),avg(nabla_grad(v)))*dS\
    + Constant(sigma)('+')/h('+')* jump(u)*jump(v)*dS \
    + Constant(sigma/2.0)('+')/h('+')*jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
    - v*inner(n,nabla_grad(u))*ds \
    - u*inner(n,nabla_grad(v))*ds \
    + Constant(sigma)/h*v*u*ds

  #  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds
  # + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \


  #define rhs functional
  rhs = div(nabla_grad(w))*div(nabla_grad(w))+2*(f-det(grad(grad(w)))) #div(grad(w))*div(grad(w))+2(-f+det(grad(grad(w))))
  L = inner(-sqrt(rhs),v)*dx - u0*dot(n,nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds


  #iterate
  u = Function(V)

  for iteration in range(0,max_it):

    # Compute solution
    solve(a == L, u)
    #solve(A,u.vector(), b)

    #examine error
    print 'Errornorm:', errornorm(u0,u)

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
#    s = 'MongeAmpereAwanou'+str(iteration)+'.pvd'
#    file = File(s)
#    file << u
    
    #update w
    w.assign(u)
    
    #get information abuot second derivatives
    if False:
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
    
      #plot(project(u,bigV), title = 'solution'+str(iteration)+' Benamou')
      #plot(project(abs(u-u_e),bigV), title = 'error'+str(iteration)+' Benamou')
  
      #plot(det(grad(grad(u))), title = 'determinant of hessian'+str(iteration)+' Benamou')
      #plot(abs(wxx-(u_e.dx(0)).dx(0)), title = 'error in dxx'+str(iteration))
  
      #plot(project(abs(f-det(grad(grad(u)))), bigV), title = 'rhs - determinant of hessian'+str(iteration))
  
  return u;
  

class Error(Expression):
  def eval(self, v, x):
    s = (x[0]-1./2)**2+(x[1]-1./2)**2
    if s < 1./10:
      v[0] = math.exp(-1./(1-100*s**2))
    else:
      v[0]=0

if __name__ == "__main__":
  # Create mesh and define function space
  deg = 2
  
  mesh = UnitSquareMesh(40, 40, 'crossed')
  V = FunctionSpace(mesh, 'DG', deg)
  bigMesh = refine(refine(mesh))
  bigV = FunctionSpace(bigMesh, 'DG', deg)


  #define poblem

  # Define boundary conditions
  #u0 = Constant(0.0) #const rhs
  #u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  #u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  #u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1


  #rhs
  #f = Constant(1.0) #const rhs
  #f = Constant(7.0) #simpleMongeAmpere
  #f = Constant(1.0) #simpleMongeAmpere2
  f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
  #f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

  #exact solution
  u_e = interpolate(u0, V)

  #start solution
  w = Function(V)
  error = Error()
  
  #choose between "identity" and disturbed exact solution
  w.assign(u_e-0.1*interpolate(error, V))
  #w.assign(u_e)  
  
  #maximum number of iterations
  max_it = 30

  u = MA_iteration_Benamou(mesh, V, u0, f, max_it,w)

 # Plot solution and mesh
  #plot(project(u,bigV), title = 'solution')
  plot(project(abs(u-u_e),bigV), title = 'error')

  #plot(det(grad(grad(u))), title = 'determinant of hessian')

  #plot(project(abs(f-det(grad(grad(u)))), bigV), title = 'rhs - determinant of hessian')

  #Hold plot
  if (iteration % 10 == 0):
    interactive()