"""
"""

from dolfin import *
import scipy.io
from convexify import convexify
from MA_iterated_benamour import MA_iteration_Benamou

import math

def MA_iteration(mesh, V, u0, f, max_it,w, sigmaB, sigmaC, sigmaG):
  #define variables

  #define exact solution
  u_e = interpolate(u0, V)
  
  #cofactor matrix of startsolution's hessian
  coeff = cofac(grad(grad(w)))
  #coeff = as_matrix([[1,0],[0,1]])

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  #define bilinear form
  a = inner(coeff*nabla_grad(v), nabla_grad(u))*dx \
    - v('+')*  dot(avg(coeff*nabla_grad(u)),n('+'))*dS \
    - v('-')*  dot(avg(coeff*nabla_grad(u)),n('-'))*dS \
    - u('+')*  dot(avg(coeff*nabla_grad(v)),n('+'))*dS \
    - u('-')*  dot(avg(coeff*nabla_grad(v)),n('-'))*dS \
    + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS \
    + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
    - v*inner(n,coeff*nabla_grad(u))*ds \
    - u*inner(n,coeff*nabla_grad(v))*ds \
    + Constant(sigmaB)/h*v*u*ds
  #  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds
  # + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \
#    - jump(v*avg(coeff*nabla_grad(u)),n)*dS\
 #   - jump(u*avg(coeff*nabla_grad(v)),n)*dS\

  #define rhs functional
  L = inner(Constant(-2.0)*f,v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigmaB)/h *u0*v*ds

  #iterate
  u = Function(V)

  for iteration in range(0,max_it):

    #update penalty
    #sigma = sigma*(iteration+1)*10;

    print sigmaG, sigmaC, sigmaB
    
    # Compute solution
    solve(a == L, u)#, solver_parameters={"linear_solver":"bicgstab", "preconditioner":"jacobi"})
    #solve(A,u.vector(), b)

    #examine error
    u_e_array = u_e.vector().array()
    u_array = u.vector().array()
    print 'Errornorm:', errornorm(u0, u)

    #plot(project(u, bigV), title = 'solution'+str(iteration))
#    plot(project(abs(u-u_e),bigV), title = 'error'+str(iteration))
#    plot(det(grad(grad(u))), title = 'determinant of hessian'+str(iteration))

    # Dump solution to file in VTK format
    s = 'MongeAmpere'+str(iteration)+'.pvd'
    file = File(s)
    file << u
    
    #damping and update w
    #w.assign(u)
    w.vector()[:] = (0.7*w.vector().array() + 0.3*u.vector().array())
    
    coeff = cofac(grad(grad(w)))    
  return u;

class Error(Expression):
  def eval(self, v, x):
    s = (x[0]-1./2)**2+(x[1]-1./2)**2
    if s < 1./10:
      v[0] = math.exp(-1./(1-100*s**2))
    else:
      v[0]=0
  

def start_iteration(mesh, V, u0, f):
  w = Function(V)
  w = interpolate(Expression('0'),V)
  u = MA_iteration_Benamou(mesh, V, u0, f, 1,w)
  return u

if __name__ == "__main__":
  # Create mesh and define function space
  deg = 3
  Nh = 2
  
  mesh = UnitSquareMesh(Nh, Nh, 'crossed')
  V = FunctionSpace(mesh, 'DG', deg)
  bigMesh = refine(mesh)
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

  
  u = start_iteration(mesh, V, u0, f)

  #start solution
  w = Function(V)
  #choose between "identity" and disturbed exact solution
  #w = interpolate(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'),V)
  #error = Expression('x[0]*x[1]*(1-x[0]) + x[0]*x[1]*(1-x[1])')
  error = Error()

  w.assign(u)

  #penalty
  sigmaB = 30.0*deg*deg
  sigmaG = 30.0*deg*deg
  sigmaC = 30.0*deg*deg
  
  #maximum number of iterations
  max_it = 10
 
  for it in range(1,8):
    Nh = Nh*2
    
    print "calculating on h=1/", Nh
        
    if False:
      plot(project(u, bigV), title = 'startsolution')
      plot(project(abs(u-u0),bigV), title = 'starterror')
      plot(det(grad(grad(u))), title = 'determinant of starthessian')
      
    u = MA_iteration(mesh, V, u0, f, max_it,w, sigmaB, sigmaC, sigmaG)

    # Plot solution and mesh
    if it > 0:
      plot(project(u,bigV), title = 'solution'+str(Nh))
      plot(project(abs(u-u0),bigV), title = 'error'+str(Nh))

    #plot(det(grad(grad(u))), title = 'determinant of hessian')
    
    #mesh = refine(mesh)
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    V = FunctionSpace(mesh, 'DG', deg)
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    w = Function(V)
    w.assign(project(u,V))
  
    #plot(project(w,bigV), title = 'projected solution')

    #Hold plot
    interactive()
