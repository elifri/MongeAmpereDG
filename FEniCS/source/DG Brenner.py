"""
"""

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from convexify import convexify
import math


if __name__ == "__main__":
  
  # Create mesh and define function space
  deg = 3
  Nh = 2
  
  mesh = UnitSquareMesh(Nh, Nh, 'crossed')
  V = FunctionSpace(mesh, 'CG', deg)
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'CG', deg)
  
  #define penalty
  sigma = 20*deg*deg

  #-------define problem------------
  
  # Define boundary conditions
  #u0 = Constant(0.0) #const rhs
  u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  #u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  #u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  #u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1
  
  #rhs
  f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
  #f = Constant(7.0)#simpleMongeAmpere
  #f = Constant(1.0) #simpleMongeAmpere2
  #f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1
  
  #define exact solution
  u_e = interpolate(u0, V)
  
  

  #define startsolution
  #choose between "identity" and disturbed exact solution
  u_ = start_iteration(mesh, V, u0, f, sigma)
  u_.assign(u_e-0.01*interpolate(error, V))
  
  #====================================
  #define brenner's iteration
  #====================================
  
  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  u  = TrialFunction(V)
  v  = TestFunction(V)
  
  coeff = cofac (grad(grad(u)))
  
  F  = (f-det(coeff) )*v*dx
  
  F = F + (dot(avg(coeff)*(nabla_grad(u)('+')) , n('+')) +  dot(avg(coeff)*(nabla_grad(u)('-')), n('-')) )*avg(v)*dS
  
  F = F - (dot(coeff*(nabla_grad(v)) , n) )*(u-u0)*ds
  
  F = F + Constant(sigma)/h*(u-u0)*v*ds
  
  
  F = action(F, u_)
  
  J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u
  
  def boundary(x, on_boundary):
      return on_boundary
  
  bc = DirichletBC(V, u0, boundary)
  
  problem = NonlinearVariationalProblem(F, u_, bc, J)
  solver  = NonlinearVariationalSolver(problem)
  
  prm = solver.parameters
  #info(prm, True)
  
  prm['newton_solver']['absolute_tolerance'] = 1E-8
  prm['newton_solver']['relative_tolerance'] = 1E-10
  prm['newton_solver']['maximum_iterations'] = 30
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
  plot(project(u_,bigV), title = 'solution')
  plot(project(abs(u_-u_e),bigV), title = 'error')
  
  plot(det(grad(grad(u_))), title = 'determinant of hessian')
  
  plot(project(abs(f-det(grad(grad(u_)))), bigV), title = 'rhs - determinant of hessian')
  
  #plot(mesh)
  
  #Hold plot
  interactive()
  
  
  # Dump solution to file in VTK format
  s = 'MongeAmpere.pvd'
  file = File(s)
  file << u_
    
