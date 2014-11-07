from dolfin import *
import scipy.io
import numpy as np
from MA_problem import *
import sys, time, math

#Find root of given a mesh, a finite Element space V, boundary conditions u0 and a initial guess u_
def Brenner_step(mesh, V, u0, sigma, u_):
  #mesh size and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  #define trial and test funcion
  u  = TrialFunction(V)
  v  = TestFunction(V)
  
  cofactor = cofac (grad(grad(u)))
  
  #initial variational formulation
  F  = (f-det(cofactor) )*v*dx
  # stabilisation terms
  F = F + (dot(avg(cofactor)*(nabla_grad(u)('+')) , n('+')) +  dot(avg(cofactor)*(nabla_grad(u)('-')), n('-')) )*avg(v)*dS
  # weak boundary conditions
  F = F + Constant(sigma)/h*(u-u0)*v*ds
  
  F = action(F, u_)
  
  # Gateaux derivative in dir. of u
  J  = derivative(F, u_, u)   
  
  #define strong boundary conditions
  def boundary(x, on_boundary):
      return on_boundary
  
  bc = DirichletBC(V, u0, boundary)
  
  #define problem and solver
  problem = NonlinearVariationalProblem(F, u_, bc, J)
  solver  = NonlinearVariationalSolver(problem)
  
  #adjust solver paramters
  prm = solver.parameters
  
  prm['nonlinear_solver']='snes'
  prm['snes_solver']['absolute_tolerance'] = 1E-6
  prm['snes_solver']['maximum_iterations'] = 20
  prm['snes_solver']['linear_solver']= 'petsc'
  prm['snes_solver']['line_search'] = 'bt' 
  set_log_level(PROGRESS)
  
  nb_iterations, converged = solver.solve()
  return u_

if __name__ == "__main__":

  #check number of command line parameters
  if len(sys.argv) != 3:
    print 'Error, please specify the problem and the polynomial degree!'
    sys.exit(-1)
  
  #read problem and polynomial degree
  problem_name = sys.argv[1]
  deg = int(sys.argv[2])

  #define quadrature degree
  parameters['form_compiler']['quadrature_degree']=2*deg
  
  # Create mesh and define function space
  Nh = 1
  initial_mesh = UnitSquareMesh(Nh, Nh, "crossed")
  mesh = adapt(initial_mesh)
  Nh = 2
  V = FunctionSpace(mesh, 'CG', deg)
  
  #define penalty (for weak boundary conditions)
  sigma = 20*deg*deg

  #-------define problem------------
  f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)

  #calculate startsolution (laplace u = -sqrt(f))
  u_ = start_iteration(mesh, V, u0, f, sigma)

  #perform nested iteration
  for it in range(1,8):
    #solve nonlinear system
    w = Brenner_step(mesh, V, u0, sigma, u_)
  
    #examine error
    print 'Errornorm:', errornorm(u0, w)
    
    #------refine grid---------
    Nh = Nh *2
    mesh = adapt(mesh)
    V = FunctionSpace(mesh, 'CG', deg)

    #interpolate old solution to refined function space
    u_ = Function(V)
    u_= interpolate(w,V)
