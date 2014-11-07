from dolfin import *
import scipy.io, numpy as np
from MA_problem import *
import math, time, sys

#define frobenius product between to matrices
def frobenius_product(a,b):
  return a[0,0]*b[0,0] + a[0,1]*b[0,1] + a[1,0]*b[1,0] + a[1,1]*b[1,1]


#perform Neilan's FE method, i.e. find root of (3.11)
def neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, u_):
  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  #init ansatz and test functions
  U  = TrialFunction(W)
  Phi  = TestFunction(W)
  
  u,w = split(U)
  v,mu = split(Phi)
  
  
  F  = (f-det(w))*v*dx
  
  #penalise jump in function
  F = F + Constant(sigmaC)('+')/avg(h)*jump(u)*jump(v)*dS

  #penalise jump in Gradient across edges
  F = F + Constant(sigmaG)('+')*avg(h)*jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS
 
  #d_dh u : mu 
  F = F + frobenius_product(w,mu)*dx
  
  #frobenius product with piecewise hessian
  F = F - frobenius_product(grad(grad(u)),mu)*dx

  #correction term for discrete hessian
  F = F + (dot(avg(mu)*nabla_grad(u)('+'),n('+'))+dot(avg(mu)*nabla_grad(u)('-'),n('-')))*dS

  #weak boundary conditions
  F = F + Constant(sigmaC)/h*(u-u0)*v*ds
  
  F = action(F, u_)
  
  # Gateaux derivative in dir. of u
  J  = derivative(F, u_, U)   
  
  problem = NonlinearVariationalProblem(F, u_, None, J)
  solver  = NonlinearVariationalSolver(problem)
  
  prm = solver.parameters
  prm['nonlinear_solver']='snes'
  prm['snes_solver']['absolute_tolerance'] = 1E-8
  prm['snes_solver']['linear_solver']= 'petsc'
  prm['snes_solver']['maximum_iterations'] = 100
  prm['snes_solver']['line_search'] = 'basic' 
  set_log_level(PROGRESS)
  
  nb_iterations, converged = solver.solve()
  
  return u_, converged
  
if __name__ == "__main__":
  Nh = 1
  
  if len(sys.argv) != 4:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  deg_hessian = int(sys.argv[3])

  parameters['form_compiler']['quadrature_degree'] = 2*deg

  initial_mesh = UnitSquareMesh(Nh, Nh, "crossed")
  mesh = adapt(initial_mesh)
  
  #define function spaces
  V = FunctionSpace(mesh, 'DG', deg)
  #space for discrete hessian
  Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))

  #combination of the two functions spaces
  W = V*Sigma

  #define penalty
  sigmaC = 20.0*deg*deg
  sigmaG = 50.0

  f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
  #get start function
  u1_ = start_iteration(mesh, V, u0, f, sigmaC)

  #calculate start functions hessian and write both into the mixed function space
  u_ = Function(W)
  
  assign(u_.sub(0), u1_) 
  #calculate current hessian
  assign(u_.sub(1), project(as_matrix((((u1_.dx(0)).dx(0), \
                     (u1_.dx(0)).dx(1)), \
                     ((u1_.dx(1)).dx(0), \
                     (u1_.dx(1)).dx(1)))), Sigma))

  #-------nested iteration--------------
  for it in range(1,8):
    #perform once the finite element method on current grid
    w, converged = neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, u_)

    if not converged:
      break
    
    #examine error
    error_norm = errornorm(u0, w.sub(0))
    print 'Errornorm:', error_norm

    #------refine grid---------
    mesh = adapt(mesh)
    Nh = Nh*2
    
    #update function spaces
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))
    W = V*Sigma
    
    #project previous solution to new grid
    u_ = Function(W)
    u_=project(w,W)
 
    #update problem information
    f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
