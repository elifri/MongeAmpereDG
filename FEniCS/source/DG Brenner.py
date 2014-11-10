"""
"""#source ~/Work/FEniCS/share/dolfin/dolfin.conf, damit dolfin geladen werden kann

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated_August import start_iteration
from DG_neilan import neilan_start
from MA_problem import *

import math
import time
import sys

def Brenner_step(mesh, V, u0, sigma, u_):
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
  
  print 'Start error ', errornorm(u0, u_)
  
  J  = derivative(F, u_, u)   # Gateaux derivative in dir. of u
  
  def boundary(x, on_boundary):
      return on_boundary
  
  bc = DirichletBC(V, u0, boundary)
  
  problem = NonlinearVariationalProblem(F, u_, bc, J)
  solver  = NonlinearVariationalSolver(problem)
  
  prm = solver.parameters
  #info(prm['snes_solver'], True)
  
  prm['nonlinear_solver']='snes'
  
  prm['snes_solver']['absolute_tolerance'] = 1E-8
  prm['snes_solver']['maximum_iterations'] = 50
  prm['snes_solver']['linear_solver']= 'petsc'
  #prm['snes_solver']['preconditioner']= 'lu'
  prm['snes_solver']['line_search'] = 'basic' 
  
  prm['newton_solver']['absolute_tolerance'] = 1E-8
  prm['newton_solver']['relative_tolerance'] = 1E-7
  prm['newton_solver']['maximum_iterations'] = 300
  prm['newton_solver']['relaxation_parameter'] = 0.2
  set_log_level(PROGRESS)
  
  nb_iterations, converged = solver.solve()
  
  #write Newton steps required to file 
  newtonStepsfile.write(str(it)+' '+str(nb_iterations)+'\n')
  
  return u_


if __name__ == "__main__":
  
  start = time.clock()
  # Create mesh and define function space
  deg = 3
  Nh = 2
  
  if len(sys.argv) != 3:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  
  parameters['form_compiler']['quadrature_degree']=2*deg+2
  
  initial_mesh = UnitSquareMesh(1, 1, "crossed")
  mesh = adapt(initial_mesh)

  V = FunctionSpace(mesh, 'CG', deg)
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'CG', deg)
  
  #define penalty
  sigma = 50

  #-------define problem------------
  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
  #define exact solution
  u_e = interpolate(u0, V)
  
  fileprefix = problem_name+'_Brenner_deg'+str(deg)+'_'
  print "processing files ", fileprefix
  
  errorfile = open('data/'+fileprefix+'l2errornorm','wa',1)
  errorfile.write('iterations l2error\n');
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa',1)
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa',1)
  newtonStepsfile.write('iterations steps\n');
  

  #define startsolution
  #choose between "identity" and disturbed exact solution
  u_ = start_iteration(mesh, V, u0, f, sigma)
  #u_=u_e
  
  newtonStepsfile.write('0 0\n')
  errorfile.write('0 '+str(errornorm(u0, u_))+'\n')
  errorfileh1.write('0 '+str(errornorm(u0, u_, norm_type='h1'))+'\n')
  
  
 # deg_hessian=deg
 # Nh, V, bigV, uTemp = neilan_start(u0, f, Nh, deg, deg_hessian, sigma, sigma, sigma, 2)
 # a,b = uTemp.split(True)
  
  #u_.vector()[:]=a.vector().array()
 # u_.assign(a)
  
  for it in range(1,8):
    print 'Starting Brenner with ', Nh
    w = Brenner_step(mesh, V, u0, sigma, u_)
  
    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm:', error_norm
    errorfile.write(str(it)+' '+str(error_norm)+'\n')
 
    error_norm = errornorm(u0, w, norm_type='H1')
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str(it)+' '+str(error_norm)+'\n')
  
    # ----Plot solution and mesh-------
    if True:
      s = 'plots/'+fileprefix+'_Nh'+str(Nh)+'.pvd'
      file = File(s)
      solution = Function(bigV, name='u')
      solution.assign(project(w,bigV))
      file << solution
    
      s = 'plots/'+fileprefix+'error_Nh'+str(Nh)+'.pvd'
      file = File(s)
      error = Function(bigV, name='error')
      error.assign(project(abs(w-u0),bigV))
      file << error
    
    #------refine grid---------
    #mesh = refine(mesh)
    Nh = Nh *2
    mesh = adapt(mesh)
    V = FunctionSpace(mesh, 'CG', deg)
    
    u_ = Function(V)
    u_= interpolate(w,V)
  
    #examine error
    #print 'Errornorm1:', errornorm(u0, w)
    #print 'Errornorm2:', errornorm(u0, interpolate(w,V))
    #print 'Errornorm3:', errornorm(u0, u_)
  
  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()
