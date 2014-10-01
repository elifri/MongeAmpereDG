from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from MA_iterated_August import start_iteration
from MA_problem import *
import math
import time
import sys

def frobenius_product(a,b):
  return a[0,0]*b[0,0] + a[0,1]*b[0,1] + a[1,0]*b[1,0] + a[1,1]*b[1,1]



def neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, sigmaB, u_):
  #====================================
  #define Neilan's iteration
  #====================================
  
  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  U  = TrialFunction(W)
  Phi  = TestFunction(W)
  
  u,w = split(U)
  v,mu = split(Phi)
  
  F  = (f-det(w))*v*dx
  
  #jump in function
  F = F + Constant(sigmaC)('+')/avg(h)*jump(u)*jump(v)*dS

  #jump  in Gradient
  #F = F + Constant(sigmaG)('+')*h('+')*dot((nabla_grad(u)('+'))-(nabla_grad(u)('-')) , (nabla_grad(v)('+'))-(nabla_grad(v)('-')))*dS
  F = F + Constant(sigmaG)('+')*avg(h)*jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS
 
  
  #jump d_dh u : mu 
  F = F + frobenius_product(w,mu)*dx
  
  #piecewise hessian
  F = F - frobenius_product(grad(grad(u)),mu)*dx

  #jump d_dh u : mu 
  #F = F + Constant(sigmaC)/avg(h)*jump(avg(mu)*nabla_grad(u),n)*dS
  F = F + (dot(avg(mu)*nabla_grad(u)('+'),n('+'))+dot(avg(mu)*nabla_grad(u)('-'),n('-')))*dS

  #boundary conditions
  F = F + Constant(sigmaC)/h*(u-u0)*v*ds
  
  
  F = action(F, u_)
  
  print 'Start error ', errornorm(u0, u_.sub(0))
  
  J  = derivative(F, u_, U)   # Gateaux derivative in dir. of u
  
  
    #------define boundary conditions--------
  def boundary(x, on_boundary):
      return on_boundary

  bc = DirichletBC(W.sub(0), u0, boundary)

  
  problem = NonlinearVariationalProblem(F, u_, None, J)
  solver  = NonlinearVariationalSolver(problem)
  
  prm = solver.parameters
  #info(prm, True)
  
  prm = solver.parameters
  
  prm['nonlinear_solver']='snes'
  prm['snes_solver']['absolute_tolerance'] = 1E-2
  prm['snes_solver']['linear_solver']= 'petsc'
  #prm['snes_solver']['preconditioner']= 'lu'
  prm['snes_solver']['maximum_iterations'] = 50
  prm['snes_solver']['line_search'] = 'basic' 
  
  prm['newton_solver']['absolute_tolerance'] = 2E-4
  prm['newton_solver']['relative_tolerance'] = 2E-4
  prm['newton_solver']['maximum_iterations'] = 1000
  prm['newton_solver']['relaxation_parameter'] = 0.5
  prm['newton_solver']['report'] = True
  
  set_log_level(PROGRESS)
  
  nb_iterations, converged = solver.solve()
  
  #write Newton steps required to file 
  newtonStepsfile.write(str(it)+' '+str(nb_iterations)+'\n')
  return u_
  
  
  
if __name__ == "__main__":
  start = time.clock()
  
  #------------ Create mesh and define function space-----------
  Nh = 2
  
  if len(sys.argv) != 4:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  deg_hessian = int(sys.argv[3])

  parameters['form_compiler']['quadrature_rule'] = 'canonical'
  parameters['form_compiler']['quadrature_degree'] = 2*deg+2


  fileprefix = problem_name+'_Neilan_GradJump_deg'+str(deg)+str(deg_hessian)+'_'
  print "processing files ", fileprefix
  
  errorfile = open('data/'+fileprefix+'l2errornorm','wa',1)
  errorfile.write('iterations l2error\n');
  fitdata_file = open('data/'+fileprefix+'fitData','wa',1)
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa',1)
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa',1)
  newtonStepsfile.write('iterations steps\n');
  newtonStepsfile.write('0 0\n');


  print "problem", problem_name, "deg ", deg, " deg_hessian", deg_hessian

  mesh = UnitSquareMesh(Nh, Nh, "crossed")
  
    #space for function
  V = FunctionSpace(mesh, 'DG', deg)
  #space for hessian entries
  Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
  #space for discrete hessian
  Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))

  #combination of functions and its discrete hessian
  W = V*Sigma

  #refined space to plot funcitons refined
  bigMesh = refine(refine(mesh))
  bigV = FunctionSpace(bigMesh, 'DG', deg, 'crossed')

  #define penalty
  sigmaC = 50.0*deg*deg
  sigmaG = 50.0*deg
  sigmaB = 50.0*deg*deg

  u1_ = Function(V)

  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
  #define start function
  
  u1_ = start_iteration(mesh, V, u0, f, sigmaC)
  #u1_ = project(Expression('x[0]*x[0]/2.0 +x[1]*x[1]/2.0',degree=1), V)

  print 'start error ', errornorm(u0, u1_)
  errorfile.write('0 '+str(errornorm(u0, u1_))+'\n')
  errorfileh1.write('0 '+str(errornorm(u0, u1_, norm_type='h1'))+'\n')
  fitdata_file.write(str(math.log(1.0/Nh))+' '+str(math.log(errornorm(u0, u1_)))+'\n')
  u_ = Function(W)
  
  assign(u_.sub(0), u1_)
 
  #calculate current hessian
  assign(u_.sub(1), project(as_matrix((((u1_.dx(0)).dx(0), \
                     (u1_.dx(0)).dx(1)), \
                     ((u1_.dx(1)).dx(0), \
                     (u1_.dx(1)).dx(1)))), Sigma))

#  assign(u_.sub(1), project(MA_exact_derivativeEx(problem_name), Sigma))

  for it in range(1,7):
    print 'Starting Neilan with ', Nh
    w = neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, sigmaB, u_)

    #examine error
    error_norm = errornorm(u0, w.sub(0))
    print 'Errornorm:', error_norm
    errorfile.write(str(it)+' '+str(error_norm)+'\n')
    fitdata_file.write(str(math.log(1.0/Nh))+' '+str(math.log(error_norm))+'\n')

    error_norm = errornorm(u0, w.sub(0), norm_type='H1')
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

    # ----Plot solution and mesh-------
    if False:
      s = 'plots/'+fileprefix+'_Nh'+str(Nh)+'.pvd'
      file = File(s)
      solution = Function(bigV, name='u')
      solution.assign(project(w.sub(0),bigV))
      file << solution
  
      s = 'plots/'+fileprefix+'error_Nh'+str(Nh)+'.pvd'
      file = File(s)
      error = Function(bigV, name='error')
      error.assign(project(w.sub(0)-u0,bigV))
      file << error

    if True:
      #w_ = Expression('1.0')
      #wTemp = Function(Sigma_single)
      #wTemp = project(w.sub(1)[0,0], Sigma_single)
      #print 'l2 error in first component', errornorm(w_, wTemp)
      
      w_ = project(MA_exact_derivativeEx(problem_name), Sigma)
      plot(mesh, title='mesh h=1/'+str(Nh))
      #plot(project(u0,bigV), title = 'solution'+str(it))
      plot(w.sub(0)-u0, title = 'error for h=1/'+str(Nh))
      #plot(w_[0,0], title = 'exact derivative'+str(Nh))
      #plot(w.sub(1)[0,0], title = 'numer derivative00 '+str(Nh))
      #plot(w.sub(1)[0,0]-w_[0,0], title = 'error der00 '+str(Nh))
      #plot(w.sub(1)[1,1], title = 'numer derivative11 '+str(Nh))
      #plot(w.sub(1)[1,1]-w_[1,1], title = 'error der11 '+str(Nh))
      #Hold plot

    #------refine grid---------
    #mesh = refine(mesh)
    Nh = Nh *2
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
    Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))
    W = V*Sigma
    
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    u_ = Function(W)
    u_=project(w,W)
      
    (u, w) = TrialFunctions(W)
    (v,mu) = TestFunctions(W)

    g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)

  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  interactive()

  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()

