# -*- coding: utf-8 -*-
"""
"""
#source ~/Work/FEniCS/share/dolfin/dolfin.conf
#, damit dolfin geladen werden kann

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
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

def frobenius_product2(a,b):
  return a[0,0]*b[0] + a[0,1]*b[1] + a[1,0]*b[2] + a[1,1]*b[3]


def determinant(a):
  return ((a[0]*a[3]) - (a[1]*a[2]))

def determinant_sym(a):
  return ((a[0]*a[3]) - ((a[1]+a[2])/2.0)**2)

def matrix_mult(A,b):
  return (as_matrix([[A[0],A[1]],[A[2],A[3]]])*b)
  #return [A[0]*b[0] + A[1]*b[1], A[2]*b[0] + A[3]*b[1]]


def neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, sigmaB, u_):

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  (u, w)= TrialFunctions(W)

  (v, mu)  = TestFunctions(W)

  #--------------define function------------------------
  F = 0 
  #(f-det(DH^2 u))*v
  F  = (f- determinant(w))*v*dx

  #jump in gradient
  F = F + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS
  #F = F + dot(nabla_grad(u)('+') - nabla_grad(u)('-'),n('+'))*v('+')*dS
  
  #jump in function
  F = F + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS

  #test with hessian
  F = F + frobenius_product(w, mu)*dx

  #piecewise hessian

  F = F - frobenius_product2(grad(grad(u)),mu)*dx

  #correction term
  F = F + dot(matrix_mult(avg(mu),nabla_grad(u)('+')) ,n('+'))*dS \
        + dot(matrix_mult(avg(mu),nabla_grad(u)('-')) ,n('-'))*dS

  #boundary conditions
  F = F + Constant(sigmaB)/h*(u-u0)*v*ds

  F = action(F, u_)


  J  = derivative(F, u_)   # Gateaux derivative in dir. of u

  #------define boundary conditions--------
  u0_boundary = Function(W)
  assign(u0_boundary.sub(0), project(u0,V))
  #assign(u0_boundary.sub(1), grad(grad(u0)))
  assign(u0_boundary.sub(1),project(Expression((('1.0','0.0','0.0','1.0'))),Sigma))

  def boundary(x, on_boundary):
      return on_boundary

  bc = DirichletBC(W, u0_boundary, boundary)

  #-----------start solver-----------

  problem = NonlinearVariationalProblem(F, u_, None, J)
  solver  = NonlinearVariationalSolver(problem)

  prm = solver.parameters
  
  prm['nonlinear_solver']='snes'
  prm['snes_solver']['absolute_tolerance'] = 1E-8
  prm['snes_solver']['maximum_iterations'] = 50
  prm['snes_solver']['linear_solver']= 'petsc'
  #prm['snes_solver']['preconditioner']= 'lu'
  prm['snes_solver']['line_search'] = 'basic'
  
  info(prm['snes_solver'], True)
  #prm['newton_solver']['absolute_tolerance'] = 1E-8
  #prm['newton_solver']['relative_tolerance'] = 1E-8
  #prm['newton_solver']['maximum_iterations'] = 1000
  #prm['newton_solver']['relaxation_parameter'] = 0.2
  #prm['newton_solver']['report'] = True

  
  set_log_level(PROGRESS)

  nb_iterations, converged = solver.solve()
  
  #write Newton steps required to file 
  newtonStepsfile.write(str(it)+' '+str(nb_iterations)+'\n')
  
  return u_


def neilan_start(u0, f, Nh, deg, deg_hessian, sigmaC, sigmaG, sigmaB, maxRefinement):
  mesh = UnitSquareMesh(Nh, Nh, "crossed")  

  #space for function
  V = FunctionSpace(mesh, 'DG', deg)
  #space for hessian entries
  Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
  #space for discrete hessian
  Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)

  #combination of functions and its discrete hessian
  W = V*Sigma

  #refined space to plot funcitons refined
  bigMesh = refine(refine(mesh))
  bigV = FunctionSpace(bigMesh, 'DG', deg, 'crossed')

  #--------define startsolution-------------
  u1_ = Function(V)
  u1_ = start_iteration(mesh, V, u0, f, sigmaC, parameters['form_compiler']['quadrature_degree'])
  
  u_ = Function(W)
  assign(u_.sub(0), u1_)
  #calculate current hessian
  assign(u_.sub(1), [project((u1_.dx(0)).dx(0),Sigma_single), \
                     project((u1_.dx(0)).dx(1),Sigma_single), \
                     project((u1_.dx(1)).dx(0),Sigma_single), \
                     project((u1_.dx(1)).dx(1),Sigma_single)])

  for it in range(1,maxRefinement):
    print 'Starting Neilan with ', Nh
    w = neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, sigmaB, u_)

    #examine error
    error_norm = errornorm(u0, w.sub(0))
    print 'Errornorm:', error_norm
 
    error_norm = errornorm(u0, w.sub(0), norm_type='H1')
    print 'Errornorm H1:', error_norm

    #------refine grid---------
    #mesh = refine(mesh)
    Nh = Nh *2
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
    Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)
    W = V*Sigma
    
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    u_=project(w,W)
  
  return Nh, V, bigV, u_

if __name__ == "__main__":
  start = time.clock()
  
  
  parameters['form_compiler']['quadrature_rule'] = 'canonical'
  
  #------------ Create mesh and define function space-----------
  Nh = 2
  
  if len(sys.argv) < 4:
    print 'Error, please specify problem and the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  deg = int(sys.argv[2])
  deg_hessian = int(sys.argv[3])

  print "deg ", deg, " deg_hessian", deg_hessian

  parameters['form_compiler']['quadrature_degree'] = 2*deg


  mesh = UnitSquareMesh(Nh, Nh, "crossed")
  interactive

  #space for function
  V = FunctionSpace(mesh, 'DG', deg)
  #space for hessian entries
  Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
  #space for discrete hessian
  Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)

  #combination of functions and its discrete hessian
  W = V*Sigma

  #refined space to plot funcitons refined
  bigMesh = refine(refine(mesh))
  bigV = FunctionSpace(bigMesh, 'DG', deg, 'crossed')

  #define penalty
  sigmaC = 20*deg*deg
  sigmaG = 50
  sigmaB = 20*deg*deg
  
  fileprefix = problem_name+'_Neilan_GradJump_deg'+str(deg)+str(deg_hessian)+'_'
  print "processing files ", fileprefix
  
  errorfile = open('data/'+fileprefix+'l2errornorm','wa',1)
  errorfile.write('iterations l2error\n');
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa',1)
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa',1)
  newtonStepsfile.write('iterations steps\n');
  
  u0, f = MA_problem(problem_name, Nh, 2*deg, mesh)

  #define exact solution
  u_e = project(u0, V)

  #--------define startsolution-------------
  u1_ = Function(V)
  #u1_ = start_iteration(mesh, V, u0, f, sigmaC, parameters['form_compiler']['quadrature_degree'])
  #u1_.assign(u_e)
  u1_ = project(Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0'), V)
  
  print 'start error ', errornorm(u0, u1_)
  errorfile.write('0 '+str(errornorm(u0, u1_))+'\n')
  errorfileh1.write('0 '+str(errornorm(u0, u1_, norm_type='h1'))+'\n')
  
  u_ = Function(W)
  
  assign(u_.sub(0), u1_)
 
  #calculate current hessian
  #assign(u_.sub(1), [project((u1_.dx(0)).dx(0),Sigma_single), \
  #                   project((u1_.dx(0)).dx(1),Sigma_single), \
  #                   project((u1_.dx(1)).dx(0),Sigma_single), \
  #                   project((u1_.dx(1)).dx(1),Sigma_single)])

  assign(u_.sub(1),project(Expression((('1.0','0.0','0.0','1.0'))),Sigma))

  #plot start solution 
  if False:
    plot(project(u_.sub(0),bigV), title = 'startsolution')
    plot(project(u_.sub(0)-u_e,bigV), title = 'start error')
    plot(determinant(u_.sub(1)), title = 'determinant of hessian')
    plot(project(determinant(u_.sub(1))-f,bigV), title = 'determinant of hessian')
  
    plot(u_.sub(1)[0], title = 'first entry of hessian')
    # plot(project(abs(u_.sub(1)[0]-u0_derivXX_e),bigV), title = 'first entry of hessian error')
    #plot(u_.sub(1)[1], title = 'second entry of hessian')
    #plot(u_.sub(1)[2], title = 'third entry of hessian')
    #plot(u_.sub(1)[3], title = 'fourth entry of hessian')
    #plot(project(abs(u_.sub(1)[3]-u0_derivYY_e),bigV), title = 'first entry of hessian error')

    interactive()
  
  for it in range(1,8):
    print 'Starting Neilan with ', Nh
    
      #-----------choose PDE------------
    u0, f = MA_problem(problem_name, Nh, 2*deg, mesh)
    
    w = neilan_step(mesh, V, Sigma, W, u0, f, sigmaC, sigmaG, sigmaB, u_)

    #examine error
    error_norm = errornorm(u0, w.sub(0), degree_rise=5)
    print 'Errornorm:', error_norm
    errorfile.write(str(it)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, w.sub(0), norm_type='H1', degree_rise=5)
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

#    u00, u01, u10, u11 = (((MA_exact_derivative('MA2'))))
#    plot(project(w.sub(1)[1]-u01,bigV), title = 'first entry of hessian error')
#    interactive()

    # ----Plot solution and mesh-------
    if True:
      s = 'plots/'+fileprefix+'_Nh'+str(Nh)+'.pvd'
      file = File(s)
      solution = Function(bigV, name='u')
      solution.assign(project(w.sub(0),bigV))
      file << solution
  
      s = 'plots/'+fileprefix+'error_Nh'+str(Nh)+'.pvd'
      file = File(s)
      error = Function(bigV, name='error')
      error.assign(project(abs(w.sub(0)-u0),bigV))
      file << error

    if False:
      #plot(mesh, title='mesh'+str(it))
      plot(project(w.sub(0),bigV), title = 'solution'+str(it))
      #plot(project(w.sub(0)-u0,bigV), title = 'error'+str(Nh))
      #plot(project(determinant(w.sub(1))-f,bigV), title = 'determinant of discr. hessian-rhs'+str(Nh))
      #plot(project(det(grad(grad(w.sub(0))))-f,bigV), title = 'determinant of hessian-rhs'+str(Nh))
      #Hold plot
      interactive()

    #------refine grid---------
    #mesh = refine(mesh)
    Nh = Nh *2
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
    Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)
    W = V*Sigma
    
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    u_=project(w,W)
  
    #plot(project(u_.sub(0),bigV), title = 'projected solution'+str(it))

  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()