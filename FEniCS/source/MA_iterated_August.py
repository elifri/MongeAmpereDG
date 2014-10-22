"""
"""

from dolfin import *
import scipy.io
from convexify import convexify
from MA_iterated_benamour import MA_iteration_Benamou
from MA_problem import *

import math
import time, sys

def MA_iteration(mesh, V, u0, f, max_it,w, sigmaB, sigmaC, sigmaG):
  
  #------define variables--------

  #define cofactor matrix of startsolution's hessian
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
    - v*inner(n,coeff*nabla_grad(u))*ds \
    - u*inner(n,coeff*nabla_grad(v))*ds \
    + Constant(sigmaB)/h*v*u*ds
#    + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \

  #define rhs functional
  L = inner(Constant(-2.0)*f,v)*dx - u0*dot(n,coeff*nabla_grad(v))*ds +Constant(sigmaB)/h *u0*v*ds

  #iterate
  u = Function(V)

  for iteration in range(0,max_it):
    # Compute solution
    solve(a == L, u)

    #examine error
    error_norm = errornorm(u0, u)
    print 'L2 errornorm h=1/'+str(Nh)+'-'+str(iteration)+':', error_norm
    errorfile.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, u, norm_type='H1')
    print 'Errornorm H1h=1/'+str(Nh)+'-'+str(iteration)+':', error_norm
    errorfileh1.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')

    #plot solution
    if False: 
      plot(project(u, bigV), title = 'solution'+str(Nh)+'-'+str(iteration))
      plot(project(abs(u-u0),bigV), title = 'error'+str(Nh)+'-'+str(iteration))
      interactive()

    #damping and update w
    w.vector()[:] = (0.7*w.vector().array() + 0.3*u.vector().array())

    coeff = cofac(grad(grad(w)))
  return u;

def start_iteration(mesh, V, u0, f, sigma):
  w = Function(V)
  w = interpolate(Expression('0'),V)
  u = MA_iteration_Benamou(mesh, V, u0, f, sigma, 1,w)
  return u

if __name__ == "__main__":
    
  #store start time
  start = time.clock()
    
  if len(sys.argv) != 3:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  
  parameters['form_compiler']['quadrature_degree']=2*deg
  
  # Create mesh and define function space
  Nh = 2
  mesh = UnitSquareMesh(Nh, Nh, 'crossed')
  V = FunctionSpace(mesh, 'CG', deg)
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'CG', deg)
  
  #define penalty
  sigma = 50

  #-------define problem------------
  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
  #define poblem



  #penalty
  sigmaB = 30.0*deg*deg
  sigmaG = 30.0*deg*deg
  sigmaC = 30.0*deg*deg

  #open files for output
  fileprefix = 'MA1_deg'+str(deg)+'_'
  errorfile = open('data/'+fileprefix+'l2errornorm','wa', 1)
  errorfile.write('iterations l2error\n');
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa', 1)
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa', 1)
  newtonStepsfile.write('iterations steps\n');
  newtonStepsfile.write('0 \n');

  #interpolate of exact solution
  u_e = interpolate(u0, V)

  #start solution
  #choose between "identity" and solution of BFO ansatz (\triangle u = -sqrt(2f)
  u = Function(V)
  #u.assign(u_e)
  u = start_iteration(mesh, V, u0, f, sigmaC)

  errorfile.write('0 '+str(errornorm(u0, u))+'\n')
  errorfileh1.write('0 '+str(errornorm(u0, u, norm_type='h1'))+'\n')

  #copy u to w
  w = Function(V)
  w.assign(u)
  
  #define maximum number of Picard iterations on each grid performed
  max_it = 10
 
  #nested iteration
  for it in range(1,8):

    print "calculating on h=1/", Nh

    #plot start solution after refinement
    if False:
      plot(project(u, bigV), title = 'startsolution')
      plot(project(abs(u-u0),bigV), title = 'starterror')
      plot(det(grad(grad(u))), title = 'determinant of starthessian')
      #Hold plot
      interactive()

    #perform max_it steps of the Picard iteration
    u = MA_iteration(mesh, V, u0, f, max_it,w, sigmaB, sigmaC, sigmaG)

    # Plot solution and mesh
    if False:
      plot(project(u,bigV), title = 'solution'+str(Nh))
      plot(project(abs(u-u0),bigV), title = 'error'+str(Nh))
      plot(mesh, title = 'mesh'+str(Nh))
      #plot(det(grad(grad(u))), title = 'determinant of hessian')
      #Hold plot
      interactive()

    #refine grid
    Nh = Nh*2
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')

    #update sapces
    V = FunctionSpace(mesh, 'CG', deg)
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)

    #project old solution to new grid
    w = Function(V)
    w.assign(project(u,V))
  
  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start))
  
  print "%.2gs" % (end-start)

