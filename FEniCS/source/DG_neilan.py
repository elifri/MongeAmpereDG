"""
"""
#source ~/Work/FEniCS/share/dolfin/dolfin.conf
#, damit dolfin geladen werden kann

from dolfin import *
import scipy.io
import numpy as np
from MA_iterated import MA_iteration
from MA_iterated_August import start_iteration
from convexify import convexify
import math
import time


def frobenius_product(a,b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

def frobenius_product2(a,b):
  return a[0,0]*b[0] + a[0,1]*b[1] + a[1,0]*b[2] + a[1,1]*b[3]


def determinant(a):
  return ((a[0]*a[3]) - ((a[1]+a[2])/2.0)**2)

def matrix_mult(A,b):
  return (as_matrix([[A[0],A[1]],[A[2],A[3]]])*b)
# return [A[0]*b[0] + A[1]*b[1], A[2]*b[0] + A[3]*b[1]]


def neilan_step(mesh, V, W, u0, sigmaC, sigmaG, sigmaB, u_):


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
  #F = F + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS
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
  assign(u0_boundary.sub(0), interpolate(u0,V))
  #assign(u0_boundary.sub(1), grad(grad(u0)))
  assign(u0_boundary.sub(1),interpolate(Expression((('1.0','0.0','0.0','1.0'))),Sigma))

  def boundary(x, on_boundary):
      return on_boundary

  bc = DirichletBC(W, u0_boundary, boundary)

  #-----------start solver-----------

  problem = NonlinearVariationalProblem(F, u_, None, J)
  solver  = NonlinearVariationalSolver(problem)

  prm = solver.parameters
  #info(prm, True)

  prm['newton_solver']['absolute_tolerance'] = 1E-8
  prm['newton_solver']['relative_tolerance'] = 1E-8
  prm['newton_solver']['maximum_iterations'] = 10
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

  nb_iterations, converged = solver.solve()
  
  #write Newton steps required to file 
  newtonStepsfile.write(str(it)+' '+str(nb_iterations)+'\n')
  
  return u_


if __name__ == "__main__":
  start = time.clock()
  
  #------------ Create mesh and define function space-----------
  Nh = 2
  deg = 2
  deg_hessian = 2


  print "deg ", deg, " deg_hessian", deg_hessian

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

  #define penalty
  sigmaC = 10*deg*deg
  sigmaG = 10*deg*deg
  sigmaB = 10*deg*deg

  #-----------choose PDE------------
  #define rhs
  f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
  #f = Constant(7.0)#simpleMongeAmpere
  #f = Constant(1.0) #simpleMongeAmpere2
  #f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1

  # Define boundary conditions
  #u0 = Constant(0.0) #const rhs
  u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  #u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  #u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  #u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1

  fileprefix = 'MA1_Neilan_deg'+str(deg)+str(deg_hessian)+'_'
  errorfile = open('data/'+fileprefix+'l2errornorm','wa')
  errorfile.write('iterations l2error\n');
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa')
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa')
  newtonStepsfile.write('iterations steps\n');
  

  #define exact solution
  u_e = project(u0, V)

  #--------define startsolution-------------
  u1_ = Function(V)
  u1_ = start_iteration(mesh, V, u0, f, sigmaC)
  #u1_.assign(u_e)
  
  errorfile.write('0 '+str(errornorm(u0, u1_))+'\n')
  errorfileh1.write('0 '+str(errornorm(u0, u1_, norm_type='h1'))+'\n')

  
  u_ = Function(W)
  
  assign(u_.sub(0), u1_)
 
  #calculate current hessian
  assign(u_.sub(1), [project((u1_.dx(0)).dx(0),Sigma_single), \
                     project((u1_.dx(0)).dx(1),Sigma_single), \
                     project((u1_.dx(1)).dx(0),Sigma_single), \
                     project((u1_.dx(1)).dx(1),Sigma_single)])

#  assign(u_.sub(1), [u0_derivXX_e, u0_derivXY_e, u0_derivXY_e, u0_derivYY_e])



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
    w = neilan_step(mesh, V, W, u0, sigmaC, sigmaG, sigmaB, u_)

    #examine error
    error_norm = errornorm(u0, u_.sub(0))
    print 'Errornorm:', error_norm
    errorfile.write(str(it)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, u_.sub(0), norm_type='H1')
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

    # ----Plot solution and mesh-------
    s = 'plots/'+fileprefix+'_Nh'+str(Nh)+'.pvd'
    file = File(s)
    solution = Function(bigV, name='u')
    solution.assign(project(u_.sub(0),bigV))
    file << solution

    s = 'plots/'+fileprefix+'error_Nh'+str(Nh)+'.pvd'
    file = File(s)
    error = Function(bigV, name='error')
    error.assign(project(abs(w.sub(0)-u0),bigV))
    file << error


    if False:
      plot(mesh, title='mesh'+str(it))
      #plot(project(w.sub(0),bigV), title = 'solution'+str(it))
      plot(project(abs(w.sub(0)-u0),bigV), title = 'error'+str(Nh))
      #plot(project(determinant(w.sub(1))-f,bigV), title = 'determinant of discr. hessian-rhs'+str(Nh))
      #plot(project(det(grad(grad(w.sub(0))))-f,bigV), title = 'determinant of hessian-rhs'+str(Nh))
      #Hold plot
      interactive()

    #------refine grid---------
    #mesh = refine(mesh)
    Nh = Nh *2
    mesh = UnitSquareMesh(Nh, Nh, 'left/right')
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
    Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)
    W = V*Sigma
    
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    u_.assign(project(w,W))
  
    #plot(project(u_.sub(0),bigV), title = 'projected solution'+str(it))

  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start))
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()