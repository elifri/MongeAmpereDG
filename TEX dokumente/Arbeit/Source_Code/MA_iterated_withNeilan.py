from dolfin import *
import scipy.io
from MA_iterated_August import start_iteration
import DG_neilan as neilan
import sys, math, time

#Calculates the cofactor matrix of the piecewise Hessian
def calc_CofacpiecewiseHessian(mesh, Sigma, u_):

  dh = Function(Sigma)
  #calculate cofac current hessian
  assign(dh, [project((u_.dx(1)).dx(1),Sigma_single), \
                     project(-(u_.dx(0)).dx(1),Sigma_single), \
                     project(-(u_.dx(1)).dx(0),Sigma_single), \
                     project((u_.dx(0)).dx(0),Sigma_single)])
  return dh

#Calculates the cofactor matrix of the discrete Hessian (Neilan's approach)
def calc_CofacdiscreteHessian(mesh, Sigma, u_):
  
  #define geometry
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  dh = TrialFunction(Sigma)
  mu = TestFunction(Sigma)
  
  #define bilinearform to caculate discrete Hessian
  
  #test with hessian
  a_DH = neilan.frobenius_product(dh, mu)*dx

  #piecewise hessian
  l_DH = neilan.frobenius_product2(grad(grad(u_)),mu)*dx

  #correction term
  l_DH = l_DH - dot(neilan.matrix_mult(avg(mu),nabla_grad(u_)('+')) ,n('+'))*dS \
        - dot(neilan.matrix_mult(avg(mu),nabla_grad(u_)('-')) ,n('-'))*dS
  
  dh = Function(Sigma)
  solve(a_DH == l_DH, dh)
  
  #calc cofactor matrix  
  cofactorM = Function(Sigma)
  assign(cofactorM, [project(dh[3],Sigma_single), \
                     project(Constant(-1.0)*((dh[1]+dh[2])/Constant(2.0)), Sigma_single), \
                     project(Constant(-1.0)*((dh[1]+dh[2])/Constant(2.0)), Sigma_single), \
                     project(dh[0],Sigma_single)])

  #plot(cofactorM[0], title = 'first entry of cofac discrete hessian')  
  #plot(cofactorM[1], title = 'second entry of cofac discrete hessian')
  #interactive()
  
  return cofactorM


def MA_iteration_withNeilan(mesh, V, Sigma, u0, f, max_it,w, sigmaB, sigmaC, sigmaG):

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  for iteration in range(0,max_it):

    #cofactor matrix of startsolution's hessian
    coeff = calc_CofacdiscreteHessian(mesh, Sigma, w)
    #coeff = calc_CofacpiecewiseHessian(mesh, Sigma, w)

    #define bilinear form
    a = inner(neilan.matrix_mult(coeff,nabla_grad(v)) , nabla_grad(u))*dx \
      - v('+')*  dot(avg(neilan.matrix_mult(coeff,nabla_grad(u))),n('+'))*dS \
      - v('-')*  dot(avg(neilan.matrix_mult(coeff,nabla_grad(u))),n('-'))*dS \
      - u('+')*  dot(avg(neilan.matrix_mult(coeff,nabla_grad(v))),n('+'))*dS \
      - u('-')*  dot(avg(neilan.matrix_mult(coeff,nabla_grad(v))),n('-'))*dS \
      + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS \
      + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
      - v*inner(n,neilan.matrix_mult(coeff,nabla_grad(u)))*ds \
      - u*inner(n,neilan.matrix_mult(coeff,nabla_grad(v)))*ds \
      + Constant(sigmaB)/h*v*u*ds
    #  + Constant(sigma)/h*inner(grad(v),n)*inner(grad(u),n)*ds
    # + Constant(sigma)('+')/h('+')* jump(grad(u),n)*jump(grad(v),n)*dS \
    #    - jump(v*avg(coeff*nabla_grad(u)),n)*dS\
   #   - jump(u*avg(coeff*nabla_grad(v)),n)*dS\

    #define rhs functional
    L = inner(Constant(-2.0)*f,v)*dx - u0*dot(n,neilan.matrix_mult(coeff,nabla_grad(v)))*ds +Constant(sigmaB)/h *u0*v*ds

    #iterate
    u_ = Function(V)

  
    #update penalty
    #sigma = sigma*(iteration+1)*10;

    print sigmaG, sigmaC, sigmaB
    #coeff = calc_discreteHessian(mesh, Sigma, w)
    # Compute solution
    solve(a == L, u_)#, solver_parameters={"linear_solver":"bicgstab", "preconditioner":"jacobi"})
    #solve(A,u.vector(), b)
    
    #examine error
    print 'Errornorm:', errornorm(u0, u_)

    if False: 
      plot(project(u_, bigV), title = 'solution'+str(iteration))
      plot(project(abs(u_-u0),bigV), title = 'error'+str(iteration))
      interactive()
#    plot(det(grad(grad(u_))), title = 'determinant of hessian'+str(iteration))

    #interactive()
    
    #damping and update w
    #w.assign(u)
    w.vector()[:] = (0.7*w.vector().array() + 0.3*u_.vector().array())
    

    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm:', error_norm
    errorfile.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, w, norm_type='H1')
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')

    
    coeff = cofac(grad(grad(w)))    
  return w;

class Error(Expression):
  def eval(self, v, x):
    s = (x[0]-1./2)**2+(x[1]-1./2)**2
    if s < 1./10:
      v[0] = math.exp(-1./(1-100*s**2))
    else:
      v[0]=0
  
if __name__ == "__main__":
  start = time.clock()
  
  #read parameters
  if len(sys.argv) != 4:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  deg_hessian = int(sys.argv[3])

  parameters['form_compiler']['quadrature_rule'] = 'canonical'
  parameters['form_compiler']['quadrature_degree'] = 2*deg+2

  fileprefix = problem_name+'_iteratedNeilan_deg'+str(deg)+str(deg_hessian)+'_'
  print "processing files ", fileprefix
  
  errorfile = open('data/'+fileprefix+'l2errornorm','wa',1)
  errorfile.write('iterations l2error\n');
  fitdata_file = open('data/'+fileprefix+'fitData','wa',1)
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa',1)
  errorfileh1.write('iterations h1error\n');
  newtonStepsfile = open('data/'+fileprefix+'newtonSteps','wa',1)
  newtonStepsfile.write('iterations steps\n');
  newtonStepsfile.write('0 0\n');
  mesh = UnitSquareMesh(Nh, Nh, 'crossed')
  V = FunctionSpace(mesh, 'DG', deg)
  
  #space for discrete hessian
  Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))
  
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'DG', deg)

  #define poblem
  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
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
  
  #u = u_e  
  u = start_iteration(mesh, V, u0, f, 50)

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
  max_it = 20
 
  for it in range(1,8):
    
    print "calculating on h=1/", Nh
        
    if False:
      plot(project(u, bigV), title = 'startsolution')
      plot(project(abs(u-u0),bigV), title = 'starterror')
      plot(det(grad(grad(u))), title = 'determinant of starthessian')
      
    u = MA_iteration_withNeilan(mesh, V, Sigma, u0, f, max_it,w, sigmaB, sigmaC, sigmaG)

    #examine error
    error_norm = errornorm(u0, u)
    print 'Errornorm:', error_norm
    #errorfile.write(str(it)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, u, norm_type='H1')
    print 'Errornorm H1:', error_norm
    #errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

    # Plot solution and mesh
    if False:
      plot(project(u,bigV), title = 'solution'+str(Nh))
      plot(project(abs(u-u0),bigV), title = 'error'+str(Nh))
      #plot(det(grad(grad(u))), title = 'determinant of hessian')

    Nh = Nh*2
    
    #mesh = refine(mesh)
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma_single = FunctionSpace(mesh, 'DG', deg_hessian)
    Sigma = VectorFunctionSpace(mesh, 'DG', deg_hessian, dim=4)

    
    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    w = Function(V)
    w.assign(project(u,V))
    
    #print 'Errornorm after projection to bigger space :', errornorm(u0, w)
    
  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()