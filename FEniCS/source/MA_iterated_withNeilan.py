from dolfin import *
import scipy.io
from MA_iterated_August import start_iteration
from MA_problem import *
import DG_neilan2 as neilan
import sys, math, time
from DG_iterated import neilan_step

#Calculates the cofactor matrix of the piecewise Hessian
def calc_CofacpiecewiseHessian(mesh, Sigma, u_):

  dh = Function(Sigma)
  #calculate cofac current hessian
  dh = project(as_matrix(((u_.dx(1,1), -u_.dx(0,1)),(-u_.dx(0,1), u_.dx(0,0)))), Sigma)
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
  l_DH = neilan.frobenius_product(grad(grad(u_)),mu)*dx

  #correction term
  l_DH = l_DH - (dot(avg(mu)*nabla_grad(u_)('+'),n('+'))+dot(avg(mu)*nabla_grad(u_)('-'),n('-')))*dS
  
  dh = Function(Sigma)
  solve(a_DH == l_DH, dh)
  
  #calc cofactor matrix  
  cofactorM = Function(Sigma)
  #assign(cofactorM, project(as_matrix(((dh[1,1],-dh[1,0]),(dh[0,1], dh[0,0]))), Sigma))
  assign(cofactorM, project(as_matrix(((dh[1,1],-dh[1,0]-dh[0,1]/2.0),(-dh[1,0]-dh[0,1]/2.0, dh[0,0]))), Sigma))
  return cofactorM


def MA_iteration_withNeilan(mesh, V, Sigma, u0, f, max_it, w, sigmaC, sigmaG, alpha):

  # ------Define variational problem---------
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  for iteration in range(0,max_it):
    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm at beginning:', error_norm

    #cofactor matrix of startsolution's hessian, choose between piecewise evaluation and Neilan's discrete version
    #cofactor = calc_CofacdiscreteHessian(mesh, Sigma, w)
    #cofactor = calc_CofacpiecewiseHessian(mesh, Sigma, w)
    cofactor = cofac(grad(grad(w)))

    #define bilinear form
    a = inner(cofactor*nabla_grad(v) , nabla_grad(u))*dx \
      - v('+')*  dot(avg(cofactor*nabla_grad(u)),n('+'))*dS \
      - v('-')*  dot(avg(cofactor*nabla_grad(u)),n('-'))*dS \
      - u('+')*  dot(avg(cofactor*nabla_grad(v)),n('+'))*dS \
      - u('-')*  dot(avg(cofactor*nabla_grad(v)),n('-'))*dS \
      + Constant(sigmaC)('+')/h('+')* jump(u)*jump(v)*dS \
      + Constant(sigmaG)('+')*h('+')* jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
      - v*inner(n,cofactor*nabla_grad(u))*ds \
      - u*inner(n,cofactor*nabla_grad(v))*ds \
      + Constant(sigmaC)/h*v*u*ds

    #define rhs functional
    L = inner(Constant(-2.0)*f,v)*dx - u0*dot(n,cofactor*nabla_grad(v))*ds +Constant(sigmaC)/h *u0*v*ds

    def u0_boundary(x, on_boundary):
      return on_boundary

    bc = DirichletBC(V, u0, u0_boundary)

    #iterate
    u_ = Function(V)
    solve(a == L, u_, bc)
    
    #plot solution
    if False: 
      plot(project(u, bigV), title = 'solution'+str(Nh)+'-'+str(iteration))
      plot(project(abs(u-u0),bigV), title = 'error'+str(Nh)+'-'+str(iteration))
      interactive()

    #damping and update w
    w.vector()[:] = ((1.0-alpha)*w.vector().array() + alpha*u_.vector().array())

    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm:', error_norm
    errorfile.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, w, norm_type='H1')
    print 'Errornorm H1:', error_norm
    errorfileh1.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')
  return w;

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
  
  Nh = 2
  
  #define mesh and functions spaces 
  mesh = UnitSquareMesh(Nh, Nh, 'crossed')
  V = FunctionSpace(mesh, 'DG', deg)
  #space for discrete hessian
  Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))
  
  bigMesh = refine(mesh)
  bigV = FunctionSpace(bigMesh, 'DG', deg)

  #define poblem
  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
  
    #penalty
  sigmaG = 30.0*deg*deg
  sigmaC = 30.0*deg*deg
  
  #damping
  alpha = 0.5
  
  #start solution
  u_ = start_iteration(mesh, V, u0, f, sigmaC)

  #maximum number of iterations per grid
  max_it = 15
 
  for it in range(1,8):
    
    print "calculating on h=1/", Nh
        
    if False:
      plot(project(u, bigV), title = 'startsolution')
      plot(project(abs(u-u0),bigV), title = 'starterror')
      plot(det(grad(grad(u))), title = 'determinant of starthessian')
      
    w = MA_iteration_withNeilan(mesh, V, Sigma, u0, f, max_it,u_, sigmaC, sigmaG, alpha)
    #w = neilan_step(mesh, V, Sigma, u_, u0, f, sigmaC, sigmaG, sigmaC, max_it, alpha)

    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm:', error_norm
    #errorfile.write(str(it)+' '+str(error_norm)+'\n')

    error_norm = errornorm(u0, w, norm_type='H1')
    print 'Errornorm H1:', error_norm
    #errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

    # Plot solution and mesh

    #------refine----------
    Nh = Nh*2
    
    mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    V = FunctionSpace(mesh, 'DG', deg)
    Sigma = TensorFunctionSpace(mesh, 'DG', deg_hessian, shape=(2,2))

    bigMesh = refine(mesh)
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    u_ = Function(V)
    u_ = project(w,V)
    error_norm = errornorm(u0, w)
    print 'Errornorm after projection:', error_norm
    
    g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh)
    if True:
      #plot(project(u,bigV), title = 'solution'+str(Nh))
      plot(project(abs(u_-u0),bigV), title = 'error'+str(Nh))
      interactive()

    
    print 'Errornorm after projection to bigger space :', errornorm(u0, u_)
    
  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()