from dolfin import *
import scipy.io
from MA_iterated_August import start_iteration
from MA_problem import *
import DG_neilan2 as neilan
import sys, math, time

#Calculates the cofactor matrix of the piecewise Hessian
def calc_pos_def(A):
  rad = A[0,0] * A[0,0] + (A[1,1] - 2 * A[0,0]) * A[1,1] + 4 * A[0,1] * A[1,0];
  #fetch numerical zeroes
  if math.fabs(rad) < 1e-10:  
    rad = 0;

  s = sqrt(rad);
  ev0 = (A[0,0] + A[1,1] - s) / 0.2e1;
  ev1 = (A[0,0] + A[1,1] + s) / 0.2e1;
 
  if ev0 < eps:
    if ev0 < 0:
      print 'found negative eigenvalue'
      A = A+as_matrix(((-ev0+eps, 0),(0,-ev0+eps)))
    else:
      print 'correcting eigenvalue'
      A = A+as_matrix(((ev0+eps, 0),(0,ev0+eps)))
  
  return A

#Calculates the cofactor matrix of the discrete Hessian (Neilan's approach)
def calc_CofacdiscreteHessian(mesh, Sigma, u_):
  
  #define geometry
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  dh_ = TrialFunctions(Sigma)
  mu_ = TestFunctions(Sigma)
  
  dh = as_tensor(((dh_[0], dh_[1]), (dh_[1], dh_[2])))
  mu = as_tensor(((mu_[0], mu_[1]), (mu_[1], mu_[2])))
  
  #define bilinearform to caculate discrete Hessian
  
  #test with hessian
  a_DH = neilan.frobenius_product(dh,mu)*dx

  #piecewise hessian
  l_DH = neilan.frobenius_product(grad(grad(u_)),mu)*dx

  l_DH = l_DH - dot(avg(mu)*nabla_grad(u_)('+') ,n('+'))*dS \
              - dot(avg(mu)*nabla_grad(u_)('-') ,n('-'))*dS
  
  dh = Function(Sigma)
  solve(a_DH == l_DH, dh)
  
  #calc cofactor matrix and symmetrise
  cofactorM = as_matrix(((dh[2],-dh[1]),(-dh[1], dh[0])))
  #assign(cofactorM, project(as_matrix(((dh[1,1],-dh[1,0]-dh[0,1]/2.0),(-dh[1,0]-dh[0,1]/2.0, dh[0,0]))), Sigma))
  
  #plot(cofactorM[0,1], title = 'numer derivative01 '+str(Nh))
  #plot(cofactorM[1,0], title = 'numer derivative10 '+str(Nh))
  #interactive()
  
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
      #- avg(v) * jump(cofactor*nabla_grad(u),n)*dS \

 

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
  #damping
  alpha = 0.3

  #parameters['form_compiler']['quadrature_rule'] = 'canonical'
  parameters['form_compiler']['quadrature_degree'] = 2*deg

  fileprefix = problem_name+'_iterated_deg'+str(deg)+'_alpha'+str(alpha)
  print "processing files ", fileprefix
  
  errorfile = open('data/'+fileprefix+'l2errornorm','wa',1)
  errorfile.write('iterations l2error\n');
  fitdata_file = open('data/'+fileprefix+'fitData','wa',1)
  errorfileh1 = open('data/'+fileprefix+'h1errornorm','wa',1)
  errorfileh1.write('iterations h1error\n');
  
  Nh = 2
  
  #define mesh and functions spaces 
  init_mesh = UnitSquareMesh(1, 1, 'crossed')
  mesh = []
  mesh.append( adapt(init_mesh))
  V = FunctionSpace(mesh[0], 'DG', deg)
  #space for discrete hessian
  Sigma = VectorFunctionSpace(mesh[0], 'DG', deg_hessian, dim=3)
  
  bigMesh = refine(mesh[0])
  bigV = FunctionSpace(bigMesh, 'DG', deg)

  #define poblem
  g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh[0])
  
    #penalty
  sigmaG = 50.0
  sigmaC = 20.0*deg*deg
  
  
  #start solution
  u_=[]
  u_.append(start_iteration(mesh[0], V, u0, f, sigmaC))

  #maximum number of iterations per grid
  max_it = 15
 
  for it in range(1,8):
    
    print "calculating on h=1/", Nh
        
    if False:
      plot(project(u, bigV), title = 'startsolution')
      plot(project(abs(u-u0),bigV), title = 'starterror')
      plot(det(grad(grad(u))), title = 'determinant of starthessian')
      
    w = MA_iteration_withNeilan(mesh[it-1], V, Sigma, u0, f, max_it,u_[it-1], sigmaC, sigmaG, alpha)
    #w = neilan_step(mesh, V, Sigma, u_, u0, f, sigmaC, sigmaG, sigmaC, max_it, alpha)
    
    #examine error
    #error_norm = errornorm(u0, w)
    #print 'Errornorm:', error_norm
    #errorfile.write(str(it)+' '+str(error_norm)+'\n')

    #error_norm = errornorm(u0, w, norm_type='H1')
    #print 'Errornorm H1:', error_norm
    #errorfileh1.write(str(it)+' '+str(error_norm)+'\n')

    # Plot solution and mesh

    #------refine----------
    Nh = Nh*2
    
    #plot(mesh)
    
    #mesh = UnitSquareMesh(Nh, Nh, 'crossed')
    mesh.append(adapt(mesh[it-1]))
    #plot(mesh)
    V = FunctionSpace(mesh[it], 'DG', deg)
    Sigma = VectorFunctionSpace(mesh[it], 'DG', deg_hessian, dim=3)

    bigMesh = refine(mesh[it])
    bigV = FunctionSpace(bigMesh, 'DG', deg)
    
    g, f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh[it])

    #u_ = project(w,V)
    u_.append(Function(adapt(w, mesh[it])))
    error_norm = errornorm(u0, w)
    print 'Errornorm after projection:', error_norm
    
    if False:
      #plot(project(u,bigV), title = 'solution'+str(Nh))
      #plot(project(abs(u_-u0),bigV), title = 'error'+str(Nh))
      interactive()
    print 'Errornorm after projection to bigger space :', errornorm(u0, u_[it])
    
  end = time.clock()
  time_file = open('data/timing','a')
  time_file.write(fileprefix+' '+str(end-start)+'\n')
  
  print "%.2gs" % (end-start)
  
  errorfile.close()
  errorfileh1.close()
  time_file.close()
