from dolfin import *
import scipy.io
from MA_problem import *
import DG_neilan as neilan
import sys, math, time

#Calculates the cofactor matrix of the discrete Hessian (Neilan's approach)
def calc_CofacdiscreteHessian(mesh, Sigma, u_): 
  #define geometry
  h = CellSize(mesh)
  n = FacetNormal(mesh)
  
  dh_ = TrialFunctions(Sigma)
  mu_ = TestFunctions(Sigma)
  
  dh = as_tensor(((dh_[0], dh_[1]), (dh_[1], dh_[2])))
  mu = as_tensor(((mu_[0], mu_[1]), (mu_[1], mu_[2])))
  
  #------define bilinearform to calculate discrete Hessian------
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
  return cofactorM


def MA_iteration_withNeilan(mesh, V, Sigma, u0, f, max_it, w, sigmaC, sigmaG, alpha):

  # ------Define variational problem---------
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  for iteration in range(0,max_it):
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

    #define rhs functional
    L = inner(Constant(-2.0)*f,v)*dx - u0*dot(n,cofactor*nabla_grad(v))*ds +Constant(sigmaC)/h *u0*v*ds

    def u0_boundary(x, on_boundary):
      return on_boundary

    bc = DirichletBC(V, u0, u0_boundary)

    #iterate
    u_ = Function(V)
    solve(a == L, u_, bc)

    #damping and update w
    w.vector()[:] = ((1.0-alpha)*w.vector().array() + alpha*u_.vector().array())

    #examine error
    error_norm = errornorm(u0, w)
    print 'Errornorm:', error_norm
    errorfile.write(str((it-1)*max_it+iteration)+' '+str(error_norm)+'\n')
  return w;

if __name__ == "__main__":
  #------read parameters--------
  if len(sys.argv) != 4:
    print 'Error, please specify the problem, the polynomial degrees of the trial and the Hessian trial fcts!'
    sys.exit(-1)
  
  problem_name = sys.argv[1]
  
  deg = int(sys.argv[2])
  deg_hessian = int(sys.argv[3])
  
  #damping
  alpha = 0.3

  parameters['form_compiler']['quadrature_degree'] = 2*deg

  Nh = 2
  
  #define mesh and functions spaces 
  init_mesh = UnitSquareMesh(1, 1, 'crossed')
  mesh = []
  mesh.append( adapt(init_mesh))
  V = FunctionSpace(mesh[0], 'DG', deg)
  #space for discrete hessian, its a 3-dimensional vector representing a symmetric matrix
  # first component is a_00, second entry is a_10=a_01, third entry is a_11
  Sigma = VectorFunctionSpace(mesh[0], 'DG', deg_hessian, dim=3)

  #define poblem
  f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh[0])
  
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
      
    w = MA_iteration_withNeilan(mesh[it-1], V, Sigma, u0, f, max_it,u_[it-1], sigmaC, sigmaG, alpha)
    #------refine----------
    Nh = Nh*2
    mesh.append(adapt(mesh[it-1]))

    V = FunctionSpace(mesh[it], 'DG', deg)
    Sigma = VectorFunctionSpace(mesh[it], 'DG', deg_hessian, dim=3)
    
    #update problem data
    f, u0 = MA_problem(problem_name, Nh, parameters['form_compiler']['quadrature_degree'], mesh[it])

	#project current solution to refined space
    u_.append(Function(adapt(w, mesh[it])))
    
