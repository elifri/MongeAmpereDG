from dolfin import *
import math

def norm(x):
  return sqrt(x[0]**2 + x[1]**2)

#-----------choose PDE------------
def MA_problem(name, Nh, degree, mesh):
  
  #define rhs
  class rhs(Expression):    
    def __init__(self, mesh, **kwargs):
        self._mesh = mesh
    def eval(self, v, x):
      if name == 'MA2':
        if math.fabs(x[0] - 1) <= 1e-12 and math.fabs(x[1] - 1) <= 1e-12:
          v[0] = 0
          print 'point', x
        else:
          val = 2-x[0]**2-x[1]**2
          v[0]= 2/(val**2)
      elif name =='MA3':
        if math.fabs(x[0] - 0.5) <= 1e-12 and math.fabs(x[1] - 0.5) <= 1e-12:
          v[0] = 0
          print 'point', x
        else:
          val = sqrt ((x[0]-0.5)**2+(x[1]-0.5)**2)
          val = 1 - 0.2/val
          if val > 0:
            v[0] = val
          else:
            v[0] = 0
      elif name =='MA4':
        if math.fabs(x[0] - 0.5) <= 1./Nh and math.fabs(x[1] - 0.5) <= 1./Nh:
          v[0] = pi/4.0*Nh*Nh
        else:
          v[0]=0
  
  class exact_sol(Expression):
    def __init__(self, mesh, **kwargs):
        self._mesh = mesh
    def eval(self, v, x):
      if name =='MA3':
        val = sqrt ( (x[0]-0.5)**2 + (x[1]-0.5)**2 )
        val = val - 0.2
        if val > 0:
          v[0] = 0.5*val**2
        else:
          v[0] = 0
  
  if name=='MA1':
    f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])', element = FiniteElement("Quadrature", triangle, degree))#MongeAmpere1
  elif name == 'MA1Reverse':
    f = Expression('(1 + (1-x[0])*(1-x[0])+x[1]*x[1]) * exp((1-x[0])*(1-x[0])+x[1]*x[1])', element = FiniteElement("Quadrature", triangle,degree))
  elif name == 'simpleMA':
    f = Constant(7.0)#simpleMongeAmpere
  elif name == 'simpleMA2':
    f = Constant(1.0) #simpleMongeAmpere2
  elif name == 'Brenner1':
    f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)', element = FiniteElement("Quadrature", triangle, degree))#BrennerEx1
  elif name == 'MA2':
    #f = rhs(element = FiniteElement("Quadrature", triangle, degree))
    element = FiniteElement("Quadrature", triangle, degree)
    f = rhs(mesh, element=element)
  elif name == 'MA3':
    f = rhs(mesh, element = FiniteElement("Quadrature", triangle, degree))
  elif name == 'MA4':
    f = rhs(mesh, element = FiniteElement("Quadrature", triangle, degree))
  elif name == 'const rhs':
    f = Constant(1.0)
   
  # Define boundary conditions
  if name=='MA1':
    g = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  elif name=='MA1Reverse':
    g = Expression('exp( (pow((1-x[0]),2)+pow(x[1],2))/2. )')#MongeAmpere1
  elif name == 'simpleMA':
    g = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  elif name == 'simpleMA2':
    g = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  elif name == 'Brenner1':
    g = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1
  elif name == 'MA2':
    g = Expression('-sqrt(2-pow(x[0],2)-pow(x[1],2))') #test 2
  elif name == 'MA3':
    g  = exact_sol(mesh)
  elif name == 'MA4':
    g  = Expression('sqrt( pow(x[0]-0.5,2) + pow(x[1]-0.5,2))') #test 4
  elif name == 'const rhs':
    g = Expression('0.0') #const rhs
  return f, g
  
def start_iteration(mesh, V, u0, f, sigma):
  
  #define variables

  #define exact solution
  u_e = interpolate(u0, V)
  
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)

  #define geometry for penalty and normals
  h = CellSize(mesh)
  n = FacetNormal(mesh)

  #define bilinear form
  a = inner(nabla_grad(v), nabla_grad(u))*dx \
    - inner(jump(v,n),avg(nabla_grad(u)))*dS\
    - inner(jump(u,n),avg(nabla_grad(v)))*dS\
    + Constant(sigma)('+')/h('+')* jump(u)*jump(v)*dS \
    + Constant(sigma/2.0)('+')/h('+')*jump(nabla_grad(u),n)*jump(nabla_grad(v),n)*dS \
    - v*inner(n,nabla_grad(u))*ds \
    - u*inner(n,nabla_grad(v))*ds \
    + Constant(sigma)/h*v*u*ds

  #define rhs functional
  rhs = 2*f 
  L = inner(-sqrt(rhs),v)*dx - u0*dot(n,nabla_grad(v))*ds +Constant(sigma)/h *u0*v*ds


  #solve Poisson problem
  u = Function(V)

  solve(a == L, u)
  
  #examine error
  print 'Errornorm:', errornorm(u0,u)
  
  return u;





