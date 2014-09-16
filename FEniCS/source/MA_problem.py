from dolfin import *
import math

#-----------choose PDE------------
def MA_problem(name):
  
  #define rhs
  class rhs(Expression):
    def eval(self, v, x):
      if name == 'MA2':
        if math.fabs(x[0] - 1) <= 1e-12 and math.fabs(x[1] - 1) <= 1e-12:
          v[0] = Nh
        else:
          val = 2-x[0]**2-x[1]**2
          v[0]= 2/(val**2)
  
  if name=='MA1':
    f = Expression('(1 + x[0]*x[0]+x[1]*x[1]) * exp(x[0]*x[0]+x[1]*x[1])')#MongeAmpere1
  elif name == 'simpleMA':
    f = Constant(7.0)#simpleMongeAmpere
  elif name == 'simpleMA2':
    f = Constant(1.0) #simpleMongeAmpere2
  elif name == 'Brenner1':
    f = Expression('2000*pow(exp(pow(x[0],6)/6+x[1]),2)*pow(x[0],4)')#BrennerEx1
  elif name == 'MA2':
    f = rhs()
  elif name == 'MA3':
    f = rhs()
  
  # Define boundary conditions
  if name=='MA1':
    u0 = Expression('exp( (pow(x[0],2)+pow(x[1],2))/2. )')#MongeAmpere1
  elif name == 'simpleMA':
    u0 = Expression('2*x[0]*x[0] + 2*x[1]*x[1] + 3*x[0]*x[1]') #simpleMongeAmpere
  elif name == 'simpleMA2':
    u0 = Expression('x[0]*x[0]/2.0 + x[1]*x[1]/2.0') #simpleMongeAmpere2
  elif name == 'Brenner1':
    u0 = Expression('20*exp(pow(x[0],6)/6.0+x[1])')#BrennerEx1
  elif name == 'MA2':
    u0 = Expression('-sqrt(2-pow(x[0],2)-pow(x[1],2))') #test 2
  elif name == 'MA3':
    u0  = Expression('sqrt( pow(x[0]-0.5,2) + pow(x[1]-0.5,2))') #test 4
  #u0 = Constant(0.0) #const rhs

  return u0, f