from dolfin import *
from meshCreator2 import *
import numpy as np
from scipy import constants as pc

#constants
e_0 = 8.85418782e-12 #vacuum permitivity

def poissonEq(electronDensity, meshArray):
  #Define Function space
  V = meshArray[0]
  
  #boundary conditions
  bcs = meshArray[1]
  
  #dielectric constant
  e_s =meshArray[2]
  
  #dopant
  nd = meshArray[3]
  
  
  # Define boundary condition
  #u0 = Constant('0') 
  #bc = DirichletBC(V, u0, boundary)
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  n = Function(V)
  
  #computing the source term -q/epsilon0 (Nd(x) - n(x))
  n.vector()[:] = np.array([i for i in electronDensity])
  
  a = e_s * inner(grad(u), grad(v))*dx
  L = (n)*v*dx

  # Compute solution
  u = Function(V)
  print 'solving Poisson\'s equation...'
  solve(a == L, u, bcs)
  #plot(u,interactive=True, axes=True)
  
  return u
