from dolfin import *
from meshCreator import *
import numpy as np
from scipy import constants as pc

#constants
e_g = Constant(12.9)
e_a = Constant(12.9-2.84*0.3)

#create mesh
nx = 20;  ny = 20
mesh = UnitSquareMesh(nx,ny)

#defining the scale of substrates
total_width = 15.0+20.0+5.0+500.0
gaas_thickness1 = 15.0/total_width
algaas_thickness1 = 20.0/total_width
algaas_thickness2 = 5.0/total_width
gaas_thickness2 = 500.0/total_width

def poissonEq(electronDensity, meshArray):
  
  # Define boundary condition
  # Create classes for defining parts of the boundaries and the interior
  # of the domain
  class Bottom(SubDomain):
      def inside(self, x, on_boundary):
          return near(x[0], 0.0)

  class Top(SubDomain):
      def inside(self, x, on_boundary):
          return near(x[0], 1.0)

  class AlGaAs1(SubDomain):
      def inside(self, x, on_boundary):
          return (between(x[0], (gaas_thickness1, gaas_thickness1 + algaas_thickness1)))
          
  class AlGaAs2(SubDomain):
      def inside(self, x, on_boundary):
          return (between(x[0], (gaas_thickness1 + algaas_thickness1, gaas_thickness1 + algaas_thickness1 + algaas_thickness2)))
  
  # Sub domain for Periodic boundary condition
  class PeriodicBoundary(SubDomain):
  
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[1] - 1.0
        y[1] = x[0]

  # Initialize sub-domain instances
  top = Top()
  bottom = Bottom()
  algaas1 = AlGaAs1()
  algaas2 = AlGaAs2()
  
  
  # Initialize mesh function for interior domains
  domains = CellFunction("size_t", mesh)
  domains.set_all(0)
  algaas1.mark(domains, 1)
  algaas2.mark(domains, 2)
  
  # Initialize mesh function for boundary domains
  boundaries = FacetFunction("size_t", mesh)
  boundaries.set_all(0)
  top.mark(boundaries, 1)
  bottom.mark(boundaries, 2)
  
  #Define Function space
  V =FunctionSpace(mesh, "Lagrange", 1, constrained_domain=PeriodicBoundary())
  
  
  #u0 = Constant('0') 
  #bc = DirichletBC(V, u0, boundary)
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  n = Function(V)
  
  #computing the source term -q/epsilon0 (Nd(x) - n(x))
  n.vector()[:] = np.array([i for i in electronDensity])
  
  # Define Dirichlet boundary conditions at top and bottom boundaries
  bcs = [DirichletBC(V, 0.0, boundaries, 1), DirichletBC(V, 1e18, boundaries, 2)]
       
  # Define new measures associated with the interior domains and
  # exterior boundaries
  dx = Measure("dx")[domains]
  ds = Measure("ds")[boundaries]
  
  # Define variational form
  a = inner(e_g*grad(u), grad(v))*dx(0) + inner(e_a*grad(u), grad(v))*dx(1) + inner(e_a*grad(u), grad(v))*dx(2)
  L =   n*v*dx(0) + n*v*dx(1)
  
  # Compute solution
  u = Function(V)
  print 'solving Poisson\'s equation...'
  solve(a == L, u, bcs)
  #plot(u,interactive=True, axes=True)
  
  return u
