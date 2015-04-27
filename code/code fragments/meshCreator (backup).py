from dolfin import *

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] < DOLFIN_EPS or x[0] > (1.0 - DOLFIN_EPS)) \
                    and on_boundary)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[1] - 1.0
        y[1] = x[0]


def meshFunc(): 
  # Create mesh and finite element
  mesh = UnitSquareMesh(20,20)
  V = FunctionSpace(mesh, "CG", 1, constrained_domain=PeriodicBoundary())
  
  # Create Dirichlet boundary condition
  u0 = Constant(0.0)
  dbc = DirichletBoundary()
  bc0 = DirichletBC(V, u0, dbc)

  # Collect boundary conditions
  bcs = [bc0]
  
  #return functionspace and boundary conditions
  return [V, bcs]
  
