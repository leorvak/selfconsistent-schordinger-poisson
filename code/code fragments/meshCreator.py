from dolfin import *
import sys, numpy, math

#**********************************************************
nx = 20;  ny = 20
mesh = UnitSquareMesh(nx,ny)

#defining the scale of substrates
total_width = 15.0+20.0+5.0+500.0
gaas_thickness1 = 15.0/total_width
algaas_thickness1 = 20.0/total_width
algaas_thickness2 = 5.0/total_width
gaas_thickness2 = 500.0/total_width

#mass of electron in MeV/c^2
m0 = 0.511 #MeV/c^2


# Define a MeshFunction over two subdomains
subdomains = MeshFunction('size_t', mesh, 2)

class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] <= gaas_thickness1 else False

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0] >= gaas_thickness1 and x[0] <= algaas_thickness1) else False
        
class Omega2(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0] >= algaas_thickness1 and x[0] <= algaas_thickness2) else False

class Omega3(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0] >= algaas_thickness2 and x[0] <= gaas_thickness2) else False
# note: it is essential to use <= and >= in the comparisons

# Mark subdomains with numbers 0 and 1
subdomain0 = Omega0()
subdomain0.mark(subdomains, 0)
subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)
subdomain2 = Omega2()
subdomain2.mark(subdomains, 2)
subdomain3 = Omega3()
subdomain3.mark(subdomains, 3)

V0 = FunctionSpace(mesh, 'DG', 0)
k = Function(V0)
nd = Function(V0)
m_eff = Function(V0)

#plot(subdomains, interactive=True)

# Loop over all cell numbers, find corresponding
# subdomain number and fill cell value in k
k_values = [12.9, (12.9-2.84*0.3), (12.9-2.84*0.3), 12.9]  # values of k in the two subdomains
nd_values = [1e16, 1e16, 0.0, 0.0]  # values of nd in the two subdomains
m_eff_values = [0.063*m0, 0.0879*m0, 0.0879*m0, 0.063*m0]
for cell_no in range(len(subdomains.array())):
    subdomain_no = subdomains.array()[cell_no]
    k.vector()[cell_no] = k_values[subdomain_no]
    nd.vector()[cell_no] = nd_values[subdomain_no]
    m_eff.vector()[cell_no] = m_eff_values[subdomain_no]

plot(k, title='dielectric', axes=True)
plot(m_eff, title='effective mass', axes=True)
plot(nd , title='doping', axes = True)

#plot(subdomains, title='subdomains')
#***************************************************************

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] < DOLFIN_EPS or x[0] > (1.0 - DOLFIN_EPS)) \
                    )

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[1] - 1.0
        y[1] = x[0]


def meshFunc(): 
  # Create mesh and finite element
  V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=PeriodicBoundary())
  
  # Create Dirichlet boundary condition
  u0 = Constant(0.0)
  dbc = DirichletBoundary()
  bc0 = DirichletBC(V, u0, dbc)

  # Collect boundary conditions
  bcs = [bc0]
  
  #return functionspace and boundary conditions
  return [V, bcs, k, nd, m_eff, mesh]
  
