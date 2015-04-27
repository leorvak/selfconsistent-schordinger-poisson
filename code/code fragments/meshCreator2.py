from dolfin import *
import sys, numpy, math

#**********************************************************
nx = 2*250 + 80;
mesh = UnitIntervalMesh(nx)

#defining the scale of substrates
total_width = nx
gaas_thickness1 = 250/total_width
algaas_thickness1 = (gaas_thickness1+80)/total_width
gaas_thickness2 = (gaas_thickness1 + algaas_thickness1)/total_width

#mass of electron in MeV/c^2
m0 = 0.511 #MeV/c^2


# Define a MeshFunction over two subdomains
subdomains = MeshFunction('size_t', mesh, 1)

class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[0] <= gaas_thickness1 else False

class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0] >= gaas_thickness1 and x[0] <= algaas_thickness1) else False
        
class Omega2(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0] >= algaas_thickness1 and x[0] <= gaas_thickness2) else False

# note: it is essential to use <= and >= in the comparisons

# Mark subdomains with numbers 0 and 1
subdomain0 = Omega0()
subdomain0.mark(subdomains, 0)
subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)
subdomain2 = Omega2()
subdomain2.mark(subdomains, 2)

V0 = FunctionSpace(mesh, 'CG', 3)
k = Function(V0)
nd = Function(V0)
m_eff = Function(V0)


# Loop over all cell numbers, find corresponding
# subdomain number and fill cell value in k
k_values = [12.9, (12.9-2.84*0.3), 12.9]  # values of k in the two subdomains
nd_values = [180.0, 0.0, 0.0]  # values of nd in the two subdomains
m_eff_values = [0.063*m0, 0.092*m0, 0.063*m0]

for cell_no in range(len(subdomains.array())):
    subdomain_no = subdomains.array()[cell_no]
    k.vector()[cell_no] = k_values[subdomain_no]
    nd.vector()[cell_no] = nd_values[subdomain_no]
    m_eff.vector()[cell_no] = m_eff_values[subdomain_no]
    
    
m_eff.vector()[:800] = m_eff_values[0]
m_eff.vector()[800:] = m_eff_values[1]

nd.vector()[:800] = nd_values[0]
nd.vector()[800:] = nd_values[1]

k.vector()[:800] = k_values[0]
k.vector()[800:] = k_values[1]


plot(m_eff, title='m_eff', interactive=True)
plot(k, title='e_s', interactive=True)
plot(nd, title='nd', interactive=True)
#***************************************************************

#NEW!!!!!!!!!!!
# Thickness of the well (in Angstroms)
dwell = 80
 
# Thickness of the barriers 
dbar = 250

# Total thickess of structure
dtot = 2*dbar + dwell
# Boundaries
class BoundaryPoint(SubDomain):
    def __init__(self,bpt):
        SubDomain.__init__(self)
        self.bpt = bpt
    def inside(self,x,on_boundary):
        return near(x[0],self.bpt)
# Layers are an interval SubDomain 
class Layer(SubDomain):
    def __init__(self,ends):
        SubDomain.__init__(self)
        self.ends = ends
       
    def inside(self,x,on_boundary):
        return between(x[0],self.ends)
        
# Initialize list of subdomain instances
layers = []
layers.append(Layer((0,dbar)))
layers.append(Layer((dbar,dbar+dwell)))
layers.append(Layer((dbar+dwell,dtot)))

# Initialize boundary instances
left = BoundaryPoint(0)
right = BoundaryPoint(dtot)

#END NEW!!!

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((x[0] < DOLFIN_EPS or x[0] > (1.0 - DOLFIN_EPS)) \
                    )

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - 1.0


def meshFunc(): 
  boundaries = FacetFunction("size_t",mesh)
  left.mark(boundaries,0)
  right.mark(boundaries,1)
  V = FunctionSpace(mesh,"CG",3)
  bcs1 = [DirichletBC(V,0.0,left),DirichletBC(V,0.0,right)]

  # Create mesh and finite element
  #V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=PeriodicBoundary())
  
  # Create Dirichlet boundary condition
  u0 = Constant(0.0)
  dbc = DirichletBoundary()
  bc0 = DirichletBC(V, u0, dbc)

  # Collect boundary conditions
  bcs = [bc0]
  
  #return functionspace and boundary conditions
  return [V, bcs1, k, nd, m_eff]
  
