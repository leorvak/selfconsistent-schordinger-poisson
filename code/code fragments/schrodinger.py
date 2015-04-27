from dolfin import *
from meshCreator import *
from scipy import constants as pc
from scipy import sqrt
import numpy as np
import math

#constants
bulk_k = 7.53e10 #N/m^2
hb2m_g = Constant(3.81/0.063) #Mev*s
hb2m_a = Constant(3.81/0.0879)
energy_gap = 5.0 #mev

#create mesh
nx = 20;  ny = 20
mesh = UnitSquareMesh(nx,ny)

#defining the scale of substrates
total_width = 15.0+20.0+5.0+500.0
gaas_thickness1 = 15.0/total_width
algaas_thickness1 = 20.0/total_width
algaas_thickness2 = 5.0/total_width
gaas_thickness2 = 500.0/total_width

def schrodingerEq(potential, meshArray):  
 
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
 
 
 #define functions
 u = TrialFunction(V)
 v = TestFunction(V)
 
 #boundary conditions
 bcs = [DirichletBC(V, 0.0, boundaries, 1), DirichletBC(V, 0.0, boundaries, 2)]
 
 # Define new measures associated with the interior domains and
 # exterior boundaries
 dx = Measure("dx")[domains]
 ds = Measure("ds")[boundaries]
  
 #define problem
 a = (inner(grad(u), grad(v)) \
      +  hb2m_g * potential*u*v)*dx(0) + (inner(grad(u), grad(v)) \
      +  hb2m_a * potential*u*v)*dx(1)
 m = hb2m_g *u*v*dx(0) + hb2m_a *u*v*dx(1)
 
 #assemble stiffness matrix
 #A = PETScMatrix()
 #assemble(a, A_tensor=A)
 #M = PETScMatrix()
 #assemble(m, A_tensor=M)
 
 A = PETScMatrix() 
 M = PETScMatrix()
 _ = PETScVector()
 L = Constant(0.)*v*dx(0) + Constant(0.)*v*dx(1) + Constant(0.)*v*dx(2) + Constant(0.)*v*dx(3)
 #assemble(a, tensor=A)
 #assemble(m, tensor=M)
 assemble_system(a, L, bcs, A_tensor=A, b_tensor=_)
 #assemble_system(m, L,bc, A_tensor=M, b_tensor=_)
 assemble_system(m, L, A_tensor=M, b_tensor=_)
 
 #create eigensolver
 eigensolver = SLEPcEigenSolver(A,M)
 eigensolver.parameters['spectrum'] = 'smallest magnitude'
 eigensolver.parameters['solver']   = 'lapack'
 eigensolver.parameters['tolerance'] = 1.e-15
 
 #solve for eigenvalues
 print 'solving Schrodinger\'s equation...'
 eigensolver.solve()
 
 
 #assign eigenvector to function
 u = Function(V)
 eigenvectors = []
 eigenvalues = []
 rx_list = []

 
 #extract first eigenpair
 r, c, rx, cx = eigensolver.get_eigenpair(0)
 
 #assign eigenvector to function
 #rx = np.array([q*xsi_factor for q in rx])
 u.vector()[:] = rx
 eigenvalues.append(r)
 eigenvectors.append([j for j in u.vector()])
     
 rx_list.append(rx) 
 f = Function(V)
 f.interpolate(potential)
 #plot(f, interactive=True)
 #plot eigenfunction and probability density
 print 'eigenvalue: ' + str(r)
 #plt  = plot(u)
 #interactive()
 
 i = 1
 while i < eigensolver.get_number_converged():
     #extract next eigenpair
     r, c, rx, cx = eigensolver.get_eigenpair(i)
     #assign eigenvector to function
     #rx = np.array([q*xsi_factor for q in rx])
     u.vector()[:] = rx
     eigenvalues.append(r)
     eigenvectors.append([j for j in u.vector()])
     
     rx_list.append(rx) 
     
     #plot eigenfunction
     print 'eigenvalue: ' + str(r)
     #plot(u).update(u)
     #interactive()
  
     #increment i
     i = i+1
 return [eigenvalues, eigenvectors,u,rx_list]
