from dolfin import *
from meshCreator2 import *
from scipy import constants as pc
from scipy import sqrt
import numpy as np
import math

#constants
bulk_k = 7.53e10 #N/m^2
hb2m =3.81 #Mev*s
#omega = sqrt(bulk_k/m_eff)
#xsi_factor = sqrt(pc.hbar/omega/m_eff)
#e_factor = pc.hbar * omega/2


def schrodingerEq(potential, meshArray):  
 #define mesh and function space 
 V = meshArray[0]
 
 #effective mass
 m_eff = meshArray[4]
 
 
 #define functions
 u = TrialFunction(V)
 v = TestFunction(V)
  
 #define problem
 a = (inner(grad(u), grad(v)) \
      + potential*u*v)*dx
 m = u*v*dx
 
 #assemble stiffness matrix
 A = PETScMatrix()
 assemble(a, tensor=A)
 M = PETScMatrix()
 assemble(m, tensor=M)
 
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
     #if(i < 3):
     #  plot(u).update(u)
     #  interactive()
  
     #increment i
     i = i+1
 return [eigenvalues, eigenvectors,u,rx_list]
