from dolfin import *
from meshCreator2 import *
from vschrodinger import schrodingerEq
from velectronDensity import electronOccupationState, electronDensityFunction
from scipy import constants as pc
import numpy as np
import matplotlib.pyplot as plt 
from vpoisson import poissonEq
from output import save_output
import time

#strat timing
start = time.clock()

#constnats
band_offset = 0.23 #MeV

meshArray = meshFunc()
V = meshArray[0]


potential = interpolate(Expression("1"), V)
n = interpolate(Expression("1"), V)

#Start of Iteration
k = 0
while k < 5:
  print str(k+1) + ' iteration'
  eigenvalues, eigenvectors, u, rx_list = schrodingerEq(potential, meshArray)
  nk = electronOccupationState(eigenvalues)
  new_n = electronDensityFunction(eigenvectors, nk)
  phi = poissonEq(new_n,meshArray)
  


  #get the new potential V(x) = -q phi(x) + band_offset
  new_potential = Function(V)
  new_potential.vector()[:] = np.array([-j + band_offset for j in phi.vector()[:]])
  

  print '\n'
  
  
  #convergence
  print 'convergence'
  v_error = []
  n_error = []

  for i in range(len(potential.vector())):
    v_error.append((new_potential.vector()[i]-potential.vector()[i])/potential.vector()[i])
    n_error.append((n.vector()[i]-new_n[i])/n.vector()[i])
  
  v_error = abs(sum(v_error)/len(v_error))
  n_error = abs(sum(n_error)/len(n_error))

  print 'v_error: ' + str(v_error)
  print 'n_error: ' + str(n_error) + '\n'


  #update potential, electron density and iteration
  potential = new_potential
  n.vector()[:] = np.array([j for j in new_n])
  plot(n)
  k = k+1
  
  if(n_error < 10e-5 and v_error < 10e-5):
    print 'Convergence occured at k: ' + str(k)
    n_average = sum(new_n)/len(new_n)
    print 'average electron density = ' + str(n_average)
    break
  
#***End of Iteration***

#end timing
end = time.clock()

#print timing information
print 'time elapsed: ' + str(end - start) + ' seconds'

#create functions for eigenvectors
#eigenfunctions = []
#i = 0
#while(i < len(eigenvectors)):
#  u = Function(V)
#  u.vector()[:] = np.array(eigenvectors[i])
#  Iu = interpolate(u,V)
#  eigenfunctions.append(Iu)
#  i = i + 1

#save output in files  
#save_output(eigenfunctions, phi, n, potential)

#plotting
v1 = Function(V)
v1.vector()[:] = new_n
Iv1 = interpolate(v1,V)


plot(Iv1, interactive=True, axes=True, title='Electron Density n(x)')

plot(new_potential, interactive=True, axes=True, title='V(x)')

plot(phi, interactive=True, axes=True, title='phi(x)')



