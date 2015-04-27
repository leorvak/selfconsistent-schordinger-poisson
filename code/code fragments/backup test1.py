from dolfin import *
from meshCreator import *
from qho import schrodingerEq
from electronDensity import electronOccupationState, electronDensityFunction
from scipy import constants as pc
import numpy as np
import matplotlib.pyplot as plt 
from poisson import poissonEq
import time

#strat timing
start = time.clock()

#constnats
band_offset = 2.88391914e-11 #joules

meshArray = meshFunc()
V = meshArray[0][1]


potential = interpolate(Expression("x[0]*x[0]/100"), V)

#Start of Iteration
k = 0
while k < 3:
  print str(k) + ' iteration'
  eigenvalues, eigenvectors, u, rx_list = schrodingerEq(potential, meshArray)
  nk = electronOccupationState(eigenvalues)
  n = electronDensityFunction(eigenvectors, nk)
  phi = poissonEq(n,meshArray)



  #get the new potential V(x) = -q phi(x) + band_offset
  new_potential = Function(V)
  new_potential.vector()[:] = np.array([-j + band_offset for j in phi.vector()[:]])

  print '\n'

  #convergence
  print 'convergence'
  v_error = []
  n_error = []

  for i in range(0,len(potential.vector())):
    v_error.append(potential.vector()[i] - new_potential.vector()[i])
    n_error.append(n[i] - new_n[i])
  
  v_error = sum(v_error)/len(v_error)
  n_error = sum(n_error)/len(n_error)

  print 'v_error: ' + str(v_error)
  print 'n_error: ' + str(n_error) + '\n'

  #get the new potential V(x) = -q phi(x) + band_offset
  n = new_n
  potential = new_potentail
#***End of Iteration***

#end timing
end = time.clock()

#print timing information
print 'time elapsed: ' + str(end - start) + ' seconds'

#plotting
v1 = Function(V)
v1.vector()[:] = new_n
Iv1 = interpolate(v1,V)

plot(potential, interactive=True)

plot(Iv1, interactive=True)

plot(new_potential, interactive=True)


