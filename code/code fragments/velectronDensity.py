import mpmath
from scipy import inf
from sympy import *
import numpy as np
import math
from scipy import constants as pc

#constants
fermiEnergy = 1e9 #for density of 2.8e15
m_eff = 0.063 * 0.511 #MeV
hb2m = 3.81
nk_factor = 1/2/hb2m/pc.pi
boltzmann = 8.6173324e-11 #mev
T = 300

#function to return the electron occupation state nk
def electronOccupationState(eigenvalues):
 print 'Finding n_k...'
 nk = []
 x = Symbol('x')
 for i in range(0,len(eigenvalues)):
  ek = float(eigenvalues[i])
  #result = nk_factor *( limit(x-ln(1+exp(x-fermiEnergy)/pc.k/300),x,oo) - (ek - ln(1+exp(ek-fermiEnergy)/pc.k/300)))
  result = mpmath.quad(lambda x: nk_factor/(1+mpmath.exp((x-fermiEnergy))), [ek, mpmath.inf])
  print float(result)
  nk.append(float(result))

 return nk

#fuction to calculate the electro density
def electronDensityFunction(eigenvectors,  nk):
 print 'Finding the electron density function n(x)'
 result = []
 for i in range(len(eigenvectors)):
   kth_term = [(j * j * nk[i]) for j in eigenvectors[i]]
   result.append(kth_term)
   
 n = np.sum(result, axis=0)
 return n
