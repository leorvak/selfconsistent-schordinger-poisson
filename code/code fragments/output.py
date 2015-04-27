from dolfin import HDF5File

#save output to files
def save_output(eigenvectors, phi, e_density, potential):
  
  #writing to HDF file
  f = HDF5File("out.h5", "w")
  
  f.write(phi, 'phi')
  f.write(e_density, 'n(x)')
  f.write(potential, 'V(x)')
  
  i = 0
  while (i < len(eigenvectors)):
    name = 'eigenvector'+str(i)
    print eigenvectors[i]
    print phi
    print e_density
    #f.write(eigenvectors[i], name)
