from dolfin import *
 
#define mesh and function space
c = Circle(0, 0, 10)
mesh = Mesh(c,10)
mesh2 = Mesh(c,50)
V = FunctionSpace(mesh, 'Lagrange', 3)
V2 = FunctionSpace(mesh2, 'Lagrange', 3)
 
#define boundary
def boundary(x, on_boundary):
    return on_boundary
 
#apply essential boundary conditions
bc = DirichletBC(V, 0, boundary)
 
#define functions
u = TrialFunction(V)
v = TestFunction(V)
 
x2 = Expression('x[0]*x[0]')
 
#define problem
a = (inner(grad(u), grad(v)) \
     + x2*u*v)*dx
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
eigensolver.solve()
 
#extract first eigenpair
r, c, rx, cx = eigensolver.get_eigenpair(0)
print 'eigenvalue:', r
 
#assign eigenvector to function
u = Function(V)
u.vector()[:] = rx
 
#plot eigenfunction and probability density
q = interpolate(u,V2)
plt  = plot(q)
interactive()
 
i = 1
while i < eigensolver.get_number_converged():
    #extract next eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    print 'eigenvalue:', r
 
    #assign eigenvector to function
    u.vector()[:] = rx
 
    #plot eigenfunction
    q = interpolate(u,V2)
    plt.update(q)
    interactive()
 
    #increment i
    i = i+1
