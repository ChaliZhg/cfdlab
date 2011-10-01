"""
This demo program solves the steady incompressible Navier-Stokes equations
for cylinder in channel problem using Taylor-Hood elements.
   Author: Praveen. C
   www   : http://math.tifrbng.res.in/~praveen
"""

# Set parameter values
Re   = 80
D    = 0.1
Uinf = 1.0
nu   = D * Uinf / Re

from dolfin import *

# Load mesh from file
mesh = Mesh("cylinder_in_channel.xml")
sub_domains = MeshFunction("uint", mesh, "subdomains.xml")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
w     = Function(W)
(u,p) = (as_vector((w[0], w[1])), w[2])

# Define boundary conditions
uinlet = Expression(("(1.0 - (x[1]/0.2)*(x[1]/0.2))", "0"))
noslip = DirichletBC(W.sub(0), (0, 0), sub_domains, 0)
inlet  = DirichletBC(W.sub(0), uinlet, sub_domains, 1)
bc     = [noslip, inlet]

# Stress tensor
T = nu*(grad(u) + grad(u).T) - p*Identity(2)
# Face normals
n = FacetNormal(mesh)

# Weak form
F =   inner(grad(u)*u, v)*dx \
    + inner(T, grad(v))*dx   \
    - q*div(u)*dx

# Derivative of weak form
dw = TrialFunction(W)
dF = derivative(F, w, dw)

pde = VariationalProblem(F, dF, bc)
#info(pde.parameters, True)
#quit()
# Set linear solver parameters
itsolver = pde.parameters["solver"]["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-10
itsolver["relative_tolerance"] = 1.0e-6
pde.solve(w)

# Save steady solution
File("steady.xml") << w.vector()

# Save vtk for visualization
(u,p) = w.split()
File("velocity.pvd") << u
File("pressure.pvd") << p

# Compute and save vorticity in vtk format
r = TrialFunction(Q)
s = TestFunction(Q)
a = r*s*dx
L = (u[0].dx(1) - u[1].dx(0))*s*dx
vv= VariationalProblem(a, L)
vort = vv.solve()
File("vorticity.pvd") << vort
