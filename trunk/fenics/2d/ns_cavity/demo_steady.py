"""
This demo program solves the steady incompressible Navier-Stokes equations
for lid-driven cavity problem using Taylor-Hood elements
"""

from dolfin import *

# Load mesh from file
mesh = UnitSquare(20,20)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define test functions
(v,q) = TestFunctions(W)

# Define trial functions
w     = Function(W)
(u,p) = (as_vector((w[0], w[1])), w[2])

# Set parameter values
nu = 0.01

# Define boundary conditions
noslip  = DirichletBC(W.sub(0), (0, 0), "x[0] < DOLFIN_EPS || x[0] > 1.0 - DOLFIN_EPS || x[1] < DOLFIN_EPS")
lid  = DirichletBC(W.sub(0), (1,0), "x[1] > 1.0 - DOLFIN_EPS")
pref = DirichletBC(W.sub(1), 0, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")

bc = [noslip, lid, pref]

# Tentative velocity step
F =   inner(grad(u)*u, v)*dx \
    + nu*inner(grad(u), grad(v))*dx \
    - div(v)*p*dx \
    + q*div(u)*dx

dw = TrialFunction(W)
dF = derivative(F, w, dw)

pde= VariationalProblem(F, dF, bc)
pde.solve(w)

(u,p) = w.split()
File("velocity.pvd") << u
