"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *

# Load mesh from file
mesh = UnitSquare(20,20)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.01
T = 3
nu = 0.01

# Define boundary conditions
noslip  = DirichletBC(V, (0, 0), "x[0] < DOLFIN_EPS || x[0] > 1.0 - DOLFIN_EPS || x[1] < DOLFIN_EPS")
lid  = DirichletBC(V, (1,0), "x[1] > 1.0 - DOLFIN_EPS")
bcu = [noslip, lid]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Create files for storing solution
ufile = File("velocity.pvd")

# Time-stepping
t = dt
p = Progress("Time-stepping")
while t < T + DOLFIN_EPS:

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "ilu")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    solve(A2, p1.vector(), b2, "gmres", "amg_hypre")
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "ilu")
    end()

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1

    # Move to next time step
    u0.assign(u1)
    p.update(t / T)
    t += dt

# Hold plot
interactive()
