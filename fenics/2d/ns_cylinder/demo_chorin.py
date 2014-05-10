"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *

# Load mesh from file
mesh = Mesh("cylinder_in_channel.xml")
sub_domains = MeshFunction("size_t", mesh, "subdomains.xml")

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
T = 100
nu = 0.01

# Define boundary conditions
noslip  = DirichletBC(V, (0, 0), sub_domains, 0)
inlet_u = DirichletBC(V, (1, 0), sub_domains, 1)
inlet_p = DirichletBC(Q, 0, sub_domains, 1)
bcu = [noslip, inlet_u]
bcp = [inlet_p]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

u0 = project(Constant((1,0)), V)

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

# vorticity
a4    = p*q*dx
A4    = assemble(a4)
vort  = Function(Q)
Lvort = (u1[0].dx(1) - u1[1].dx(0))*q*dx;

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
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", "ilu")
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "ilu")
    end()

    b4 = assemble(Lvort)
    solve(A4, vort.vector(), b4, "gmres", "ilu")
    
    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)
    plot(vort, title="Vorticity", rescale=True)

    # Save to file
    #ufile << u1

    # Move to next time step
    u0.assign(u1)
    p.update(t / T)
    t += dt

# Hold plot
interactive()
