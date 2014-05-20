"""
This demo program solves the incompressible Navier-Stokes equations
for lid-driven cavity problem using Chorin's splitting method.
"""

from dolfin import *

# Velocity at inflow boundary
class inlet_velocity(Expression):
   def __init__(self, t=0.0):
      self.t = t
   def eval(self, value, x):
      yc = x[1]/0.4
      value[0] = 6.0*yc*(1.0-yc)
      value[1] = 0.0
      return value
   def value_shape(self):
      return (2,)

# Load mesh from file
mesh = Mesh("cyl.xml")
boundaries = MeshFunction("size_t", mesh, "cyl_facet_region.xml")

File("boundaries.pvd") << boundaries

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
T  = 100.0
nu = 1.0/150.0

# Define boundary conditions

# Velocity bc
vin      = inlet_velocity(0)
inletbc  = DirichletBC(V, vin, boundaries, 1)
noslipbc = DirichletBC(V, (0,0), boundaries, 3)
# velocity control boundary
vconbc1  = DirichletBC(V, (0,0), boundaries, 4)
vconbc2  = DirichletBC(V, (0,0), boundaries, 5)
bcu      = [noslipbc, inletbc, vconbc1, vconbc2]

inlet_p = DirichletBC(Q, 0, boundaries, 2)
bcp = [inlet_p]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

u0 = project(vin, V)

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
ufile    = File("velocity.pvd")
pfile    = File("pressure.pvd")
vortfile = File("vorticity.pvd")

# Time-stepping
it= 0
t = 0.0
p = Progress("Time-stepping")
while t < T + DOLFIN_EPS:
    print "----------------------------------------------------------"
    vin.t = t + dt
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
    #plot(vort, title="Vorticity", rescale=True)

    # Move to next time step
    u0.assign(u1)
    p.update(t / T)
    t += dt
    it += 1

    # Save to file
    if it%10 == 0:
      ufile << u1
      pfile << p1
      vortfile << vort
