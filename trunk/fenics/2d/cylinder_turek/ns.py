"""
Flow over cylinder in channel
Test case from Turek

Convective term approximated by extrapolation
BDF1 in first step, BDF2 subsequently
"""
from dolfin import *

class inlet_velocity(Expression):
   def __init__(self, t=0.0):
      self.t = t
   def eval(self, value, x):
      yc = x[1]/0.41
      value[0] = 6.0*yc*(1.0-yc)
      value[1] = 0.0
      return value
   def value_shape(self):
      return (2,)

class initial_condition(Expression):
   def eval(self, value, x):
      yc = x[1]/0.41
      value[0] = 6.0*yc*(1.0-yc)
      value[1] = 0.0
      value[2] = 0.0
      return value
   def value_shape(self):
      return (3,)

mesh = Mesh("cylinder.xml")
boundaries = MeshFunction("size_t", mesh, "cylinder_facet_region.xml")

# Function space
udeg = 2
pdeg = udeg - 1

V = VectorFunctionSpace(mesh, 'CG', udeg)
W = FunctionSpace(mesh, 'CG', pdeg)
X = V * W

print "Velocity dofs = ", V.dim()
print "Pressure dofs = ", W.dim()
print "Total    dofs = ", X.dim()

# Solution variables
up0 = Function(X)
up1 = Function(X)
up2 = Function(X)

# Trial functions
up  = TrialFunction(X)
u   = as_vector((up[0],up[1]))
p   = up[2]

# Test functions
vq  = TestFunction(X)
v   = as_vector((vq[0],vq[1]))
q   = vq[2]

# Boundary condition
vin  = inlet_velocity(0.0)
bcin = DirichletBC(X.sub(0), vin,   boundaries, 1)
bccyl= DirichletBC(X.sub(0), (0,0), boundaries, 2)
bcwal= DirichletBC(X.sub(0), (0,0), boundaries, 4)
bcs  = [bcin, bccyl, bcwal]

Re = 50.0
nu = Constant(1.0/Re)
dt = 0.01
idt= Constant(1.0/dt)

up0.interpolate(initial_condition())
u0 = as_vector((up0[0], up0[1]))
u1 = as_vector((up1[0], up1[1]))

t  = 0.0
Tf = 10.0
it = 0
fu = File("u.pvd")

# First time step: BDF1
# Predicted velocity
us = u0

F1 = idt*inner(u - u0, v)*dx       \
   + inner(grad(us)*us, v)*dx      \
   - p*div(v)*dx                   \
   + nu*inner(grad(u), grad(v))*dx \
   - q*div(u)*dx

a  = lhs(F1)
L  = rhs(F1)

A  = assemble(a)
solver = LUSolver(A)
solver.parameters['reuse_factorization'] = True

b  = assemble(L)
[bc.apply(A,b) for bc in bcs]
solver.solve(up1.vector(), b)
up0.assign(up1)
t += dt
it+= 1

# Now switch to BDF2
# Predicted velocity
us = 2.0*u1 - u0

F2 = idt*inner(1.5*u - 2.0*u1 + 0.5*u0, v)*dx       \
   + inner(grad(us)*us, v)*dx      \
   - p*div(v)*dx                   \
   + nu*inner(grad(u), grad(v))*dx \
   - q*div(u)*dx

a  = lhs(F2)
L  = rhs(F2)

A  = assemble(a)
solver = LUSolver(A)
solver.parameters['reuse_factorization'] = True

while t < Tf:
    b  = assemble(L)
    [bc.apply(A,b) for bc in bcs]
    solver.solve(up2.vector(), b)
    up0.assign(up1)
    up1.assign(up2)
    t += dt
    it+= 1
    print "it = %6d,     t = %12.6e" % (it,t)
    if it%10 == 0:
        u,p = up2.split()
        fu << u
