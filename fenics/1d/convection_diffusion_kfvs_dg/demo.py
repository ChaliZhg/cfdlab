from dolfin import *

nc   = 20
xmin = -1.0
xmax =  1.0

c      = 1.0
mu     = 0.01
degree = 1
dt     = 0.001
T      = 100*dt

mesh = Interval(nc, xmin, xmax)

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return bool(x[0]-xmin <  DOLFIN_EPS and \
                    x[0]-xmin > -DOLFIN_EPS and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - xmax + xmin

Vh = FunctionSpace(mesh, "CG", degree)

u = Function(Vh)
uo= Function(Vh)
v = TrialFunction(Vh)
w = TestFunction(Vh)

uinit = Expression("-sin(pi*x[0])")
u = project(uinit, Vh)


# Weak formulation
a =   v*w*dx

R =   c*u*Dx(w,0)*dx        \
    - c*u*w*ds              \
    - mu*Dx(u,0)*Dx(w,0)*dx \
    + mu*Dx(u,0)*w*ds
R = u*w*dx + dt*R

pbc= PeriodicBoundary()
bc = PeriodicBC(Vh, pbc)

#problem = LinearVariationalProblem(B==R, u, bc)
#solver  = LinearVariationalSolver(problem)

ff = File("u.pvd", "compressed")

t = 0
it= 0
while t < T:
   uo.assign(u)
   solve(a == R, u, bc)
   t = t + dt
   it= it + 1
   print "Iter=", it, ", t=", t
   ff << u
