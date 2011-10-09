"""
2-D oil reservoir problem using DGFEM and DFLU flux
Author: Praveen. C
DOES NOT WORL YET
"""

from dolfin import *

# Material properties
mu_o   = 1.0
mu_w   = 1.0
pinlet = 1.0
poutlet= 0.0
sinlet = 1.0
cinlet = 0.0

class Inlet(SubDomain):
   def inside(self, x, on_boundary):
      return ((x[0] < DOLFIN_EPS and x[1]-0.1 < DOLFIN_EPS) or \
              (x[1] < DOLFIN_EPS and x[0]-0.1 < DOLFIN_EPS)) and \
             on_boundary

class Outlet(SubDomain):
   def inside(self, x, on_boundary):
      return ((x[0]-1 > -DOLFIN_EPS and x[1]-0.9 > -DOLFIN_EPS) or \
              (x[1]-1 > -DOLFIN_EPS and x[0]-0.9 > -DOLFIN_EPS)) and \
             on_boundary

def mobility_water(s):
   return s**2/mu_w

def mobility_oil(s):
   return (1-s)**2/mu_o

def dflu1(sl, sr, cl, cr, vn):
   if vn > 0.0:
      m_w = mobility_water (sl)
      m_o = mobility_oil (sl)
      return vn * m_w / (m_w + m_o)
   else:
      m_w = mobility_water (sr)
      m_o = mobility_oil (sr)
      return vn * m_w / (m_w + m_o)

def dflu2(sl, sr, cl, cr, vn):
   flux = dflu1(sl, sr, cl, cr, vn)
   if vn > 0.0:
      return cl*flux
   else:
      return cr*flux

# Numerical parameters
np    = 50
theta = 0.5
dt    = 0.001

mesh = UnitSquare(np, np)
n    = FacetNormal(mesh)

sub_domains = MeshFunction("uint", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(100)
inlet = Inlet()
inlet.mark(sub_domains, 0)
outlet = Outlet()
outlet.mark(sub_domains, 1)

V = FunctionSpace(mesh, "DG", 0)
W = V * V

ts, tc = TestFunctions(W)

z = Function(W)
zold = Function(W)

# Set initial condition
ic = Expression(("0.0", "0.0"))
z  = project(ic, W)

# Saturation and concentration
s    = z[0]
c    = z[1]
sold = zold[0]
cold = zold[1]

stheta = (1-theta)*sold + theta*s
ctheta = (1-theta)*cold + theta*c

# Space for pressure
Q = FunctionSpace(mesh, "CG", 1)
r = TestFunction(Q)
q = TrialFunction(Q)
p = Function(Q)

lw_old = sold**2/mu_w
lo_old = (1-sold)**2/mu_o
lt_old = lw_old + lo_old

pa = lt_old*inner(grad(q), grad(r))*dx
pL = Constant(0)*r*dx
pbc_inlet  = DirichletBC(Q, pinlet,  sub_domains, 0)
pbc_outlet = DirichletBC(Q, poutlet, sub_domains, 1)
pbc = [pbc_inlet, pbc_outlet]

v = -lt_old*grad(p)
vn= avg(dot(v, n))

lw = stheta**2/mu_w
lo = (1-stheta)**2/mu_o
f  = lw/(lo + lw)
F  = v*f
cF = ctheta*F
H  = dflu1(stheta('+'), stheta('-'), ctheta('+'), ctheta('-'), vn)
cH = dflu2(stheta('+'), stheta('-'), ctheta('+'), ctheta('-'), vn)
Hi = dflu1(sinlet, sinlet, cinlet, cinlet, vn)
cHi= dflu2(sinlet, sinlet, cinlet, cinlet, vn)

dss = Measure("ds")[sub_domains]

L = (1/dt)*(s - sold)*ts*dx + (1/dt)*(s*c - sold*cold)*tc*dx \
    - inner(F, grad(ts))*dx - inner(cF, grad(tc))*dx \
    + H*jump(ts)*dS + cH*jump(tc)*dS \
    + inner(F, n)*ts*dss(1) + inner(cF, n)*tc*dss(1)

#   + inner(F, n)*ts*dss(0) + inner(cF, n)*tc*dss(0) \
#    + Hi*ts*dss(0) + cHi*tc*dss(0) \

sbc = DirichletBC(W.sub(0), sinlet, sub_domains, 0, "geometric")
cbc = DirichletBC(W.sub(1), cinlet, sub_domains, 0, "geometric")
bc  = [sbc, cbc]

dz  = TrialFunction(W)
dL  = derivative(L, z, dz)

problem = NonlinearVariationalProblem(L, z, bc, dL)
solver  = NonlinearVariationalSolver(problem)

t = 0
T = dt
while t < T:
   zold.assign(z)
   solve(pa == pL, p, pbc)
   File("p.pvd") << p
   #AA = assemble(L)
   #print AA.norm("l2")
   #quit()
   solver.solve()
   t = t + dt
