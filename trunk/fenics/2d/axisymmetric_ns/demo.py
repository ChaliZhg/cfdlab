from dolfin import *

R     = 287.0
gamma = 1.4
mu    = 0.01
Pr    = 0.72

Cp    = gamma*R/(gamma-1)
k     = mu*Cp/Pr
gamma1= gamma/(gamma-1)

omg   = 100 # angular speed in rad/sec
r1    = 1   # Radius of inner cylinder
r2    = 2   # Radius of outer cylinder
ht    = 2   # Height of cylinder
Tbc   = 300 # Boundary temperature

npr   = 20
npy   = 40
dt    = 0.0001
Tf    = 1000*dt

def Boundary(x, on_boundary):
   return on_boundary

#mesh = Rectangle(r1, 0, r2, ht, npr, npy)
mesh = Mesh("annulus.xml")

n    = FacetNormal(mesh)
h    = CellSize(mesh)
hmin = mesh.hmin()

V = FunctionSpace(mesh, "CG", 1)
Vh= MixedFunctionSpace([V, V, V, V, V])

v = Function(Vh)
vo= Function(Vh)
w = TestFunction(Vh)

# Initial condition
v_init = Expression(("1", "0", "0", "omg*x[0]", "Tbc"), omg=omg, Tbc=Tbc)
v = interpolate(v_init, Vh)
vo.assign(v)

rho = v[0]
ur  = v[1]
uy  = v[2]
ut  = v[3]
T   = v[4]

p   = rho*R*T
e   = p/(gamma-1) + 0.5*rho*(ur**2 + uy**2 + ut**2)
H   = gamma1*R*T + 0.5*(ur**2 + uy**2 + ut**2)

U   = as_vector([rho, rho*ur, rho*uy, rho*ut, e])
Uo  = replace(U, {v:vo})

r   = Expression("x[0]", element=V.ufl_element())

# Inviscid flux matrix
F = as_matrix([[rho*ur,      rho*uy     ], \
               [p+rho*ur**2, rho*ur*uy  ], \
               [rho*ur*uy,   p+rho*uy**2], \
               [rho*ur*ut,   rho*uy*ut  ], \
               [(e+p)*ur,    (e+p)*uy   ]])

# Inviscid boundary flux
Fb = as_vector([0, p*n[0], p*n[1], 0, 0])

# Divergence
rdivu   = ur   + r*ur.dx(0) + r*uy.dx(1)
divu    = ur/r + ur.dx(0)   + uy.dx(1)

# Stress tensor
rtau_rr = 2.0*mu*r*ur.dx(0) - (2.0/3.0)*mu*rdivu
rtau_ry = mu*r*(ur.dx(1) + uy.dx(0))
rtau_yr = rtau_ry
rtau_yy = 2.0*mu*r*uy.dx(1) - (2.0/3.0)*mu*rdivu
rtau_tt = 2.0*mu*ur - (2.0/3.0)*mu*rdivu
rtau_tr = mu*r*ut.dx(0) - mu*ut
rtau_rt = rtau_tr
rtau_ty = mu*r*ut.dx(1)
rtau_yt = rtau_ty

tau_tt = 2.0*mu*ur/r - (2.0/3.0)*mu*divu
tau_tr = mu*ut.dx(0) - mu*ut/r

# Heat flux vector
rq_r = -k*r*T.dx(0)
rq_y = -k*r*T.dx(1)

# Energy flux
er = rtau_rr*ur + rtau_ry*uy + rtau_rt*ut - rq_r
ey = rtau_yr*ur + rtau_yy*uy + rtau_yt*ut - rq_y

# Viscous flux matrix, includes radius
G = as_matrix([[0,       0,     ], \
               [rtau_rr, rtau_ry], \
               [rtau_yr, rtau_yy], \
               [rtau_tr, rtau_ty], \
               [er,      ey     ]])

# Source term
S = as_vector([0,                      \
               p + rho*ut**2 + tau_tt, \
               0,                      \
               -rho*ur*ut + tau_tr,    \
               0])

Ao = as_matrix([ \
     [1,     0,      0,      0,      0               ], \
     [ur,    rho,    0,      0,      0               ], \
     [uy,    0,      rho,    0,      0               ], \
     [ut,    0,      0,      rho,    0               ], \
     [e/rho, rho*ur, rho*uy, rho*ut, rho*R/(gamma-1) ]  \
     ])
Ao = replace(Ao, {v:vo})

# Weak form
B_GAL = (1/dt)*r*inner(U-Uo,w)*dx \
        - r*F[i,j]*Dx(w[i], j)*dx      \
        + G[i,j]*Dx(w[i], j)*dx        \
        - S[i]*w[i]*dx
B_BCS = r*Fb[i]*w[i]*ds - w[i]*G[i,j]*n[j]*ds

Ar = as_matrix([ \
      [ur,        rho,           0,         0,         0              ], \
      [R*T+ur**2, 2*rho*ur,      0,         0,         rho*R          ], \
      [ur*uy,     rho*uy,        rho*ur,    0,         0              ], \
      [ur*ut,     rho*ut,        0,         rho*ur,    0              ], \
      [H*ur,      rho*(H+ur**2), rho*ur*uy, rho*ur*ut, gamma1*R*rho*ur]  \
      ])

Ay = as_matrix([ \
      [uy,        0,         rho,           0,         0              ], \
      [ur*uy,     rho*uy,    rho*ur,        0,         0              ], \
      [R*T+uy**2, 0,         2*rho*uy,      0,         rho*R          ], \
      [uy*ut,     0,         rho*ut,        rho*uy,    0              ], \
      [H*uy,      rho*uy*ur, rho*(H+uy**2), rho*uy*ut, gamma1*R*rho*uy] \
      ])

RES   = as_vector(r*Ao[i,j]*(v[j]-vo[j])/dt + Dx(r*F[i,j],j) - Dx(G[i,j],j) - S[i], i)

# Ar and Ay must be transposed
PSUP  = as_vector(Ar[j,i]*Dx(w[j],0) + Ay[j,i]*Dx(w[j],1), i)
delta = h/(sqrt(ur**2 + uy**2) + sqrt(gamma*R*T))
PSUP  = delta*PSUP
PSUP  = replace(PSUP, {v:vo})

# For derivative, we consider SUPG terms as constant
B_SUP = PSUP[i]*RES[i]*dx
B     = B_GAL + B_BCS + B_SUP

dv = TrialFunction(Vh)
dB = derivative(B, v, dv)

# Boundary condition
ur_bc_value = Expression("0")
uy_bc_value = Expression("0")
ut_bc_value = Expression("omg*x[0]", omg=omg)
T_bc_value  = Expression("Tbc", Tbc=Tbc)
ur_bc = DirichletBC(Vh.sub(1), ur_bc_value, Boundary)
uy_bc = DirichletBC(Vh.sub(2), uy_bc_value, Boundary)
ut_bc = DirichletBC(Vh.sub(3), ut_bc_value, Boundary)
T_bc  = DirichletBC(Vh.sub(4), T_bc_value,  Boundary)
bc    = [ur_bc, uy_bc, ut_bc, T_bc]

problem = NonlinearVariationalProblem(B, v, bc, dB)
solver  = NonlinearVariationalSolver(problem)

solver.parameters["linear_solver"] = "gmres"
itsolver = solver.parameters["newton_solver"]
itsolver["absolute_tolerance"] = 1.0e-8

frho = File("rho.pvd", "compressed")
fur  = File("ur.pvd",  "compressed")
fuy  = File("uy.pvd",  "compressed")
fut  = File("ut.pvd",  "compressed")
fT   = File("T.pvd",   "compressed")

iter = 0
t    = 0
while t < Tf:
   vo.assign(v)
   solver.solve()
   t    = t + dt
   iter = iter + 1
   print "Iter=", iter, ", t=", t
   if iter % 5 == 0:
      rho,ur,uy,ut,T = v.split()
      rho.rename("rho", "Density")
      frho << rho
      fur  << ur
      fuy  << uy
      fut  << ut
      fT   << T
