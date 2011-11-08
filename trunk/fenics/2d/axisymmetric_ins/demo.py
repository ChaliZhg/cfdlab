from dolfin import *

mu    = 0.01

omg   = 10   # angular speed in rad/sec

mesh = Mesh("annulus.xml")
sub_domains = MeshFunction("uint", mesh, "subdomains.xml")

n    = FacetNormal(mesh)
h    = CellSize(mesh)

Vry = VectorFunctionSpace(mesh, "CG", 2)
Vt  = FunctionSpace(mesh, "CG", 2)
Q   = FunctionSpace(mesh, "CG", 1)
X   = MixedFunctionSpace([Vry, Vt, Q])

v = Function(X)
w = TestFunction(X)

# Initial condition
v_init = Expression(("0", "0", "omg*x[0]", "0"), omg=omg)
v = interpolate(v_init, X)

ur  = v[0]
uy  = v[1]
ut  = v[2]
p   = v[3]

wr  = w[0]
wy  = w[1]
wt  = w[2]
q   = w[3]

W   = as_vector([w[0], w[1], w[2]])

r   = Expression("x[0]", element=Vt.ufl_element())

# Inviscid flux matrix
F = as_matrix([[ur**2, ur*uy], \
               [ur*uy, uy**2], \
               [ur*ut, uy*ut]])

# Divergence
rdivu   = ur   + r*ur.dx(0) + r*uy.dx(1)
rdivw   = wr   + r*wr.dx(0) + r*wy.dx(1)

# Stress tensor
rtau_rr = 2.0*mu*r*ur.dx(0)
rtau_ry = mu*r*(ur.dx(1) + uy.dx(0))
rtau_yr = rtau_ry
rtau_yy = 2.0*mu*r*uy.dx(1)
rtau_tt = 2.0*mu*ur
rtau_tr = mu*r*ut.dx(0) - mu*ut
rtau_rt = rtau_tr
rtau_ty = mu*r*ut.dx(1)
rtau_yt = rtau_ty

tau_tt = 2.0*mu*ur/r
tau_tr = mu*ut.dx(0) - mu*ut/r

# Viscous flux matrix, includes radius
G = as_matrix([[rtau_rr, rtau_ry], \
               [rtau_yr, rtau_yy], \
               [rtau_tr, rtau_ty]])

# Source term
S = as_vector([ut**2 + tau_tt, \
               0,                      \
               -ur*ut + tau_tr])

# Weak form
B =  Dx(r*F[i,j],j)*W[i]*dx  \
   + G[i,j]*Dx(W[i], j)*dx  \
   - S[i]*W[i]*dx           \
   - p*rdivw*dx             \
   - q*rdivu*dx

dv = TrialFunction(X)
dB = derivative(B, v, dv)

# Boundary condition
ut_bc_value = Expression("omg*x[0]", omg=omg)

ury_inner_bc  = DirichletBC(X.sub(0), (0,0), sub_domains, 0)
ury_outer_bc  = DirichletBC(X.sub(0), (0,0), sub_domains, 1)
ury_bottom_bc = DirichletBC(X.sub(0), (0,0), sub_domains, 2)
ury_top_bc    = DirichletBC(X.sub(0), (0,0), sub_domains, 3)

ut_inner_bc  = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 0)
ut_outer_bc  = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 1)
ut_bottom_bc = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 2)
ut_top_bc    = DirichletBC(X.sub(1), ut_bc_value, sub_domains, 3)

p_bc  = DirichletBC(X.sub(2), 0, "x[0]-1<DOLFIN_EPS && x[1]<DOLFIN_EPS", "pointwise")

bc    = [ury_inner_bc, ury_outer_bc, ury_bottom_bc, ury_top_bc, \
         ut_inner_bc, ut_outer_bc, ut_bottom_bc, ut_top_bc, \
         p_bc]

problem = NonlinearVariationalProblem(B, v, bc, dB)
solver  = NonlinearVariationalSolver(problem)

#solver.parameters["linear_solver"] = "gmres"
#itsolver = solver.parameters["newton_solver"]
#itsolver["absolute_tolerance"] = 1.0e-8

fury = File("ury.pvd", "compressed")
fut  = File("ut.pvd",  "compressed")
fp   = File("p.pvd",   "compressed")

solver.solve()
ury,ut,p = v.split()
fury << ury
fut  << ut
fp   << p
