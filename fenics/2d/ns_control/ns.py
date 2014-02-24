from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.io as sio
import scipy.sparse.linalg as la

#-----------------------------------------------------------------------------
class HeatSource(Expression):
   def eval(self, value, x):
      value[0] = 7.0*sin(2.0*pi*x[0])*cos(2.0*pi*x[1])
      return value
#-----------------------------------------------------------------------------
class HeatFlux(Expression):
   def __init__(self, amp):
      self.amp = amp
   def eval(self, value, x):
      ff = ((x[0]-0.4)*(x[0]-0.6))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[0]>0.4 and x[0]<0.6:
         value[0] = exp(-0.00001/ff)
      else:
         value[0] = 0.0
      value[0] *= self.amp
      return value
#-----------------------------------------------------------------------------
# Velocity at inflow boundary
class velocity(Expression):
   def __init__(self, amp):
      self.amp = amp
   def eval(self, value, x):
      ff = ((x[1]-0.7)*(x[1]-0.9))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[1]>0.7 and x[1]<0.9:
         value[0] = exp(-0.0001/ff)
      else:
         value[0] = 0.0
      value[0] *= self.amp
      return value
#-----------------------------------------------------------------------------
# temperature at inflow boundary
class temperature(Expression):
   def __init__(self, amp):
      self.amp = amp
   def eval(self, value, x):
      ff = ((x[1]-0.7)*(x[1]-0.9))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[1]>0.7 and x[1]<0.9:
         value[0] = 0.2 * exp(-0.00001/ff)
      else:
         value[0] = 0.0
      value[0] *= self.amp
      return value
#-----------------------------------------------------------------------------
# all variables at inflow boundary
# we only need x velocity and temperature at inflow boundary
# which is the control boundary
class allvar(Expression):
   def __init__(self, vamp, tamp):
      self.vamp = vamp
      self.tamp = tamp
   def eval(self, value, x):
      ff = ((x[1]-0.7)*(x[1]-0.9))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
         value[2] = 0.0
      elif x[1]>0.7 and x[1]<0.9:
         value[0] = exp(-0.0001/ff)
         value[2] = 0.2 * exp(-0.00001/ff)
      else:
         value[0] = 0.0
         value[2] = 0.0
      value[0] *= self.vamp  # x velocity
      value[1] = 0.0         # y velocity
      value[2] *= self.tamp  # temperature
      value[3] = 0.0         # pressure
      return value
   def value_shape(self):
      return (4,)
#-----------------------------------------------------------------------------
def nonlinear_form(Re,Gr,Pr,hfamp,ds,up,vp):
   nu    = 1/Re
   k     = 1/(Re*Pr)
   G     = Gr/Re**2
   hs    = HeatSource()
   hf    = HeatFlux(hfamp)
   # Trial function
   u  = as_vector((up[0],up[1]))   # velocity
   T  = up[2]                      # temperature
   p  = up[3]                      # pressure
   # Test function
   v  = as_vector((vp[0],vp[1]))   # velocity
   S  = vp[2]                      # temperature
   q  = vp[3]                      # pressure
   Fns   =   inner(grad(u)*u, v)*dx      \
          + nu*inner(grad(u),grad(v))*dx \
          - div(v)*p*dx                  \
          - G*T*v[1]*dx
   Ftemp =   inner(grad(T),u)*S*dx       \
           + k*inner(grad(T),grad(S))*dx \
           - hs*S*dx                     \
           + hf*S*ds(3)
   Fdiv  =  - q*div(u)*dx
   return Fns + Ftemp + Fdiv

#-----------------------------------------------------------------------------
def linear_form(Re,Gr,Pr,us,Ts,u,T,p,v,S,q):
   nu    = 1/Re
   k     = 1/(Re*Pr)
   G     = Gr/Re**2
   Aform = - inner( grad(u)*us, v )*dx \
           - inner( grad(us)*u, v )*dx \
           + div(v)*p*dx \
           - nu*inner( grad(u), grad(v) )*dx \
           + G*T*v[1]*dx \
           + div(u)*q*dx \
           - inner(grad(Ts),u)*S*dx \
           - inner(grad(T),us)*S*dx \
           - k*inner(grad(T),grad(S))*dx
   return Aform
#-----------------------------------------------------------------------------
class NSProblem():
   def __init__(self, udeg, Re, Gr, Pr):
      mesh = Mesh("square.xml")
      sub_domains = MeshFunction("size_t", mesh, "subdomains.xml")

      self.udeg = udeg
      self.tdeg = udeg - 1
      self.pdeg = udeg - 1
      self.Re = Re
      self.Gr = Gr
      self.Pr = Pr
      self.ds = Measure("ds")[sub_domains]

      # velocity control amplitude
      self.vamp  = Constant(0)
      # temperature control amplitude
      self.tamp  = Constant(0)
      # heat flux control amplitude
      self.hfamp = Constant(0)

      self.V = VectorFunctionSpace(mesh, "CG", self.udeg)
      self.W = FunctionSpace(mesh, "CG", self.tdeg)
      self.Q = FunctionSpace(mesh, "CG", self.pdeg)
      self.X = MixedFunctionSpace([self.V, self.W, self.Q])

      # Velocity bc
      noslipbc1 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 0)
      noslipbc2 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 3)
      # velocity control boundary
      gs        = velocity(self.vamp)
      vconbc    = DirichletBC(self.X.sub(0).sub(0), gs, sub_domains, 2)
      yvbc      = DirichletBC(self.X.sub(0).sub(1), 0,  sub_domains, 2)

      # Temperature bc
      tbc1    = DirichletBC(self.X.sub(1), 0.0, sub_domains, 0)
      ts      = temperature(self.tamp)
      tbc2    = DirichletBC(self.X.sub(1), ts, sub_domains, 2)

      self.bc = [noslipbc1, noslipbc2, vconbc, yvbc, tbc1, tbc2]
      self.vbc= vconbc
      self.tbc= tbc2

   # solves steady state equations
   def steady_state(self, Relist):
      up = Function(self.X)
      vp = TestFunction(self.X)
      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hfamp,self.ds,up,vp)

      dup = TrialFunction(self.X)
      dF  = derivative(F, up, dup)
      problem = NonlinearVariationalProblem(F, up, self.bc, dF)
      solver  = NonlinearVariationalSolver(problem)

      for R in Relist:
         print "-----------------------------------------------------"
         print "Reynolds number = ", R 
         print "-----------------------------------------------------"
         Re.assign(R)
         solver.solve()

      # Save FE solution
      print "Saving FE solution into steady.xml"
      File("steady.xml") << up.vector()
      # Save vtk format
      u,T,p = up.split()
      print "Saving vtk files steady_u.pvd, steady_p.pvd, steady_T.pvd"
      File("steady_u.pvd") << u
      File("steady_p.pvd") << p
      File("steady_T.pvd") << T

   # Generate linear state representation
   def linear_system(self):
      parameters.linear_algebra_backend = "uBLAS"
       
      # Load Stationary solution from file
      ups = Function(self.X)
      File("steady.xml") >> ups.vector()
      us= as_vector((ups[0],ups[1]));
      Ts= ups[2]

      u,T,p = TrialFunctions(self.X)
      v,S,q = TestFunctions(self.X)

      # Mass matrix
      Ma = assemble(inner(u,v)*dx + T*S*dx)

      rows, cols, values = Ma.data()
      Ma = sps.csc_matrix((values, cols, rows))
      print "Size of Ma =",Ma.shape[0], Ma.shape[1]

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      Aform = linear_form(Re,Gr,Pr,us,Ts,u,T,p,v,S,q)
      Aa = assemble(Aform)

      # Convert to sparse format
      rows, cols, values = Aa.data()
      Aa = sps.csc_matrix((values, cols, rows))
      print "Size of Aa =",Aa.shape[0],Aa.shape[1]

      # Collect all dirichlet boundary dof indices
      bcinds = []
      for b in self.bc:
         bcdict = b.get_boundary_values()
         bcinds.extend(bcdict.keys())

      # indices of free nodes
      N = Aa.shape[0]
      innerinds = np.setdiff1d(range(N),bcinds).astype(np.int32)

      # indices of velocity control
      vinds = []
      bcdict = self.vbc.get_boundary_values()
      vinds.extend(bcdict.keys())

      # indices of temperature control
      tinds = []
      bcdict = self.tbc.get_boundary_values()
      tinds.extend(bcdict.keys())

      # mass matrix
      M = Ma[innerinds,:][:,innerinds]
      print "Size of M =",M.shape[0],M.shape[1]

      # stiffness matrix
      A = Aa[innerinds,:][:,innerinds]
      print "Size of A =",A.shape[0],A.shape[1]

      # velocity control operator
      ua = interpolate(allvar(1,1), self.X)
      ua = ua.vector().array()
      Bv = Aa[innerinds,:][:,vinds].dot(ua[vinds])
      print "Size of Bv =",Bv.shape[0]
      Bt = Aa[innerinds,:][:,tinds].dot(ua[tinds])
      print "Size of Bt =",Bt.shape[0]

      # heat flux control operator
      Bh = assemble(HeatFlux(1)*S*self.ds(3))
      Bh = Bh.array()
      Bh = Bh[innerinds]
      print "Size of Bh =",Bh.shape[0]

      # Save matrices in matlab format
      print "Saving linear system into linear.mat"
      sio.savemat('linear.mat', mdict={'M':M, 'A':A, 'Bv':Bv, 'Bt':Bt, 'Bh':Bh})

      # Compute eigenvalues/vectors
      vals, vecs = la.eigs(A, k=2, M=M, sigma=0, which='LR')
      for val in vals:
         print np.real(val), np.imag(val)
      
      # TODO: eigenvectors are complex
      ua = Function(self.X)

      ua.vector()[innerinds] = vecs[:,0]
      File("evec1.xml") << ua.vector()
      u,T,p = ua.split()
      File("evec1_u.pvd") << u
      File("evec1_T.pvd") << T

      ua.vector()[innerinds] = vecs[:,1]
      File("evec2.xml") << ua.vector()
      u,T,p = ua.split()
      File("evec2_u.pvd") << u
      File("evec2_T.pvd") << T

   # Runs nonlinear model
   def run(self):
      up = Function(self.X)
      vp = TestFunction(self.X)

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hfamp,self.ds,up,vp)

      # Set initial condition
      upold = Function(self.X)
      File("steady.xml") >> upold.vector()
      uppert = Function(self.X)
      File("evec1.xml") >> uppert.vector()
      upold.vector()[:] += 0.01 * uppert.vector().array()

      fu = File("u.pvd")
      ft = File("T.pvd")

      u,T,p = upold.split()
      fu << u
      ft << T

      dt = 0.01
      final_time = dt*100

      # First time step, we do backward euler
      B1 = (1/dt)*inner(up[0] - upold[0], vp[0])*dx     \
         + (1/dt)*inner(up[1] - upold[1], vp[1])*dx     \
         + (1/dt)*inner(up[2] - upold[2], vp[2])*dx + F

      dup = TrialFunction(self.X)
      dB1 = derivative(B1, up, dup)
      problem1 = NonlinearVariationalProblem(B1, up, self.bc, dB1)
      solver1  = NonlinearVariationalSolver(problem1)

      time, iter = 0, 0
      up.assign(upold)
      solver1.solve()
      iter += 1
      time += dt

      u,T,p = up.split()
      fu << u
      ft << T

      while time < final_time:
         up.assign(upold)
         solver1.solve()
         iter += 1
         time += dt
         u,T,p = up.split()
         fu << u
         ft << T
