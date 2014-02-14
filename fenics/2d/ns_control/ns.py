from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.io as sio

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
      value[0] = self.amp
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
         value[0] = -exp(-0.0001/ff)
      else:
         value[0] = 0.0
      value[0] = self.amp * value[0]/exp(-0.0001/0.01**2)
      value[1] = 0.0
   def value_shape(self):
      return (2,)
#-----------------------------------------------------------------------------
def nonlinear_form(Re,Gr,Pr,hfamp,ds,u,T,p,v,S,q):
   nu    = 1/Re
   k     = 1/(Re*Pr)
   G     = Gr/Re**2
   hs    = HeatSource()
   hf    = HeatFlux(hfamp)
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
      self.pdeg = udeg - 1
      self.Re = Re
      self.Gr = Gr
      self.Pr = Pr
      self.ds = Measure("ds")[sub_domains]

      # velocity control amplitude
      self.vamp  = Constant(0)
      # heat flux control amplitude
      self.hfamp = Constant(0)

      self.V = VectorFunctionSpace(mesh, "CG", self.udeg)
      self.W = FunctionSpace(mesh, "CG", self.udeg)
      self.Q = FunctionSpace(mesh, "CG", self.pdeg)
      self.X = MixedFunctionSpace([self.V, self.W, self.Q])

      # Velocity bc
      noslipbc1 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 0)
      noslipbc2 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 3)
      # velocity control boundary
      gs        = velocity(self.vamp)
      conbc     = DirichletBC(self.X.sub(0), gs,    sub_domains, 2)

      # Temperature bc
      tbc1    = DirichletBC(self.X.sub(1), 0.0, sub_domains, 0)
      tbc2    = DirichletBC(self.X.sub(1), 0.0, sub_domains, 2)

      self.bc = [noslipbc1, noslipbc2, conbc, tbc1, tbc2]
      self.vbc= conbc

   def steady_state(self, Relist):
      up = Function(self.X)
      u  = as_vector((up[0],up[1]))   # velocity
      T  = up[2]                      # temperature
      p  = up[3]                      # pressure
      # Define test functions
      (v,S,q) = TestFunctions(self.X)
      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hfamp,self.ds,u,T,p,v,S,q)

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
      File("steady.xml") << up.vector()
      # Save vtk format
      u,T,p = up.split()
      File("steady_u.pvd") << u
      File("steady_p.pvd") << p
      File("steady_T.pvd") << T

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
      print Ma

      rows, cols, values = Ma.data()
      Ma = sps.csc_matrix((values, cols, rows))
      N = Ma.shape[0]
      print Ma.shape[0], Ma.shape[1]

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      Aform = linear_form(Re,Gr,Pr,us,Ts,u,T,p,v,S,q)
      Aa = assemble(Aform)
      print Aa

      # Convert to sparse format
      rows, cols, values = Aa.data()
      Aa = sps.csc_matrix((values, cols, rows))
      N = Aa.shape[0]
      print "Size of Aa =",Aa.shape[0],Aa.shape[1]

      # Collect all dirichlet boundary dof indices
      bcinds = []
      for b in self.bc:
         bcdict = b.get_boundary_values()
         bcinds.extend(bcdict.keys())

      # indices of free nodes
      innerinds = np.setdiff1d(range(N),bcinds).astype(np.int32)

      # indices of velocity control
      vinds = []
      bcdict = self.vbc.get_boundary_values()
      vinds.extend(bcdict.keys())

      Mc = Ma[innerinds,:][:,innerinds]
      Ac = Aa[innerinds,:][:,innerinds]
      print "Size of Mc =",Mc.shape[0],Mc.shape[1]
      print "Size of Ac =",Ac.shape[0],Ac.shape[1]

      Bc = Aa[innerinds,:][:,vinds]

      # Save matrices in matlab format
      sio.savemat('Mc.mat', mdict={'Mc': Mc})
      sio.savemat('Ac.mat', mdict={'Ac': Ac})
