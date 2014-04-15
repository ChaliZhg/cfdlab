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
         value[0] = 0.4 * exp(-0.00001/ff)
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
           - hf*S*ds(3)
   Fdiv  =  - q*div(u)*dx
   return Fns + Ftemp + Fdiv

#-----------------------------------------------------------------------------
# contribution of heat flux term is not included here
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
      self.tdeg = udeg
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

      self.V = VectorFunctionSpace(mesh, "CG", self.udeg)  # velocity
      self.W = FunctionSpace(mesh, "CG", self.tdeg)        # temperature
      self.Q = FunctionSpace(mesh, "CG", self.pdeg)        # pressure
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
      #solver.parameters['newton_solver']['linear_solver'] = 'gmres'
      #solver.parameters['newton_solver']['absolute_tolerance'] = 1.0e-2
      #solver.parameters['newton_solver']['relative_tolerance'] = 1.0e-1
      #info(solver.parameters, True)

      for R in Relist:
         print "-----------------------------------------------------"
         print "Reynolds number = ", R 
         print "-----------------------------------------------------"
         Re.assign(R)
         solver.solve()
         u = as_vector((up[0],up[1]))
         d = assemble(div(u)*div(u)*dx)
         print "Divergence L2 norm = ", np.sqrt(d)

      # Save FE solution
      print "Saving FE solution into steady.xml"
      File("steady.xml") << up.vector()
      # Save vtk format
      u,T,p = up.split()
      print "Saving vtk files steady_u.pvd, steady_p.pvd, steady_T.pvd"
      File("steady_u.pvd") << u
      File("steady_p.pvd") << p
      File("steady_T.pvd") << T

   # Returns dof indices which are free
   # freeinds = free indices of velocity, temperature, pressure
   # pinds    = free indices of pressure
   def get_indices(self):
      # Collect all dirichlet boundary dof indices
      bcinds = []
      for b in self.bc:
         bcdict = b.get_boundary_values()
         bcinds.extend(bcdict.keys())

      # total number of dofs
      N = self.X.dim()

      # indices of free nodes
      freeinds = np.setdiff1d(range(N),bcinds).astype(np.int32)

      # pressure indices
      pinds = self.X.sub(2).dofmap().dofs()

      return freeinds, pinds

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
      print "Size of Ma =",Ma.shape

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      Aform = linear_form(Re,Gr,Pr,us,Ts,u,T,p,v,S,q)
      Aa = assemble(Aform)

      # Convert to sparse format
      rows, cols, values = Aa.data()
      Aa = sps.csc_matrix((values, cols, rows))
      print "Size of Aa =",Aa.shape

      freeinds,pinds = self.get_indices()

      print "Writing free indices into freeinds.txt"
      f = open('freeinds.txt','w')
      for item in freeinds:
          f.write("%d\n" % item)
      f.close()

      print "Writing pressure indices into pinds.txt"
      f = open('pinds.txt','w')
      for item in pinds:
          f.write("%d\n" % item)
      f.close()

      # indices of velocity control
      vinds = self.vbc.get_boundary_values().keys()

      # indices of temperature control
      tinds = self.tbc.get_boundary_values().keys()

      # mass matrix
      M = Ma[freeinds,:][:,freeinds]
      print "Size of M =",M.shape

      # stiffness matrix
      A = Aa[freeinds,:][:,freeinds]
      print "Size of A =",A.shape

      # velocity control operator
      ua = interpolate(allvar(1,1), self.X)
      ua = ua.vector().array()
      Bv = Aa[freeinds,:][:,vinds].dot(ua[vinds])
      print "Size of Bv =",Bv.shape[0]
      Bt = Aa[freeinds,:][:,tinds].dot(ua[tinds])
      print "Size of Bt =",Bt.shape[0]

      # heat flux control operator
      Bh = assemble(HeatFlux(1)*S*self.ds(3))
      Bh = Bh.array()
      Bh = Bh[freeinds]
      print "Size of Bh =",Bh.shape[0]

      B = np.column_stack((Bv,Bt,Bh))
      print "Size of B = ",B.shape

      # Save matrices in matlab format
      print "Saving linear system into linear.mat"
      sio.savemat('linear.mat', mdict={'M':M, 'A':A, 'B':B}, oned_as='column')

      # Compute eigenvalues/vectors
      vals, vecs = la.eigs(A, k=2, M=M, sigma=0, which='LR')
      for val in vals:
         print np.real(val), np.imag(val)
      
      # TODO: eigenvectors are complex
      ua = Function(self.X)

      # Save real part of eigenvector. << outputs only real part
      ua.vector()[freeinds] = vecs[:,0]
      File("evec1.xml") << ua.vector()
      u,T,p = ua.split()
      File("evec1_u.pvd") << u
      File("evec1_T.pvd") << T
      File("evec1_p.pvd") << p

      # Save imaginary part of eigenvector. << outputs only real part
      ua.vector()[freeinds] = vecs[:,0] * (-1j)
      File("evec2.xml") << ua.vector()
      u,T,p = ua.split()
      File("evec2_u.pvd") << u
      File("evec2_T.pvd") << T
      File("evec2_p.pvd") << p

   # Runs nonlinear model
   def run(self,with_control):
      up = Function(self.X)
      vp = TestFunction(self.X)

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hfamp,self.ds,up,vp)

      if with_control:
         # compute indices of velocity and temperature
         freeinds,pinds = self.get_indices()
         vTinds = np.setdiff1d(freeinds,pinds).astype(np.int32)
         gain = sio.loadmat('gain.mat')

      fhist = open('history.dat','w')

      ups = Function(self.X)
      File("steady.xml") >> ups.vector()
      us,Ts,ps = ups.split()
      KEs = assemble(0.5*inner(us,us)*dx)
      print 'Kinetic energy of steady state =', KEs

      # Set initial condition
      up1 = Function(self.X)
      up1.assign(ups)

      # Add perturbation using unstable eigenvector
      uppert = Function(self.X)
      File("evec1.xml") >> uppert.vector()
      up1.vector()[:] += 0.01 * uppert.vector().array()

      fu = File("u.pvd")
      ft = File("T.pvd")

      u,T,p = up1.split()
      fu << u
      ft << T
      KE  = assemble(0.5*inner(u,u)*dx)
      dKE = assemble(0.5*inner(u-us,u-us)*dx)
      print 'Kinetic energy =', KEs, KE, dKE
      fhist.write(str(0)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+"\n")

      dt = 0.01
      final_time = dt*2000
      time, iter = 0, 0

      if with_control:
        dy = up1.vector().array() - ups.vector().array()
        a = -np.dot(gain['Kt'], dy[vTinds])
        print a
        self.vamp.assign(a[0])
        self.tamp.assign(a[1])
        self.hfamp.assign(a[2])

      # First time step, we do backward euler
      B1 = (1/dt)*inner(up[0] - up1[0], vp[0])*dx     \
         + (1/dt)*inner(up[1] - up1[1], vp[1])*dx     \
         + (1/dt)*inner(up[2] - up1[2], vp[2])*dx + F

      dup = TrialFunction(self.X)
      dB1 = derivative(B1, up, dup)
      problem1 = NonlinearVariationalProblem(B1, up, self.bc, dB1)
      solver1  = NonlinearVariationalSolver(problem1)

      up.assign(up1)
      solver1.solve()
      iter += 1
      time += dt
      print 'Iter = {:5d}, t = {:f}'.format(iter, time)

      u,T,p = up.split()
      fu << u
      ft << T
      KE  = assemble(0.5*inner(u,u)*dx)
      dKE = assemble(0.5*inner(u-us,u-us)*dx)
      print 'Kinetic energy =', KEs, KE, dKE
      fhist.write(str(time)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+"\n")
      print '--------------------------------------------------------------'

      # From now on use BDF2
      up2 = Function(self.X)
      B2 = (1/dt)*inner(1.5*up[0] - 2.0*up1[0] + 0.5*up2[0], vp[0])*dx     \
         + (1/dt)*inner(1.5*up[1] - 2.0*up1[1] + 0.5*up2[1], vp[1])*dx     \
         + (1/dt)*inner(1.5*up[2] - 2.0*up1[2] + 0.5*up2[2], vp[2])*dx + F
      dB2 = derivative(B2, up, dup)
      problem2 = NonlinearVariationalProblem(B2, up, self.bc, dB2)
      solver2  = NonlinearVariationalSolver(problem2)

      while time < final_time:
         up2.assign(up1)
         up1.assign(up)
         if with_control:
            dy = up1.vector().array() - ups.vector().array()
            a = -np.dot(gain['Kt'], dy[vTinds])
            print a
            self.vamp.assign(a[0])
            self.tamp.assign(a[1])
            self.hfamp.assign(a[2])
         solver2.solve()
         iter += 1
         time += dt
         print 'Iter = {:5d}, t = {:f}'.format(iter, time)
         u,T,p = up.split()
         if iter%10 == 0:
            fu << u
            ft << T
         KE  = assemble(0.5*inner(u,u)*dx)
         dKE = assemble(0.5*inner(u-us,u-us)*dx)
         print 'Kinetic energy =', KEs, KE, dKE
         fhist.write(str(time)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+"\n")
         fhist.flush()
         print '--------------------------------------------------------------'

      fhist.close()
