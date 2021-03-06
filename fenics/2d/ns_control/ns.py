from dolfin import *
from bmark import *
import numpy as np
import scipy.sparse as sps
import scipy.io as sio
import scipy.sparse.linalg as la

# Position of some boundary zones
x1,x2,y1,y2 = 0.4,0.6,0.1,0.4

#-----------------------------------------------------------------------------
def g(s):
    if s < -1.0:
        return 0.0
    elif s > 1.0:
        return 1.0
    else:
        return 0.5 + s*(0.9375 - s*s*(0.625 - 0.1875*s*s))
#-----------------------------------------------------------------------------
# Smooth ramp from 0 to 1; middle of ramp is at T, width dt
# Put large negative value to disable ramp
def f(t):
    T, dt = -20.0, 2.0
    t1 = (t - T)/dt
    return g(t1)
#-----------------------------------------------------------------------------
class UnitNormal(Expression):
   def eval(self, value, x):
       if near(x[0],0.0):
           value[0] = -1.0
           value[1] =  0.0
       elif near(x[0],1.0):
           value[0] =  1.0
           value[1] =  0.0
       elif near(x[1],0.0):
           value[0] =  0.0
           value[1] = -1.0
       elif near(x[1],1.0):
           value[0] =  0.0
           value[1] =  1.0
       else:
           value[0] =  0.0
           value[1] =  0.0
       return value
   def value_shape(self):
      return (2,)
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
      ff = ((x[0]-x1)*(x[0]-x2))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[0]>x1 and x[0]<x2:
         value[0] = 0.4 * exp(-0.00001/ff)
      else:
         value[0] = 0.0
      value[0] *= self.amp
      return value
#-----------------------------------------------------------------------------
# Velocity at inflow boundary
class velocity(Expression):
   def __init__(self, amp, y3, y4):
      self.amp = amp
      self.y3  = y3
      self.y4  = y4
   def eval(self, value, x):
      ff = ((x[1]-self.y3)*(x[1]-self.y4))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[1]>self.y3 and x[1]<self.y4:
         value[0] = exp(-0.0001/ff)
      else:
         value[0] = 0.0
      value[0] *= self.amp
      return value
#-----------------------------------------------------------------------------
# temperature at inflow boundary
class temperature(Expression):
   def __init__(self, amp, y3, y4):
      self.amp = amp
      self.y3  = y3
      self.y4  = y4
   def eval(self, value, x):
      ff = ((x[1]-self.y3)*(x[1]-self.y4))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
      elif x[1]>self.y3 and x[1]<self.y4:
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
   def __init__(self, vamp, tamp, y3, y4):
      self.vamp = vamp
      self.tamp = tamp
      self.y3   = y3
      self.y4   = y4
   def eval(self, value, x):
      ff = ((x[1]-self.y3)*(x[1]-self.y4))**2
      if ff < DOLFIN_EPS:
         value[0] = 0.0
         value[2] = 0.0
      elif x[1]>self.y3 and x[1]<self.y4:
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
# Returns symmetric part of strain tensor
def epsilon(u):
   return 0.5*(nabla_grad(u) + nabla_grad(u).T)
#-----------------------------------------------------------------------------
def nonlinear_form(Re,Gr,Pr,hf,ds,up,vp):
   nu    = 1/Re
   k     = 1/(Re*Pr)
   G     = Gr/Re**2
   hs    = HeatSource()
   # Trial function
   u  = as_vector((up[0],up[1]))   # velocity
   T  = up[2]                      # temperature
   p  = up[3]                      # pressure
   tau= 2*nu*epsilon(u)
   # Test function
   v  = as_vector((vp[0],vp[1]))   # velocity
   S  = vp[2]                      # temperature
   q  = vp[3]                      # pressure
   Fns   =   inner(grad(u)*u, v)*dx      \
          + inner(tau,epsilon(v))*dx     \
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
   tau   = 2*nu*epsilon(u)
   Aform = - inner( grad(u)*us, v )*dx \
           - inner( grad(us)*u, v )*dx \
           + div(v)*p*dx \
           - inner( tau, epsilon(v) )*dx \
           + G*T*v[1]*dx \
           + div(u)*q*dx \
           - inner(grad(Ts),u)*S*dx \
           - inner(grad(T),us)*S*dx \
           - k*inner(grad(T),grad(S))*dx
   return Aform
#-----------------------------------------------------------------------------
class NSProblem():
   def __init__(self, udeg, Re, Gr, Pr, y3, y4):
      mesh = Mesh("square.xml")
      sub_domains = create_subdomains(mesh,x1,x2,y1,y2,y3,y4)

      self.udeg = udeg
      self.tdeg = udeg
      self.pdeg = udeg - 1
      self.Re = Re
      self.Gr = Gr
      self.Pr = Pr
      self.ds = Measure("ds")[sub_domains]
      self.y3 = y3
      self.y4 = y4

      self.V = VectorFunctionSpace(mesh, "CG", self.udeg)  # velocity
      self.W = FunctionSpace(mesh, "CG", self.tdeg)        # temperature
      self.Q = FunctionSpace(mesh, "CG", self.pdeg)        # pressure
      self.X = MixedFunctionSpace([self.V, self.W, self.Q])

      print 'Number of degrees of freedom = ', self.X.dim()

      # Velocity bc
      noslipbc1 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 0)
      noslipbc2 = DirichletBC(self.X.sub(0), (0,0), sub_domains, 3)
      # velocity control boundary
      self.gs   = velocity(0.0,self.y3,self.y4)
      vconbc    = DirichletBC(self.X.sub(0).sub(0), self.gs, sub_domains, 2)
      yvbc      = DirichletBC(self.X.sub(0).sub(1), 0,  sub_domains, 2)

      # Temperature bc
      tbc1    = DirichletBC(self.X.sub(1), 0.0, sub_domains, 0)
      self.ts = temperature(0.0,self.y3,self.y4)
      tbc2    = DirichletBC(self.X.sub(1), self.ts, sub_domains, 2)

      self.bc = [noslipbc1, noslipbc2, vconbc, yvbc, tbc1, tbc2]
      self.vbc= vconbc
      self.tbc= tbc2
      self.hf = HeatFlux(0.0)

   # solves steady state equations
   def steady_state(self, Relist):
      up = Function(self.X)
      vp = TestFunction(self.X)
      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hf,self.ds,up,vp)

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
      u.rename("v","velocity"); T.rename("T","temperature"); p.rename("p","pressure");
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
      freeinds = np.setdiff1d(range(N),bcinds,assume_unique=True).astype(np.int32)

      # pressure indices
      pinds = self.X.sub(2).dofmap().dofs()

      return freeinds, pinds

   # Generate linear state representation
   # Also compute k eigenvalues/vectors
   def linear_system(self, k=0):
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

      print 'size of pinds =', len(pinds)
      print 'size of vinds =', len(vinds)
      print 'size of tinds =', len(tinds)

      # mass matrix
      M = Ma[freeinds,:][:,freeinds]
      print "Size of M =",M.shape

      # stiffness matrix
      A = Aa[freeinds,:][:,freeinds]
      print "Size of A =",A.shape

      # velocity control operator
      ua = interpolate(allvar(1,1,self.y3,self.y4), self.X)
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

      if k>0:
        # Compute eigenvalues/vectors of (A,M)
        print "Computing eigenvalues/vectors ..."
        sigma = -1.0
        vals, vecs = la.eigs(A, k=k, M=M, sigma=sigma, which='LR')
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

        # Compute eigenvalues/vectors of (A^T,M^T)
        # First transpose A; M is symmetric
        A.transpose()
        vals, vecs = la.eigs(A, k=k, M=M, sigma=sigma, which='LR')
        for val in vals:
            print np.real(val), np.imag(val)
      
        for e in range(0,k,2):
            filename = "evec"+str(e+1)+"a"
            print "Writing into file ", filename
            # Save real part of eigenvector. << outputs only real part
            ua.vector()[freeinds] = vecs[:,e]
            File(filename+".xml") << ua.vector()
            u,T,p = ua.split()
            File(filename+"_u.pvd") << u
            File(filename+"_T.pvd") << T
            File(filename+"_p.pvd") << p

            filename = "evec"+str(e+2)+"a"
            print "Writing into file ", filename
            # Save imaginary part of eigenvector. << outputs only real part
            ua.vector()[freeinds] = vecs[:,e] * (-1j)
            File(filename+".xml") << ua.vector()
            u,T,p = ua.split()
            File(filename+"_u.pvd") << u
            File(filename+"_T.pvd") << T
            File(filename+"_p.pvd") << p

   # Compute controllability term
   def ctrb(self,m):
      eigvec = Function(self.X)
      u = as_vector((eigvec[0], eigvec[1]))
      T = eigvec[2]
      p = eigvec[3]
      n = UnitNormal()
      nu  = 1.0/self.Re
      tau = 2*nu*epsilon(u)

      for k in range(1,m+1):
         f1 = "evec"+str(k)+"a.xml"
         # first eigenvector
         File(f1) >> eigvec.vector()
         # velocity control
         sigma = project(-p*n + tau*n, self.V)
         f2 = "ctrb"+str(k)+"_u.pvd"
         File(f2) << sigma

         # temperature control
         #hf = project(inner(grad(T),n), self.Q)
         hf = project(inner(grad(T),n)*n, self.V)
         f3 = "ctrb"+str(k)+"_T.pvd"
         File(f3) << hf

         # heat flux control
         t = project(T*n, self.V)
         f4 = "ctrb"+str(k)+"_h.pvd"
         File(f4) << t
         print "Wrote files ", f2, f3, f4

   # Runs nonlinear model
   def run(self,with_control=False):
      up = Function(self.X)
      vp = TestFunction(self.X)

      Re = Constant(self.Re)
      Gr = Constant(self.Gr)
      Pr = Constant(self.Pr)
      F = nonlinear_form(Re,Gr,Pr,self.hf,self.ds,up,vp)

      fcont = open('control.dat','w')
      if with_control:
         # compute indices of velocity and temperature
         freeinds,pinds = self.get_indices()
         vTinds = np.setdiff1d(freeinds,pinds,assume_unique=True).astype(np.int32)
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
      up1.vector()[:] += 0.01 * uppert.vector()

      fu = File("u.pvd")
      ft = File("T.pvd")

      uppert.vector()[:] = up1.vector() - ups.vector()
      u,T,p = uppert.split()
      fu << u
      ft << T
      KE  = assemble(0.5*inner(u,u)*dx)
      dKE = assemble(0.5*inner(u,u)*dx)
      dHE = assemble(0.5*inner(T,T)*dx)
      print 'Kinetic energy =', KEs, KE, dKE, dHE
      fhist.write(str(0)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+" "+str(dHE)+"\n")

      dt = 0.1
      final_time = dt*2000
      time, iter = 0, 0

      if with_control:
        dy = up1.vector().array() - ups.vector().array()
        a = -f(time)*np.dot(gain['Kt'], dy[vTinds])
        self.gs.amp = a[0]
        self.ts.amp = a[1]
        self.hf.amp = a[2]
        fcont.write(str(time)+" "+str(a[0])+" "+str(a[1])+" "+str(a[2])+"\n")

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

      uppert.vector()[:] = up.vector() - ups.vector()
      u,T,p = uppert.split()
      fu << u
      ft << T
      KE  = assemble(0.5*inner(u,u)*dx)
      dKE = assemble(0.5*inner(u,u)*dx)
      dHE = assemble(0.5*inner(T,T)*dx)
      print 'Kinetic energy =', KEs, KE, dKE, dHE
      fhist.write(str(time)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+" "+str(dHE)+"\n")
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
         # initial guess by extrapolation
         up.vector()[:] = 2 * up1.vector() - up2.vector()
         if with_control:
            dy = up.vector().array() - ups.vector().array()
            a = -f(time)*np.dot(gain['Kt'], dy[vTinds])
            self.gs.amp = a[0]
            self.ts.amp = a[1]
            self.hf.amp = a[2]
            fcont.write(str(time)+" "+str(a[0])+" "+str(a[1])+" "+str(a[2])+"\n")
            fcont.flush()
         solver2.solve()
         iter += 1
         time += dt
         print 'Iter = {:5d}, t = {:f}'.format(iter, time)
         uppert.vector()[:] = up.vector() - ups.vector()
         u,T,p = uppert.split()
         if iter%10 == 0:
            fu << u
            ft << T
         KE  = assemble(0.5*inner(u,u)*dx)
         dKE = assemble(0.5*inner(u,u)*dx)
         dHE = assemble(0.5*inner(T,T)*dx)
         print 'Kinetic energy =', KEs, KE, dKE, dHE
         fhist.write(str(time)+" "+str(KEs)+" "+str(KE)+" "+str(dKE)+" "+str(dHE)+"\n")
         fhist.flush()
         print '--------------------------------------------------------------'

      fhist.close()
