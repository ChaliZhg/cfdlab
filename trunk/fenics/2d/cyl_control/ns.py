from dolfin import *
import numpy as np
import scipy.sparse as sps
import scipy.io as sio
import scipy.sparse.linalg as la

#-----------------------------------------------------------------------------
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
#-----------------------------------------------------------------------------
# Returns symmetric part of strain tensor
def epsilon(u):
   return 0.5*(nabla_grad(u) + nabla_grad(u).T)
#-----------------------------------------------------------------------------
def nonlinear_form(Re,up,vp):
   nu    = 1/Re
   # Trial function
   u  = as_vector((up[0],up[1]))   # velocity
   p  = up[2]                      # pressure
   tau= 2*nu*epsilon(u)
   # Test function
   v  = as_vector((vp[0],vp[1]))   # velocity
   q  = vp[2]                      # pressure
   Fns   =   inner(grad(u)*u, v)*dx      \
          + inner(tau,epsilon(v))*dx     \
          - div(v)*p*dx
   Fdiv  =  - q*div(u)*dx
   return Fns + Fdiv
#-----------------------------------------------------------------------------
class NSProblem():
   def __init__(self, udeg, Re):
      mesh = Mesh("cyl.xml")
      boundaries = MeshFunction("size_t", mesh, "cyl_facet_region.xml")

      self.udeg = udeg
      self.pdeg = udeg - 1
      self.Re = Re

      self.V = VectorFunctionSpace(mesh, "CG", self.udeg)  # velocity
      self.Q = FunctionSpace(mesh, "CG", self.pdeg)        # pressure
      self.X = MixedFunctionSpace([self.V, self.Q])

      print "Number of dofs =", self.X.dim()

      # Velocity bc
      self.vin = inlet_velocity(0)
      inletbc  = DirichletBC(self.X.sub(0), self.vin, boundaries, 0)
      noslipbc = DirichletBC(self.X.sub(0), (0,0), boundaries, 2)
      # velocity control boundary
      vconbc1   = DirichletBC(self.X.sub(0), (0,0), boundaries, 3)
      vconbc2   = DirichletBC(self.X.sub(0), (0,0), boundaries, 4)

      self.bc = [noslipbc, inletbc, vconbc1, vconbc2]

   # solves steady state equations
   def steady_state(self, Relist):
      up = Function(self.X)
      vp = TestFunction(self.X)
      Re = Constant(self.Re)
      F = nonlinear_form(Re,up,vp)

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
      u,p = up.split()
      u.rename("v","velocity"); p.rename("p","pressure");
      print "Saving vtk files steady_u.pvd, steady_p.pvd"
      File("steady_u.pvd") << u
      File("steady_p.pvd") << p
