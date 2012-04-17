"""
Test grid adaptation
"""
from dolfin import *
import numpy

# Initial mesh
n = 10
mesh = UnitSquare(n,n)

flag = CellFunction("bool", mesh)
for c in cells(mesh):
   if c.index() == 50:
      flag[c] = True
   else:
      flag[c] = False

mesh_new = refine(mesh, flag)
plot(mesh)
plot(mesh_new)
interactive()
