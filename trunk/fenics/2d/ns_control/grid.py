"""
Mark boundaries
"""

from dolfin import *

set_log_level(1)

n  = 50

# Read mesh
mesh = UnitSquareMesh(n,n)
print 'Number of cells    =', mesh.num_cells()
print 'Number of vertices =', mesh.num_vertices()

# Save mesh
File("square.xml") << mesh
