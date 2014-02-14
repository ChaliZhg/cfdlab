"""
Mark boundaries
"""

from dolfin import *

set_log_level(1)

# These values must be same as those in cylinder_in_channel.geo file
x1 = 0.4
x2 = 0.6

y1 = 0.1
y2 = 0.4

y3 = 0.7
y4 = 0.9

n  = 50

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
class OutFlow(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],0.0) and between(x[1], (y1, y2)) and on_boundary

# Sub domain for inflow (right)
class InFlow(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],1.0) and between(x[1],(y3,y4)) and on_boundary

# Sub domain for inflow (right)
class Heat(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],0.0) and between(x[0],(x1,x2)) and on_boundary

# Read mesh
mesh = UnitSquareMesh(n,n)

# Create mesh functions over the cell facets
sub_domains = FacetFunction("size_t", mesh)

# Mark all facets as sub domain 4, includes interior edges
sub_domains.set_all(4)

# Mark all boundary edges
boundary = Boundary()
boundary.mark(sub_domains, 0)

# Mark outflow
outflow = OutFlow()
outflow.mark(sub_domains, 1)

# Mark control area
inflow = InFlow()
inflow.mark(sub_domains, 2)

# Mark heating boundary
heat = Heat()
heat.mark(sub_domains, 3)

# Save mesh
File("square.xml") << mesh

# Save sub domains to file
File("subdomains.xml") << sub_domains

# Save sub domains to VTK files for visualization
File("subdomains.pvd") << sub_domains
