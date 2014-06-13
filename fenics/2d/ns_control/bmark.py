"""
Mark boundaries
"""

from dolfin import *

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Sub domain for outflow (left wall)
class OutFlow(SubDomain):
    def __init__(self,y1,y2):
        self.y1 = y1
        self.y2 = y2
        SubDomain.__init__(self) # call base class constructor
    def inside(self, x, on_boundary):
        return (near(x[0],0.0) and
               between(x[1], (self.y1, self.y2)) and
               on_boundary)

# Sub domain for inflow (right wall)
class InFlow(SubDomain):
    def __init__(self,y3,y4):
        self.y3 = y3
        self.y4 = y4
        SubDomain.__init__(self) # call base class constructor
    def inside(self, x, on_boundary):
        return (near(x[0],1.0) and
                between(x[1],(self.y3,self.y4)) and 
                on_boundary)

# Sub domain for heat flux (bottom wall)
class Heat(SubDomain):
    def __init__(self,x1,x2):
        self.x1 = x1
        self.x2 = x2
        SubDomain.__init__(self) # call base class constructor
    def inside(self, x, on_boundary):
        return (near(x[1],0.0) and
                between(x[0],(self.x1,self.x2)) and
                on_boundary)

def create_subdomains(mesh,x1,x2,y1,y2,y3,y4):
    # Create mesh functions over the cell facets
    sub_domains = FacetFunction("size_t", mesh)

    # Mark all facets as sub domain 4, includes interior edges
    sub_domains.set_all(4)

    # Mark all boundary edges
    boundary = Boundary()
    boundary.mark(sub_domains, 0)

    # Mark outflow
    outflow = OutFlow(y1,y2)
    outflow.mark(sub_domains, 1)

    # Mark control area
    inflow = InFlow(y3,y4)
    inflow.mark(sub_domains, 2)

    # Mark heating boundary
    heat = Heat(x1,x2)
    heat.mark(sub_domains, 3)

    # Save sub domains to file
    File("subdomains.xml") << sub_domains

    # Save sub domains to VTK files for visualization
    File("subdomains.pvd") << sub_domains

    return sub_domains
