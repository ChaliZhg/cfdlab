from ns import *
from param import *
import sys

# k = number of eigenvalues/vectors to compute
k = 0
if len(sys.argv) == 2:
   k = int(sys.argv[1])

problem = NSProblem(udeg, Re, Gr, Pr)
problem.linear_system(k=k)
