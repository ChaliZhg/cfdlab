from param import *
from ns import *

# How many eigenfunctions to check
k=6

problem = NSProblem(udeg, Re, Gr, Pr, y3, y4)
problem.ctrb(k)
