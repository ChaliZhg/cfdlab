from ns import *
from param import *

problem = NSProblem(udeg, Re, Gr, Pr, y3, y4)
problem.steady_state([50,60,70,80,85,90,95,100])
