from param import *
from ns import *

problem = NSProblem(udeg, Re, Gr, Pr, y3, y4)
problem.run(with_control=True)
