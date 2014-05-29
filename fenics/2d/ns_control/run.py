from param import *
from ns import *

problem = NSProblem(udeg, Re, Gr, Pr)
problem.run(with_control=True)
