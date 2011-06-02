#!/usr/bin/python

import sys
import math

mode = sys.argv[1]

x1 = 0.75

def a(x):
   if x<=0.6:
      return math.sin(math.pi*x)
   else:
      return 0.5 * math.sin(math.pi*x)

fin = open("solver.in", "r")
x2 = float(fin.readline())
fin.close()

if mode == "1":
   # Read primal from file
   fin = open("primal.dat", "r")
   u = float(fin.readline())
   fin.close()
   # Compute adjoint solution: v = (dJ/du) / (dR/du), J = u**2
   v = 2*u/(x1 - u)
   # Write adjoint solution to file
   fout = open("adjoint.dat", "w")
   fout.write("%.15e" % v)
   fout.write("\n")
   fout.close()
elif mode == "2":
   # Read primal from file
   fin = open("primal.dat", "r")
   u = float(fin.readline())
   fin.close()
   # Read adjoint from file
   fin = open("adjoint.dat", "r")
   v = float(fin.readline())
   fin.close()
   # Residual = dJ/du + v * dR/du
   AR = 2*u + v*(u - x1)
   # Write adjoint residual to file
   fout = open("a_residual.dat", "w")
   fout.write("%.15e" % AR)
   fout.write("\n")
   fout.close()
