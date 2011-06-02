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
   # State solution
   u = x1 + math.sqrt(x1**2 + a(x2))
   # Write primal solution to file
   fout = open("primal.dat", "w")
   fout.write("%.15e" % u)
   fout.write("\n")
   fout.close()
elif mode == "2":
   # Read primal from file
   fin = open("primal.dat", "r")
   u = float(fin.readline())
   fin.close()
   # Residual
   R = 0.5*u**2 - x1*u - 0.5*a(x2)
   # Write primal residual to file
   fout = open("p_residual.dat", "w")
   fout.write("%.15e" % R)
   fout.write("\n")
   fout.close()
   # Read adjoint from file
   fin = open("adjoint.dat", "r")
   v = float(fin.readline())
   fin.close()
   # Adjoint correction
   VdotR = v * R
   # Write adjoint correcion to file
   fout = open("VdotR.dat", "w")
   fout.write("%.15e" % VdotR)
   fout.write("\n")
   fout.close()

# Functional
J = u**2

# Write function to file
fout = open("obj.dat", "w")
fout.write("%.15e" % J)
fout.write("\n")
fout.close()
