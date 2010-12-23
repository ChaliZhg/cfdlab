#!/usr/bin/python

import sys
import math

mode = sys.argv[1]

fin = open("solver.in", "r")
x1 = float(fin.readline())
x2 = float(fin.readline())
fin.close()

if mode == "1":
   # State solution
   u = x1 + math.sqrt(x1**2 + math.exp(-50*(x2-0.5)**2))
   # Write primal solution to file
   fout = open("primal.dat", "w")
   fout.write(str(u))
   fout.write("\n")
   fout.close()
elif mode == "2":
   # Read primal from file
   fin = open("primal.dat", "r")
   u = float(fin.readline())
   fin.close()
   # Residual
   R = 0.5*u**2 - x1*u - 0.5*math.exp(-50*(x2-0.5)**2)
   # Write primal residual to file
   fout = open("p_residual.dat", "w")
   fout.write(str(R))
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
   fout.write(str(VdotR))
   fout.write("\n")
   fout.close()

# Functional
J = u**2

# Write function to file
fout = open("obj.dat", "w")
fout.write(str(J))
fout.write("\n")
fout.close()
