#!/usr/bin/python

# Read primal residual from file
fin = open("p_residual.dat", "r")
R = float(fin.readline())
fin.close()

# Read adjoint residual from file
fin = open("a_residual.dat", "r")
AR = float(fin.readline())
fin.close()

# Read primal2-primal1 from file
fin = open("dprimal.dat", "r")
du = float(fin.readline())
fin.close()

# Read adjoint2-adjoint1 from file
fin = open("dadjoint.dat", "r")
dv = float(fin.readline())
fin.close()

# Remaining error
re = (dv * R) + (du * AR)

# Write primal solution to file
fout = open("RE.dat", "w")
fout.write(str(re))
fout.write("\n")
fout.close()
