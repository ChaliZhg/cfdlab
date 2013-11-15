import os
os.system("grep Scalar channel.log > scalar.dat")

lines = file('scalar.dat').readlines()
smin     = []
smax     = []
for line in lines:
   s = line.split()
   smin.append(float(s[5]))
   smax.append(float(s[6]))

print "Scalar min =", min(smin)
print "Scalar max =", max(smax)
