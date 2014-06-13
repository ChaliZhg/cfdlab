from ns import *
from param import *
import os

# NOTE: THIS CODE IS SPECIFIC TO 50x50 mesh
n  = 50 
nw = 5
dy = 1.0/n
print 'dy = ', dy

# yc in [nw*dy, 1-nw*dy]; arange does not include endpoint
yc = np.arange(nw*dy, 1.0-nw*dy+dy, dy)
print 'yc = ', yc

m = len(yc)
e = np.zeros(m)
for i in range(m):
    y3 = yc[i] - nw*dy
    y4 = yc[i] + nw*dy
    problem = NSProblem(udeg, Re, Gr, Pr, y3, y4)
    problem.linear_system()
    # Run gain.m
    os.system("matlab -nodesktop -nosplash < gain.m > log")
    f = open('maxeig.dat','r')
    e[i] = float(f.readline())
    f.close()
    print yc[i], e[i]

for i in range(m):
    print yc[i], e[i]
