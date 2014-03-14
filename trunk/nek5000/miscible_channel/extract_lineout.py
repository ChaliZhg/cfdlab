'''
Interpolates solution along a line
'''

DB = "avgchannel.nek5000"
Mesh = "mesh"

# End points of the line
p0 = (7.0, 0.0, 0)
p1 = (7.0, 1.0, 0)

nPoints = 200   # No. of sample points

OpenDatabase(DB)
AddPlot("Pseudocolor", "s1")
AddPlot("Pseudocolor", "x_velocity")
lastDB = TimeSliderGetNStates() - 1
SetTimeSliderState(lastDB)
DrawPlots()
DefineScalarExpression("xc", "coord(%s)[0]" % Mesh)
DefineScalarExpression("yc", "coord(%s)[1]" % Mesh)
 
# Do a lineout on all 4 variables to produce 4 curves.
Lineout(p0, p1, ("s1", "x_velocity", "xc", "yc"), nPoints)

# Get the data
SetActiveWindow(2)

SetActivePlots(0)
s1 = GetPlotInformation()["Curve"]

SetActivePlots(1)
ux = GetPlotInformation()["Curve"]

SetActivePlots(2)
xc = GetPlotInformation()["Curve"]

SetActivePlots(3)
yc = GetPlotInformation()["Curve"]

f=open('data.txt','w')
f.write("#x  y  s1  ux\n")
for i in range(len(s1)/2):
   idx = 2*i + 1
   f.write("%g\t%g\t%g\t%g\n" % (xc[idx], yc[idx], s1[idx], ux[idx]))

f.close()
quit()
