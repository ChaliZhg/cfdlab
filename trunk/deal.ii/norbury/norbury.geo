Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

n = 50;
c = Pi/n;

h = 0.01;
h1 = 0.05;
h2 = 0.2;

R = 5.0*2;

a0 = 0.3991;
a1 = -0.0050;
a2 = -0.0379;
a3 = 0.0059;
a4 = 0.0018;
a5 = -0.0009;
a6 = 0.0001;
a7 = 0.0001;

//Function to calculate s
Function Calculate_s
s = a0 + a1 * Cos(alpha) + a2 * Cos(2 * alpha) + a3 * Cos(3 *alpha)
		+ a4 * Cos(4 * alpha) + a5 * Cos(5 * alpha) + a6 * Cos(6 * alpha) 
		+ a7 * Cos(7 * alpha);
Return

Point(1) = {0, 0, 0,h1};
Point(2) = {R,0,0,h2};
Point(3) = {0,R,0,h2};

alpha = 0;
Call Calculate_s;
Point(4) = {1+s*Cos(alpha),s*Sin(alpha),0,h};


alpha = Pi;
Call Calculate_s;
Point(5) = {1+s*Cos(alpha),s*Sin(alpha),0,h};

//Function to draw a parametric curve.
Function DrawCurve
For i In {0:n}
	alpha = c*i;

	Call Calculate_s;
		
	newPointX = 1+s * Cos(alpha);
	newPointY = s * Sin(alpha);

	newPoint = newp;    // Automatically selects new point number
	Point(newPoint) = {newPointX, newPointY, 0, h};
EndFor

BSpline(newID) = {4, newPoint-n:newPoint, 5};
Return

newID = 1;
Call DrawCurve;

//Printf("newPoint %f", newPoint);

Line(2) = {5,4};
Line(3) = {1,5};
Line(4) = {4,2};

Circle(5) = {2,1,3};

Line(6) = {3,1};

Line Loop(1) = {2,1};
Line Loop(2) = {3,-1,4,5,6};

Plane Surface(1) = {1};
Plane Surface(2) = {2}; // A = core region

Physical Surface(1) = {2};
Physical Surface(2) = {1}; // A = core region

Physical Line(3) = {2,3,4}; // neumann
Physical Line(4) = {6}; // zero bc
Physical Line(5) = {5}; // farfield
