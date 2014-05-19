Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

r = 0.05;
theta1 = 75.0*Pi/180.0;
theta2 = 65.0*Pi/180.0;

xmin = -1.5;
xmax = 2.2;
ymin = 0;
ymax = 0.4;

xc = 0.25;
yc = 0.2;

h1 = 0.01;
h2 = 2*Pi*r/100;
h3 = h2/4.0;

Point(1) = {xmin, ymin, 0, h1};
Point(2) = {xmax, ymin, 0, h1};
Point(3) = {xmax, ymax, 0, h1};
Point(4) = {xmin, ymax, 0, h1};
Point(5) = {xc-r, yc, 0, h2};
Point(6) = {xc, yc+r, 0, h2};
Point(7) = {xc+r, yc, 0, h2};
Point(8) = {xc, yc-r, 0, h2};
Point(9) = {xc+r*Cos(theta1),yc+r*Sin(theta1),0, h3};
Point(10) = {xc+r*Cos(theta2),yc+r*Sin(theta2),0, h3};
Point(11) = {xc+r*Cos(theta2),yc-r*Sin(theta2),0, h3};
Point(12) = {xc+r*Cos(theta1),yc-r*Sin(theta1),0, h3};

Point(100) = {xc, yc, 0, h2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {5, 100, 6};
Circle(6) = {6, 100, 9};
Circle(7) = {9, 100, 10};
Circle(8) = {10, 100, 7};
Circle(9) = {7, 100, 11};
Circle(10) = {11, 100, 12};
Circle(11) = {12, 100, 8};
Circle(12) = {8, 100, 5};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8,9,10,11,12};
Plane Surface(1) = {1,2};

Physical Line(1) = {4};   // inlet
Physical Line(2) = {2};   // outlet
Physical Line(3) = {1,3,5,6,8,9,11,12}; // bottom, top, cylinder
Physical Line(4) = {7};   // top slot
Physical Line(5) = {10};   // bottom slot

Physical Surface(1) = {1};
