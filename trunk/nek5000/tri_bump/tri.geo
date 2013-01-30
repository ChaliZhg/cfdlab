Mesh.RecombineAll=1;           //recombine all defined surfaces
//Mesh.Algorithm=8;              //delquad mesher
//Mesh.RecombinationAlgorithm=1; //blossom

xmin=-5;
xmax=+10;
Ht=1.0;
theta=Pi/4;
Bt=Ht*Cos(theta);
ymax=3*Ht;

h1=0.2;
h2= Ht;

alpha=2*Pi - (Pi/2+theta);
x1=Bt + h1*Tan(alpha);
x2=h1*Tan(alpha);
y2=Ht+h1;

cl=0.1;

Point(1) = {xmin, 0, 0, cl};
Point(2) = {-Bt, 0, 0, cl};
Point(3) = {0, Ht, 0, cl};
Point(4) = {Bt, 0, 0, cl};
Point(5) = {xmax, 0, 0, cl};
Point(6) = {xmax, Ht, 0, cl};
Point(7) = {xmax, ymax, 0, cl};
Point(8) = {0, ymax, 0, cl};
Point(9) = {xmin, ymax, 0, cl};
Point(10) = {xmin, Ht, 0, cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,9};
Line(10) = {9,10};
Line(11) = {10,1};
Line(12) = {10,3};
Line(13) = {3,6};
Line(14) = {3,8};

n1 = 15;
n2 = 15;
n3 = 15;
n4 = 50;

r1 = 1.5;
r2 = 0.05;
r3 = 1.4;
r4 = 1.1;

Line Loop(1) = {1,2,-12,11};
Plane Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,10};

Line Loop(2) = {12,14,9,10};
Plane Surface(2) = {2};
Transfinite Surface(2) = {10,3,8,9};

Line Loop(3) = {4,5,6,-13};
Plane Surface(3) = {3};
Transfinite Surface(3) = {3,4,5,6};

Line Loop(4) = {13,7,8,-14};
Plane Surface(4) = {4};
Transfinite Surface(4) = {3,6,7,8};


Transfinite Line{-1,-12} = n1 Using Progression r1;
Transfinite Line{9} = n1 Using Progression r1;

Transfinite Line{11,2,4,6} = n2 Using Bump r2;

Transfinite Line{14,-10,7} = n3 Using Progression r3;

Transfinite Line{5,13,-8} = n4 Using Progression r4;

Physical Surface(1000000) = {1,2,3,4};
