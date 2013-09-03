Mesh.RecombineAll=1;           //recombine all defined surfaces
Mesh.Algorithm=8;              //delquad mesher
Mesh.RecombinationAlgorithm=1; //blossom

r = 0.4;
R = 5;

h1 = 0.02;
h2 = 0.01;
h3 = 0.1;  //farfield

Point(1) = {0, 0, 0,h1};
Point(2) = {R,0,0,h3};
Point(3) = {0,R,0,h3};

Point(4) = {1-r, 0, 0, h2};
Point(5) = {1,0,0, h2};
Point(6) = {1+r,0,0, h2};


Circle(1) = {2,1,3};
Line(2) = {3,1};
Line(3) = {1, 4};
Circle(4) = {6,5,4};
Line(5) = {6,2};

Line(6) = {4,5};
Line(7) = {5,6};

Line Loop(1) = {1,2,3,-4,5};
Line Loop(2) = {6,7,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Surface(1) = {1};
Physical Surface(2) = {2};

Physical Line(3) = {3,6,7,5}; // neumann
Physical Line(4) = {2}; // zero bc
Physical Line(5) = {1}; // farfield
