L1=0.6;
L2=2.4;
h1=0.2;
h2=0.8;

cl2=0.02;
cl1=cl2/10.0;

Point(1) = {0, 0, 0, cl2};
Point(2) = {L1, 0, 0, cl1};
Point(3) = {L1, h1, 0, cl1};
Point(4) = {L1+L2, h1, 0, cl2};
Point(5) = {L1+L2, h1+h2, 0, cl2};
Point(6) = {0, h1+h2, 0, cl2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};

Physical Surface(100000) = {1};

Physical Line(100001) = {6};        // in flow
Physical Line(100002) = {1,2,3,5};  // bottom and top
Physical Line(100003) = {4};        // outlet
