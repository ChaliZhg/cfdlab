delta=10; // in degrees
x0=0;
x1=0.5;
x2=1.5;
H=1;

h = 0.01;

y2 = (x2-x1)*Tan(delta*Pi/180);
Point(1) = {x0,0,0,h};
Point(2) = {x1,0,0,h};
Point(3) = {x2,y2,0,h};
Point(4) = {x2,H,0,h};
Point(5) = {x0,H,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};

Line Loop(1) = {1,2,3,4,5};
Plane Surface(1) = {1};

Physical Surface(100000) = {1};

Physical Line(100001) = {1,2}; // bottom
Physical Line(100002) = {3};   // outlet
Physical Line(100003) = {4};   // top
Physical Line(100004) = {5};   // inlet
