//Square
h = .02;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

Line(5) = {1, 2}; 
Line(6) = {2, 3}; 
Line(7) = {3, 4}; 
Line(8) = {4, 1}; 

Line Loop(9) = {5,6,7,8}; 

Plane Surface(10) = 9;

Physical Line(1) = {5};
Physical Line(2) = {6};
Physical Line(3) = {7};
Physical Line(4) = {8};

Physical Surface(10) = {10};

