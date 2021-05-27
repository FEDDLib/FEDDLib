//Square
h = .25;
L = 3.;
obstBegin = 1.;
obstLength = h;
obstHeight = 0.5;

Point(1) = {0, 0, 0, h};
Point(2) = {obstBegin, 0, 0, h};
Point(3) = {obstBegin, obstHeight, 0, h};
Point(4) = {obstBegin + obstLength, obstHeight, 0, h};
Point(5) = {obstBegin + obstLength, 0, 0, h};
Point(6) = {L, 0, 0, h};
Point(7) = {L, 1, 0, h};
Point(8) = {0, 1, 0, h};

Line(10) = {1, 2}; 
Line(11) = {2, 3}; 
Line(12) = {3, 4}; 
Line(13) = {4, 5}; 
Line(14) = {5, 6}; 
Line(15) = {6, 7};
Line(16) = {7, 8};
Line(17) = {8, 1};

//Solid bottom
Line(18) = {5, 2};

//Fluid
Line Loop(101) = {10,11,12,13,14,15,16,17}; 
Plane Surface(201) = {101};
//Solid
Line Loop(102) = {11,12,13,18};
Plane Surface(202) = {102};
