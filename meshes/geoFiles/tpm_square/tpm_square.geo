//Square
h = .1;
L = 1.;
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

Lp2 = 0.4;
pp2 = 0.3;
Point(6) = {pp2, L - Lp2, 0, h};
Point(7) = {pp2 + Lp2, L - Lp2, 0, h};
Point(8) = {pp2 + Lp2, L, 0, h};
Point(9) = {pp2, L, 0, h};

Line(10) = {1, 2}; 
Line(11) = {2, 3}; 
Line(12) = {3, 8}; 
Line(13) = {8, 7}; 
Line(14) = {7, 6}; 
Line(15) = {6, 9};
Line(16) = {9, 4};
Line(17) = {4, 1};

Line Loop(100) = {16,17,10,11,12,13,14,15}; 

Plane Surface(101) = {100};

Line(18) = {9, 8};
Line Loop(102) = {13,14,15,18}; 

Plane Surface(103) = {102};


Physical Line(1) = {10};
Physical Line(2) = {17};
Physical Line(3) = {11};
Physical Line(4) = {16};
Physical Line(6) = {20};
Physical Line(5) = {12};


Physical Surface(10) = {101, 103};
