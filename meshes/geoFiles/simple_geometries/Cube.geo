//Square
h = 0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0, 1, h};
Point(6) = {1, 0, 1, h};
Point(7) = {1, 1, 1, h};
Point(8) = {0, 1, 1, h};


Line(9) = {1, 2}; 
Line(10) = {2, 3}; 
Line(11) = {3, 4}; 
Line(12) = {4, 1}; 

Line(13) = {5, 6}; 
Line(14) = {6, 7}; 
Line(15) = {7, 8}; 
Line(16) = {8, 5}; 

Line(17) = {1, 5}; 
Line(18) = {2, 6}; 
Line(19) = {3, 7}; 
Line(20) = {4, 8}; 

Curve Loop(1) = {11, 20, -15, -19};
Curve Loop(2) = {10, 11, 12, 9};
Curve Loop(3) = {17, 13, -18, -9};
Curve Loop(4) = {13, 14, 15, 16};
Curve Loop(5) = {14, -19, -10, 18};
Curve Loop(6) = {16, -17, -12, 20};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Physical Surface("1") = {6};
Physical Surface("2") = {3};
Physical Surface("3") = {2};
Physical Surface("4") = {4};
Physical Surface("5") = {1};
Physical Surface("6") = {5};

Surface Loop(1) = {4, 3, 6, 2, 5, 1};
Volume(1) = {1};
Physical Volume("10") = {1};
