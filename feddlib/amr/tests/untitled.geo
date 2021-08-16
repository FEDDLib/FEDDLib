//+
Point(1) = {-1, 1, 0, 1.0};
//+
Point(2) = {-1, 0, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {0, -1, 0, 1.0};
//+
Point(5) = {4, -1, 0, 1.0};
//+
Point(6) = {4, 1, 0, 1.0};
//+
Point(7) = {4, 1, 1, 1.0};
//+
Point(8) = {4, -1, 1, 1.0};
//+
Point(9) = {0, -1, 1, 1.0};
//+
Point(10) = {0, 0, 1, 1.0};
//+
Point(11) = {-1, 0, 1, 1.0};
//+
Point(12) = {-1, 1, 1, 1.0};
//+
Line(1) = {12, 11};
//+
Line(2) = {11, 10};
//+
Line(3) = {10, 9};
//+
Line(4) = {9, 8};
//+
Line(5) = {8, 7};
//+
Line(6) = {7, 12};
//+
Line(7) = {12, 1};
//+
Line(8) = {1, 6};
//+
Line(9) = {6, 5};
//+
Line(10) = {5, 8};
//+
Line(11) = {7, 6};
//+
Line(12) = {4, 9};
//+
Line(13) = {4, 3};
//+
Line(14) = {3, 10};
//+
Line(15) = {3, 2};
//+
Line(16) = {2, 11};
//+
Line(17) = {2, 1};
//+
Line(18) = {4, 5};
//+
Curve Loop(1) = {8, -11, 6, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 10, 5, 11};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -4, -12, 18};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, -12, 13, 14};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, -2, -16, -15};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {16, -1, 7, -17};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {8, 9, -18, 13, 15, 17};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {5, 6, 1, 2, 3, 4};
//+
Plane Surface(8) = {8};
//+
Surface Loop(1) = {8, 2, 7, 1, 6, 5, 4, 3};
//+
Volume(1) = {1};
//+
Physical Surface("1") = {8, 4, 3, 7, 1, 5};
//+
Physical Surface("2") = {6};
//+
Physical Surface("3") = {2};
//+
Physical Volume("10") = {1};
