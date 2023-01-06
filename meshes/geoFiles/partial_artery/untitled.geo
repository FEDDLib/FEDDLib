n_elem = 32;
n_elem_solid = 30;

nb_layers = 45;

L=1.5;

Point(0) = {1, 1, 1.5, 1.0};
Point(1) = {1, 0, 0, 1.0};
Point(2) = {1, 0, 1.5, 1.0};
Point(3) = {0, 1, 1.5, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {-0.25, 1, 0, 1.0};
Point(6) = {-0.25, 1, 1.5, 1.0};
Point(7) = {1, -0.25, 1.5, 1.0};
Point(8) = {1, -0.25, 0, 1.0};
Point(9) = {1,1,0,1.0};

Circle(1) = {2,0,3};
Circle(2) = {4,9,1};
Circle(3) = {6,0,7};
Circle(4) = {8,9,5};

Line(5) = {3, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line(8) = {6, 3};
Line(9) = {2, 7};
Line(10) = {7, 8};
Line(11) = {8, 1};
Line(12) = {1, 2};

Curve Loop(1) = {6, 7, 8, 5};
Curve Loop(2) = {2, 12, 1, 5};
Curve Loop(3) = {6, -4, 11, -2};
Curve Loop(4) = {10, 11, 12, 9};
Curve Loop(5) = {3, -9, 1, -8};
Curve Loop(6) = {4, 7, 3, 10};


Surface(1) = {4};
Surface(2) = {1};
Surface(3) = {3};
Surface(4) = {5};
Surface(5) = {6};
Surface(6) = {2};

Surface Loop(1) = {1, 2, 5, 6, 4, 3};
Volume(1) = {1};

// ############ TRANSFINITE LINES ############
Transfinite Line{1,2}=n_elem+1; //Number of points on a 1/4 of inner circle length: 1.5708
Transfinite Line{3,4}=n_elem*1.25; //Number of points on a 1/4 of outer circle: 1.9634
Transfinite Line{6,8,9,11}=n_elem_solid/6+1; //Number of points on a 1/4 of circle
Transfinite Line{5,7,10,12}=n_elem_solid+1;      //Number of structure layers+1 

Physical Surface("1", 13) = {3};
Physical Surface("2", 14) = {5};
Physical Surface("3", 15) = {4, 2};
Physical Surface("4", 16) = {1};
Physical Surface("5", 17) = {6};
Physical Curve("6", 18) = {1, 3, 4, 2};
Physical Curve("7", 19) = {9, 11};
Physical Curve("8", 20) = {8, 6};
