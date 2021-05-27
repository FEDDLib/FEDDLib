//Cube
h = .5;
L = 1.;
layer = 2;

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

outOuter[] = Extrude {0, 0, 1} {
  Surface{101};
  Layers{layer};
};

outInner[] = Extrude {0, 0, 1} {
    Surface{103};
    Layers{layer};
};

//Physical Volume(10) = {1,2};

//Physical Line(1) = {10,17,11};
//Physical Line(4) = {16};
//Physical Line(6) = {20};
//Physical Line(5) = {12};

//Physical Surface(1) = {101, 103};
//Physical Surface(1) = {outOuter[0]};//front of outer
//Physical Surface(2) = {outOuter[2]};//top left of outer
//Physical Surface(3) = {outOuter[3]};//left of outer
//Physical Surface(4) = {outOuter[4]};//bottom of outer
//Physical Surface(5) = {outOuter[5]};//right of outer
//Physical Surface(6) = {outOuter[6]};//top right of outer
//Physical Surface(7) = {outOuter[7]};//right inner of outer
//Physical Surface(8) = {outOuter[8]};//bottom inner of outer
//Physical Surface(9) = {outOuter[9]};//left inner of outer

//Physical Surface(11) = {outInner[0]};//back of inner
//Physical Surface(12) = {outInner[2]};//right of inner
//Physical Surface(13) = {outInner[3]};//bottom of inner
//Physical Surface(14) = {outInner[4]};//left of inner
//Physical Surface(15) = {outInner[5]};//top of inner

Physical Surface(1) = {101, 103, outOuter[0], outInner[0], outOuter[3], outOuter[4], outOuter[5]};
Physical Surface(4) = {outOuter[2]};//top left of outer
Physical Surface(5) = {outInner[5]};//top of inner
Physical Surface(6) = {outOuter[6]};//top right of outer

Physical Volume(10) = { outOuter[1], outInner[1] };


