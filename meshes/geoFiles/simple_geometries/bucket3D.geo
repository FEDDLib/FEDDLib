h = .5;
L = 1.;
height = 10.;
layer = height/h;

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

Line(10) = {1, 2}; 
Line(11) = {2, 3}; 
Line(12) = {3, 4};
Line(13) = {4, 1};

Line Loop(100) = {10,11,12,13};

//outL1[] = Extrude {0, 0, height} {
//    Line{10};
//    Layers{layer};
//};

//outL2[] = Extrude {0, 0, height} {
//    Line{12};
//    Layers{layer};
//};
//Physical Line(4) = {outL2[2],outL2[3], outL1[2],outL1[3]}; //x=0 or x=1 and y=0 or y=1


Plane Surface(101) = {100};

out[] = Extrude {0, 0, height} {
  Surface{101};
  Layers{layer};
};

Physical Point(1) = {1,2,3,4};
Physical Line(2) = {11, 13}; //bottom ring x=0 or x=1
Physical Line(3) = {10, 12}; //bottom ring y=0 or y=1
Physical Surface(4) = {101}; //bottom
Physical Line(5) = {108, 109, 113, 117};// side lines
Physical Surface(6) = {out[3], out[5]}; //x=0 or x=1
Physical Surface(7) = {out[2], out[4]}; //y=0 or y=1
Physical Surface(8) = {out[0]}; //top

Physical Volume(10) = { out[1] };


