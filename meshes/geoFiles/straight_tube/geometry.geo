//-*- C++ -*-

lc = 0.15;
h = 0.08;
R= 0.40;
L = 5.0;

nb_layers = 45;
n_elem_quarter = 16;
n_elem_solid = 5;

// ############# POINTS ##############
//Base
Point(0) = {0, 0, 0, lc};
Point(1) = {R, 0, 0, lc};
Point(2) = {0, R, 0, lc};
Point(3) = {-R, 0, 0, lc};
Point(4) = {0, -R, 0, lc};
Point(5) = {R+h, 0, 0, lc};
Point(6) = {0, R+h, 0, lc};
Point(7) = {-R-h, 0, 0, lc};
Point(8) = {0, -R-h, 0, lc};


//############### LINES ##############
//Base
Circle(1) = {1,0,2};
Circle(2) = {2,0,3};
Circle(3) = {3,0,4};
Circle(4) = {4,0,1};
Circle(5) = {5,0,6};
Circle(6) = {6,0,7};
Circle(7) = {7,0,8};
Circle(8) = {8,0,5};

Line (9) = {1,5};
Line (10) = {2,6};
Line (11) = {3,7};
Line (12) = {4,8};


// ############# LINE LOOPS ############
Line Loop(1) = {-1,9,5,-10}; 
Line Loop(2) = {-2,10,6,-11};
Line Loop(3) = {-3,11,7,-12};
Line Loop(4) = {-4,12,8,-9};
Line Loop(5) = {1,2,3,4};


// ############ SURFACES ############
Plane Surface(1) = {1}; 
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};


// ############ TRANSFINITE LINES ############
Transfinite Line{1,2,3,4,5,6,7,8}=n_elem_quarter; //Number of points on a 1/4 of circle
Transfinite Line{9,10,11,12}=n_elem_solid;      //Number of structure layers+1 

// ############ TRANSFINITE SURFACES ############
Transfinite Surface { 1 } = { 1,5,6,2} ; 
Transfinite Surface { 2 } = {2,6,7,3} ;
Transfinite Surface { 3 } = {3,7,8,4} ;
Transfinite Surface { 4 } = {4,8,5,1} ;

//############# EXTRUSION TO 3D ###############
nb_extrusions = 4; // fixed, do not change

For i In {1:nb_extrusions}
 out[] = Extrude{0, 0, L}{Surface{i}; Layers{nb_layers}; }; 
 IndexWallFaceOutlet[i-1] = out[0];
 IndexInnerWall[i-1] = out[2];
 IndexOuterWall[i-1] = out[nb_extrusions];
EndFor

out[] = Extrude{0, 0, L}{Surface{5}; Layers{nb_layers}; }; 
IndexOutlet = out[0];
