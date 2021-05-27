// all measurements in meters
h = 5.;
xBottom = 47.;
yBottom = 102.;
zBottom = 12.;

ySlopeTopLeft = 14.;
ySlopeTopRight = yBottom - 14.;
ySlopeMidLeft = 14.+7.;
ySlopeMidRight = yBottom - 14. - 7.;

xMiddle = 23.;
zMiddle = zBottom + 9.;
xTop = 18.;
//layerBottom = zB/h;

//Bottom part.
Point(10) = {0, 0, 0, h};
Point(11) = {xBottom, 0, 0, h};
Point(12) = {xBottom, yBottom, 0, h};
Point(13) = {0, yBottom, 0, h};

Line(100) = {10, 11};
Line(101) = {11, 12};
Line(102) = {12, 13};
Line(103) = {13, 10};

Line Loop(1000) = {100,101,102,103};
Plane Surface(2000) = {1000};

//Left, front view
Point(14) = {xBottom, 0, zBottom, h};
Point(15) = {xTop, 0, zMiddle, h};
Point(16) = {0, 0, zMiddle, h};
Line(104) = {11, 14};
Line(105) = {14, 15};
Line(106) = {15, 16};
Line(107) = {16, 10};
Line Loop(1001) = {100,104,105,106,107};
Plane Surface(2001) = {1001};

//Right, front view
Point(18) = {xBottom, yBottom, zBottom, h};
Point(19) = {xTop, yBottom, zMiddle, h};
Point(20) = {0, yBottom, zMiddle, h};
Line(109) = {12, 18};
Line(110) = {18, 19};
Line(111) = {19, 20};
Line(112) = {20, 13};
Line Loop(1002) = {-102,109,110,111,112};
Plane Surface(2002) = {1002};




//some intermediate points and lines
// for front parts
Point(22) = {xBottom, ySlopeMidLeft, zBottom, h};
Point(24) = {xBottom, ySlopeMidRight, zBottom, h};
Line(119) = {22, 14};
Line(122) = {24, 18};
Line(127) = {22, 24};

// for top parts
Point(21) = {xTop, ySlopeTopLeft, zMiddle, h};
Point(23) = {xTop, ySlopeTopRight, zMiddle, h};
Line(117) = {15, 21};
Line(120) = {19, 23};
Line(129) = {23, 21};

//Front
//Line(114) = {18, 14};
Line Loop(1003) = {101,109,-122,-127,119,-104};
Plane Surface(2003) = {1003};

//Back
Line(115) = {20, 16};
Line Loop(1004) = {-103,-112,115,107};
Plane Surface(2004) = {1004};

//Top
//Line(116) = {15, 19};
Line Loop(1005) = {-106,117,-129,-120,111,115};
Plane Surface(2005) = {1005};

//Left slope
Line(118) = {21, 22};

Line Loop(1006) = {105,117,118,119};
Plane Surface(2006) = {1006};

//Right slope
Line(121) = {23, 24};
Line Loop(1007) = {110,120,121,122};
Plane Surface(2007) = {1007};

//Left slope inner
Point(25) = {xMiddle, ySlopeMidLeft, zBottom, h};
Line(123) = {22, 25};
Line(124) = {25, 21};
Line Loop(1008) = {118,123,124};
Plane Surface(2008) = {1008};

//Right slope inner
Point(26) = {xMiddle, ySlopeMidRight, zBottom, h};
Line(125) = {24, 26};
Line(126) = {26, 23};
Line Loop(1009) = {121,125,126};
Plane Surface(2009) = {1009};

//Slope bottom
Line(128) = {26, 25};
Line Loop(1010) = {127,125,128,-123};
Plane Surface(2010) = {1010};

//Slope
Line Loop(1011) = {-128,126,129,-124};
Plane Surface(2011) = {1011};

Surface Loop(10000) = {2000, 2001, 2006, 2008, 2011, 2010, 2003, 2009, 2007, 2002, 2005, 2004};

Volume(20000) = {10000};
Physical Surface(1) = {2000,2001,2002,2003,2004};
Physical Surface(2) = {2005,2006,2007,2008,2009,2010,2011};
Physical Volume(10) = {20000};

Mesh.Algorithm3D = 4;
