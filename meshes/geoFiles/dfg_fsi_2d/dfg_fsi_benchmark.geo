// Gmsh project created on Fri Jun 01 15:23:36 2018
// Punkte
h = 0.03; // Element-size-factor ist 0.05. Zu finden unter: Tools->Options->Mesh->General->Element size factor
hObst = h/2;
hSolid = h/3;
hOutflow = 1.5*h;
Point(1) = {0, 0, 0, h};
Point(2) = {2.5, 0, 0, hOutflow};
Point(3) = {2.5, 0.41, 0, hOutflow};
Point(4) = {0, 0.41, 0, h};
Point(5) = {0.15, 0.2, 0, hObst};
Point(6) = {0.2, 0.2, 0, hObst};
Point(7) = {0.249, 0.21, 0, hSolid};
Point(8) = {0.249, 0.19, 0, hSolid};
Point(9) = {0.6, 0.19, 0, hSolid};
Point(10) = {0.6, 0.2, 0, hSolid};
Point(11) = {0.6, 0.21, 0, hSolid};

// Alle Liniensegmente; Eingabeargument ist Punktindex
Line(1) = {1, 2}; // Ab hier Wand
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {8, 9}; // Ab hier Struktur
Line(6) = {9, 10};
Line(7) = {10, 11};
Line(8) = {11, 7};

// Alle drei Kreissegmente des Kreises/ Zylinders;
// Eingabearguemnt ist Punkteindex
Circle(9) = {7, 6, 5};
Circle(10) = {5, 6, 8};
Circle(11) = {8, 6, 7}; // Verbindungsstueck in der Struktur

// Wenn dann Physical Line

// Definiere verschiedene Teilgebiete;
// Eingabeargument ist der Line- bzw. Circle-Index
//Physical Curve("wall") = {3, 1};
//Physical Curve("inflow") = {4};
//Physical Curve("outflow") = {2};
//Physical Curve("circle") = {9, 10, 11};
//Physical Curve("structure") = {5, 6, 7, 8, 11}; 

// Definiere komplette Raender;
// Eingabeargument erneut Line- bzw. Circle-Index
Line Loop(12) = {4, 1, 2, 3}; // kompletter Kanal
Line Loop(13) = {9, 10, 11}; // Kreisscheibe
//-11 anstatt 11, da sonst falsche Orientierung.
Line Loop(14) = {5, 6, 7, 8, -11}; // Balken (Struktur); 

// Gebe nun die einzelnen Gebiete an;
// Eingabeargumente ist ist Curve Loop Index;
// 1. Argument: Komplettes Gebiet (surface)
// 2. - n. Argument: Loecher im Gebiet
Plane Surface(15) = {12, 13, 14}; // Fluid
Plane Surface(16) = {14}; // Struktur

// Siehe oben
//Physical Curve("wall") += {3, 1};
//Physical Curve("inflow") += {4};
//Physical Curve("outflow") += {2};
//Physical Curve("circle") += {9, 10, 11};
//Physical Curve("structure") += {5, 6, 7, 8, 11}; 
//Physical Surface("inner") = {15, 16};
