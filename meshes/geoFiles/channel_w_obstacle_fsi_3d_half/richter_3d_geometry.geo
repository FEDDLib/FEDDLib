// Gmsh project created on Fri Jun 01 15:23:36 2018
// Punkte
h = 0.035; // Element-size-factor ist/war bei 0.05


// vorne (z = 0)
Point(1) = {0, 0, 0, 2*h};
Point(2) = {0.4, 0, 0, h};
Point(3) = {0.5, 0, 0, h};
Point(4) = {1.5, 0, 0, 3*h};
Point(5) = {1.5, 0.4, 0, 3*h};
Point(6) = {0, 0.4, 0, 2*h};
Point(7) = {0.4, 0.2, 0, h};
Point(8) = {0.5, 0.2, 0, h};

// mitte (Spiegelachse; z = 0.2)
Point(9) = {0, 0, 0.2, 2*h};
Point(10) = {1.5, 0, 0.2, 3*h};
Point(11) = {1.5, 0.4, 0.2, 3*h};
Point(12) = {0, 0.4, 0.2, 2*h};
Point(13) = {0.4, 0, 0.2, h};
Point(14) = {0.5, 0, 0.2, h};
Point(15) = {0.4, 0.2, 0.2, h};
Point(16) = {0.5, 0.2, 0.2, h};

// hinten (z = 0.4)
Point(17) = {0, 0, 0.4, 2*h};
Point(18) = {1.5, 0, 0.4, 3*h};
Point(19) = {1.5, 0.4, 0.4, 3*h};
Point(20) = {0, 0.4, 0.4, 2*h};

// Alle Liniensegmente; Eingabeargument ist Punktindex
// ########## vorne ##########
// Ab hier Wand (oder Kasten)
Line(1) = {1, 2}; // unten1
Line(2) = {2, 3}; // unten2 (Struktur-Boden)
Line(3) = {3, 4}; // unten3
Line(4) = {4, 5}; // rechts
Line(5) = {5, 6}; // oben
Line(6) = {6, 1}; // links
// Ab hier Struktur (ohne Boden)
Line(7) = {3, 8}; // rechts
Line(8) = {8, 7}; // oben
Line(9) = {7, 2}; // links

// ########## Spiegelachse (= mitte; z = 0.2) ##########
Line(10) = {9, 13}; // unten1
Line(11) = {10, 11}; // rechts
Line(12) = {11, 12}; // oben
Line(13) = {12, 9}; // links
Line(14) = {13, 14}; // unten2
Line(15) = {14, 10}; // unten3
Line(16) = {14, 16}; // innen-rechts
Line(17) = {16, 15}; // innen-oben
Line(18) = {15, 13}; // innen-links

// hinten
Line(19) = {17, 18};
Line(20) = {18, 19};
Line(21) = {19, 20};
Line(22) = {20, 17};

// von vorne zur mitte
Line(23) = {1, 9};
Line(24) = {2, 13};
Line(25) = {3, 14};
Line(26) = {4, 10};
Line(27) = {5, 11};
Line(28) = {6, 12};
Line(29) = {7, 15};
Line(30) = {8, 16};

// von mitte nach hinten
Line(31) = {9, 17};
Line(32) = {10, 18};
Line(33) = {11, 19};
Line(34) = {12, 20};


// Definiere komplette Raender;
// Eingabeargument erneut Line-Index

// vorne (= Symmetrieachse)
Line Loop(101) = {6, 1, -9, -8, -7, 3, 4, 5}; // Fluid: links-unten1-links_S-oben_S-rechts_S-unten3-rechts-oben
Line Loop(102) = {9, 2, 7, 8}; // Struktur: links-unten2-rechts-oben

// mitte; evtl. braucht man Fluid nicht
Line Loop(103) = {13, 10, -18, -17, -16, 15, 11, 12}; // Fluid: siehe oben
Line Loop(104) = {16, 17, 18, 14}; // Struktur: siehe oben => das hier ist Teil des Interfaces

// hinten
Line Loop(105) = {19, 20, 21, 22}; 

// inflow
Line Loop(106) = {23, 31, -22, -34, -28, 6};

// outflow
Line Loop(107) = {26, 32, 20, -33, -27, -4};

// oben
Line Loop(108) = {28, 34, -21, -33, -27, 5};

// unten
Line Loop(109) = {23, 31, 19, -32, -26, -3, 25, -14, -24, -1}; // Fluid
Line Loop(110) = {24, 14, -25, -2}; // Struktur => homogener Dirichletrand

// Rest der Struktur
Line Loop(111) = {29, -17, -30, 8}; // oben
Line Loop(112) = {24, -18, -29, 9}; // vorne (nahe am Inflow)
Line Loop(113) = {25, 16, -30, -7}; // hinten (nahe am Outflow)


// Gebe nun die einzelnen Gebiete an;
// Eingabeargumente ist Curve Loop Index;
// 1. Argument: Komplettes Gebiet (surface)
// 2. - n. Argument: Loecher im Gebiet
// Fluid:
Plane Surface(201) = {108}; // wand oben
Plane Surface(202) = {109}; // wand unten
Plane Surface(203) = {106}; // inflow
Plane Surface(204) = {107}; // outflow
Plane Surface(205) = {101}; // wand vorne (z = 0; Symmetrieachse)
Plane Surface(206) = {105}; // wand hinten (z = 0.4)

// Struktur bzw. Interface:
Plane Surface(207) = {111}; // interface oben
Plane Surface(208) = {112}; // interface vorne
Plane Surface(209) = {113}; // interface hinten
Plane Surface(210) = {104}; // interface seite
Plane Surface(211) = {102}; // Struktur vorne (z = 0; Symmetrieachse => Normale ist Null)
Plane Surface(212) = {110}; // Struktur unten (homogener Dirichletrand)

// Berandung der Domains:
// Visualierung der Orientierung der Surfaces via: Tools->Options->Geometry->Visibility->Normals box (also die Linke) auf 100 Stellen
// Wir orientieren alle Normalen nach aussen.
// Berandung des Fluids; 
Surface Loop(301) = {203, 206, -204, -205, 201, -202, -208, 209, -207, -210};

// Berandung der Struktur:
Surface Loop(302) = {208, -209, -211, 210, 207, -212};

Volume(401) = {301}; // Fluid
Volume(402) = {302}; // Struktur


Point(21) = {0.45, 0.15, 0.15, h};

Point{21} In Volume{402};

/*
Physical Surface("back") = {surfaceVectorFluid[0], surfaceVectorStruc[0]}; // this is commonly "front" for me, b/c z=0
Physical Surface("front") = {15, 16}; // von PlaneSurface(15) und (16); this is commonly "back" for me, b/c z=1
Physical Surface("wall") = {surfaceVectorFluid[3], surfaceVectorFluid[5]}; // oben und unten die Wand
Physical Surface("inflow") = {surfaceVectorFluid[2]};
Physical Surface("outflow") = {surfaceVectorFluid[4]};
Physical Surface("circle") = {surfaceVectorFluid[6], surfaceVectorFluid[7], surfaceVectorFluid[8]};
Physical Surface("interface") = {surfaceVectorStruc[2], surfaceVectorStruc[3], surfaceVectorStruc[4], surfaceVectorStruc[5]}; 
Physical Surface("structureDir") = {surfaceVectorStruc[6]};
Physical Volume("fluid") = {surfaceVectorFluid[1]};
Physical Volume("structure") = {surfaceVectorStruc[1]};
*/

