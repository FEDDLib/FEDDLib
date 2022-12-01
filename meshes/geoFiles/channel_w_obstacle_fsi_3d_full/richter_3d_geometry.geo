// Gmsh project created on Fri Jun 01 15:23:36 2018
// Punkte
h = 0.01; 

// Die Richtungen beziehen sich, wenn man wie bei Richter S. 286 die Geometrie betrachtet

// Fluid vorne (z = -0.4)
Point(1) = {0, 0, -0.4, 2*h};
Point(2) = {1.5, 0, -0.4, 3*h};
Point(3) = {1.5, 0.4, -0.4, 3*h};
Point(4) = {0, 0.4, -0.4, 2*h};

// Fluid hinten (z = 0.4)
Point(5) = {0, 0, 0.4, 2*h};
Point(6) = {1.5, 0, 0.4, 3*h};
Point(7) = {1.5, 0.4, 0.4, 3*h};
Point(8) = {0, 0.4, 0.4, 2*h};

// Struktur vorne (z = -0.2)
Point(9) = {0.4, 0, -0.2, h};
Point(10) = {0.5, 0, -0.2, h};
Point(11) = {0.5, 0.2, -0.2, h};
Point(12) = {0.4, 0.2, -0.2, h};

// Struktur hinten (z = 0.2)
Point(13) = {0.4, 0, 0.2, h};
Point(14) = {0.5, 0, 0.2, h};
Point(15) = {0.5, 0.2, 0.2, h};
Point(16) = {0.4, 0.2, 0.2, h};

// Alle Liniensegmente; Eingabeargument ist Punktindex
// ########## Fluid vorne ##########
Line(1) = {1, 2}; // unten
Line(2) = {2, 3}; // rechts
Line(3) = {3, 4}; // oben
Line(4) = {4, 1}; // links

// ########## Fluid hinten ##########
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// ########## Fluid von vorne nach hinten ##########
Line(9) = {1, 5}; // links-unten
Line(10) = {2, 6}; // rechts-unten
Line(11) = {3, 7}; // rechts-oben
Line(12) = {4, 8}; // links-oben

// ########## Struktur vorne ##########
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

// ########## Struktur hinten ##########
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};

// ########## Struktur von vorne nach hinten ##########
Line(21) = {9, 13}; // links-unten
Line(22) = {10, 14}; // rechts-unten
Line(23) = {11, 15}; // rechts-oben
Line(24) = {12, 16}; // links-oben


// Definiere komplette Raender;
// Eingabeargument erneut Line-Index
// ################################
// Fluid vorne
Line Loop(101) = {1, 2, 3, 4};

// Fluid hinten
Line Loop(102) = {5, 6, 7, 8};

// Fluid oben
Line Loop(103) = {-3, 11, 7, -12};

// Fluid unten
Line Loop(104) = {1, 10, -5, -9};

// Fluid inflow (links)
Line Loop(105) = {-4, 12, 8, -9};

// Fluid outflow (rechts)
Line Loop(106) = {2, 11, -6, -10};

// Solid vorne
Line Loop(107) = {13, 14, 15, 16};

// Solid hinten
Line Loop(108) = {17, 18, 19, 20};

// Solid oben
Line Loop(109) = {-15, 23, 19, -24};

// Solid unten
Line Loop(110) = {13, 22, -17, -21};

// Solid links
Line Loop(111) = {-16, 24, 20, -21};

// Solid rechts
Line Loop(112) = {14, 23, -18, -22};

// Gebe nun die einzelnen Flaechen an;
// Eingabeargumente ist Line Loop Index.
// Erstes Argument: Gesamt; Zweites bis n-tes Argument: Was weggenommen wird
// Fluid:
Plane Surface(201) = {103}; // wand oben
Plane Surface(202) = {104, 110}; // wand unten
Plane Surface(203) = {105}; // inflow (links)
Plane Surface(204) = {106}; // outflow (rechts)
Plane Surface(205) = {101}; // wand vorne (z = -0.4)
Plane Surface(206) = {102}; // wand hinten (z = 0.4)

// Interface bzw. Solid:
Plane Surface(207) = {109}; // interface oben = solid oben
Plane Surface(208) = {107}; // interface vorne = solid vorne
Plane Surface(209) = {108}; // interface hinten = solid hinten
Plane Surface(210) = {111}; // interface links = solid links
Plane Surface(211) = {112}; // interface rechts = solid rechts
Plane Surface(212) = {110}; // Struktur unten (homogener Dirichletrand)

// Berandung der Domains:
// Die Normalen werden ausgehend von der ersten Line aus der vom Plane Surface zugehoerigen Line Loop immer nach rechts orientiert.
// Visualierung der Orientierung der Surfaces via: Tools->Options->Geometry->Visibility->Normals box (also die Linke) auf 100 Stellen
// Wir orientieren alle Normalen nach aussen.
// Berandung des Fluids; 
// Ersten 6 sind Fluid Rand und dann noch Struktur
Surface Loop(301) = {-201, 202, -203, 204, 205, -206, 207, -208, 209, 210, -211};

// Berandung der Struktur:
Surface Loop(302) = {-207, 208, -209, -210, 211, 212};

// Erstes Argument: Gesamt; Zweites bis n-tes Argument: Was weggenommen wird
Volume(401) = {301}; // Fluid
Volume(402) = {302}; // Solid

// Auswertungspunkt A = (0.45, 0.15, 0.15) in die Geometrie der Struktur aufnehmen
Point(17) = {0.45, 0.15, 0.15, h};
Point{17} In Volume{402};

















