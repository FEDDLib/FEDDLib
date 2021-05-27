Merge "richter_3d_geometry.msh";

// Fuer die Nummern, siehe richter_3d_geometry.geo und dort unter visibility
// oder in der richter_3d_geometry.geo Datei.

// Bei Surface(1) wird ueberall homogener Dirichletrand gesetzt
Physical Surface(1) = {201, 202, 205, 206}; // Wand: oben, unten, vorne, hinten
Physical Surface(2) = {203}; // inflow
Physical Surface(5) = {204}; // outflow
Physical Surface(6) = {207, 208, 209, 210, 211}; // Interface, muss dieselbe Flag haben wie das Interface vom Solid haben

Physical Volume(10) = {401};
