Merge "richter_3d_geometry.msh";

// Fuer die Nummern, siehe richter_3d_geometry.geo und dort unter visibility
// oder in der richter_3d_geometry.geo Datei.

// Bei Surface1 wird ueberall homogener Dirichletrand gesetzt
Physical Surface(1) = {202, 201, 206}; // wand unten (202) und oben (201) und hinten (z = 0.4; keine Struktur; 206)
Physical Surface(2) = {203}; // inflow
Physical Surface(4) = {205}; // Symmetrieachse => Dirichlet Null in Normalenrichtung
Physical Surface(5) = {204}; // outflow
Physical Surface(6) = {207, 208, 209, 210}; // Interface, muss dieselbe Flag haben wie das Interface vom Solid haben
Physical Line(3) = {9,7,8}; // Line, die Interface und an der Symmetrieachse liegt
Physical Volume(10) = {401};
