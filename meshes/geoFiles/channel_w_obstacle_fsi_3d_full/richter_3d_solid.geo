Merge "richter_3d_geometry.msh";

// Fuer die Nummern, siehe richter_3d_geometry.geo und dort unter visibility
// oder in der richter_3d_geometry.geo Datei.

// Bei Surface(1) wird ueberall homogener Dirichletrand gesetzt
Physical Surface(1) = {212}; // unten, d.h. homogener Dirichlet
Physical Surface(6) = {207, 208, 209, 210, 211}; // Interface, muss dieselbe Flag wie das Interface vom Fluid haben

Physical Volume(10) = {402};
