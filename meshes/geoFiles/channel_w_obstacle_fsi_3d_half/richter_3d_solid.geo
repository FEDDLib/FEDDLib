Merge "richter_3d_geometry.msh";

Physical Surface(1) = {212}; // unten, d.h. homogener Dirichlet
Physical Surface(4) = {211}; // an der Symmetrieachse, d.h. Normale zur Flaeche ist Null
Physical Surface(6) = {207, 208, 209, 210}; // Interface, muss dieselbe Flag wie das Interface vom Fluid haben
Physical Line(3) = {9,7,8}; // Line, die Interface und an der Symmetrieachse liegt
Physical Volume(10) = {402};
