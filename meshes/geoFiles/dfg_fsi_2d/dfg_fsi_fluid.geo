Merge "dfg_fsi_benchmark.msh";

Physical Line(1) = {3, 1}; // wand (3,1) + circle (=obstacle) (9,10)  => Homogener Dirichletrand
Physical Line(2) = {4}; // inflow
Physical Line(3) = {2}; // outflow
Physical Line(4) = {9, 10}; // obstacle
Physical Line(5) = {5, 6, 7, 8}; // Interface, muss dieselbe Flag haben wie das Interface vom Solid haben

Physical Surface(10) = {15};
