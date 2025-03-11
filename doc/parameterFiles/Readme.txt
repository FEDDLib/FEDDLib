**** Description of parameter Files used inside FEDDLib ****
@2024

There are several parameter files used in the FEDDLib that play a crucial role in defining and solving physical problems. These parameter files allow the user to customize various aspects of the problem and solution parameters.
Here is a breakdown of the three main parameter files:

1. parametersProblem.xml: This file contains all the parameters related to the specific details of the definition of a problem. It includes physical constants, like e.g. density, or solver related information, like e.g. solver tolerances, as well as options for visualization export. 
                          By modifying this file, the user can adjust the problem-specific settings.

2. parametersSolver.xml:  This file includes parameters that can be used to specify the setting of the linear solver and nonlinear solver. Here, the linear solver type is often Belos, which indicates the utilization of the Belos package. 
                          If the user sets in the parametersProblem.xml file, the linearization to be NOX they can specify here NOX-related parameters, such as globalization strategies and related parameters.

3. parametersPrec.xml:    This file is dedicated to parameters related to the usage of preconditioners. If the parametersProblem.xml file specifies the Monolithic option, the preconditioner type should be set to Frosch, as it is implemented in the Frosch package. 
                          On the other hand, if the Blockpreconditioner Teko was specified, the corresponding Teko-related parameters should be set in this file.

We include here four different files that provide an overview of the meaning and options of different parameters inside these files. The  parametersPrec.xml File is here further divided into parametersPrecMonolithic.xml and parametersPrecTeko.xml files, as they have distinct keywords. 
It's important to note that this overview is not intended to be copied directly but rather to provide a comprehensive understanding of the parameter files and their purpose. If any new relevant keywords are introduced, they should be added to the list along with an explanation of their usage and significance.


-----------------------------------------------------------------------------

Collection of detailed explanations
===================================


-----------------------------------
RGDSW coarse space -- Coarse space options

RGDSW: reduced-dimension GDSW
GDSW: generalized Dryja-Smith-Widlund

The implemented coarse spaces of the RGDSWCoarseOperator are based on the following paper:

(*) Clark R. Dohrmann, Olof B. Widlund, 2017, On the Design of Small Coarse Spaces for Domain Decomposition Algorithms
    https://doi.org/10.1137/17M1114272

Dohrmann/Widlund have proposed several variants of energy-minimizing coarse spaces. 
(Technically, they are only energy-minimizing for symmetric, positive definite problems.) 
These prescribe a partition of unity on the interface, which is subsequently extended to the interior of the subdomains.
Not all variants are implemented in FROSch (which FEDDLib uses). 
Furthermore, the variant coming from adaptive coarse spaces is also not implemented; cf. Fig. 1 in
    A. Heinlein, A. Klawonn, J. Knepper, O. Rheinbach, O. B. Widlund, 2022, Adaptive GDSW Coarse Spaces of Reduced Dimension for Overlapping Schwarz Methods
    https://doi.org/10.1137/20M1364540

The variants of the RGDSW coarse spaces in (*) are named Option 1, Option 2.1, Option 2.2, and Option 2.1+2. 
Thereof, Option 1 and Option 2.2 are implemented.

Option 1: Multiplicity scaling on the interface
   Prescribe the inverse of the "coarse node multiplicity" to a node/dof on the interface.
   Eq. (1) in (*)
   Coarse node multiplicity: Let n be a node/dof. How many coarse nodes are ancestors of n?
      Ancestor: The set of adjacent subdomains of a node/dof is a subset of (not equal to) the set of adjacent subdomains of the ancestor.
      Detailed description: RGDSW generates interface components that overlap.
      Each interface component is (generally) associated with a coarse node. 
      For each coarse node, adjacent coarse edges and coarse faces are added to obtain a "coarse star".
      This coarse star is the RGDSW interface component and overlaps with neighboring coarse stars.
      The question now is: To which coarse stars does a node belong? This is the sought multiplicity.

Option 2.2: Distance scaling on the interface
   Prescribe sth. like (1/d1(n)) / (1/d1(n) + ... + 1/dk(n)) to a node/dof n.
   di(n) is the euclidian distance of n to its ancestor coarse node i in {1,...,k}.
   Eq. (5) in (*)


Example for choosing the coarse space variants in the parameter file:

Current operator implementation (via interface partition of unity coarse operator):

<Parameter name="CoarseOperator Type" type="string" value="IPOUHarmonicCoarseOperator"/>
<ParameterList name="IPOUHarmonicCoarseOperator">
    <ParameterList name="Blocks">
        <ParameterList name="1">
            <ParameterList name="InterfacePartitionOfUnity">
                <Parameter name="Type" type="string" value="RGDSW" />
                <ParameterList name="RGDSW">
                    <Parameter name="Distance Function" type="string" value="Constant" />           <!-- Option 1
                    or
                    <Parameter name="Distance Function" type="string" value="Inverse Euclidean" />  <!-- Option 2.2

Legacy operator implementation (via RGDSW coarse operator or GDSW coarse operator):

<Parameter name="CoarseOperator Type" type="string" value="RGDSWCoarseOperator"/>
<ParameterList name="RGDSWCoarseOperator">
    <ParameterList name="Blocks">
        <ParameterList name="1">
            <Parameter name="Option" type="string" value="1"/>
            or 
            <Parameter name="Option" type="string" value="2.2"/>

-----------------------------------
Mesh file format

* The FEDDLib can read MESH (INRIA) files (file ending .mesh). 
  - Binary files are not supported.
  - Not all features of the MESH format are supported.
  
  Overview of the format:
     [header]
	 Vertices
	 [number of nodes]
     [node coordinates, flags]
	 Edges
	 [number of edges]
	 [list of edges (n x 2) matrix, flags]
	 Triangles
	 [number of triangles]
	 [list of triangles (m x 3) matrix, flags]
	 Tetrahedra
	 [number of tetrahedra]
	 [list of tetrahedra (s x 4) matrix, flags]
  
  - The FEDDLib also ignores some of the header that is required in INRIA files. 
    This is actually an advantage, since it makes the reader more robust.
    For example: It does not matter if the coordinates are written with single or double precision.
  - During import, a mesh file is read multiple times. 
    This is, of course, not strictly necessary but simplified the initial implementation of the importer. 
    For the examples considered so far, the additional time required is negligible.
* The FEDDLib also has an importer for a specific format of Gmsh files (MSH file extension). 
  It was not written by the developers of the FEDDLib, and its code base was last updated in 2014. 
  The importer is not currently used and was also not used extensively in the past. 
  1) MESH format is simple (easy to maintain), works, and does all we need.
  2) Gmsh format is more complicated than the MESH format. An importer is more difficult and costly to maintain.
  3) Both formats are not widely supported.
     But the MESH format is more geared towards simple meshes, while the Gmsh format is more geared towards geometries (e.g., supporting parametrization), i.e., features that we do not require.
* VTK support in FEDDLib would be good, since VTK is more widely supported (e.g., (partially) by Gmsh and ParaView). 
  It could also simplify the process of setting flags, since these can be implemented in VTK as fields (i.e., scalar function defined on the mesh nodes). 
  Unfortunately, the specifications for the latest VTK format have not been published (Gmsh only supports the older version). 
  Ideally, there exists a C++ importer for VTK from the developers of VTK. 
  However, the VTK format is not overly complicated if only a small subset of features is supported. 
  A self-written importer might, thus, be feasible.
  VTK could also be used as an alternative for exporting variables from the FEDDLib (instead of HDF5). (It is much easier to deal with an ASCII VTK file than with an HDF5 file.)


How to prepare a mesh / prerequisites that need to be satisfied:

* Mesh file must be in INRIA-MESH format.
* Mesh must only contain P1 element nodes (i.e., a triangle has three nodes; a P2 triangle would have 6 nodes).
* Fluid-structure interaction: fluid mesh and structure mesh must (currently) be separate files.
* Each node is associated with a node flag. 
  The flag indicates what type of boundary condition is set (on the surface). 
  For interior nodes, a default flag is used (see below).
* Providing only flagged nodes is not sufficient to set boundary conditions (e.g., if a P2 mesh is generated, it may be undefined at the interface of two boundary conditions what condition should be prescribed). 
  Thus, elements that make up the interface between two boundary conditions need to be specified in the mesh file as well, along with the appropriate flag.
  
  Example:
     Consider a tetrahedralized tube with an inflow and outflow. 
     We set a velocity profile at the inflow, a zero flow at the wall, and a do-nothing boundary condition at the outflow. 
     So we need to differentiate between three types of triangles on the surface of the tube. 
     Each triangle receives a specific flag s.t. the FEDDLib knows which boundary condition to set. 
     However, it is unclear what condition should be set at the boundary nodes of the outlet, where Dirichlet and Neumann condition both meet. 
     Should it be a no-slip condition (yes!) or a do-nothing boundary condition (no! In this case fluid would leak out to the side of the tube). 
     We need to specify all edges that make up the boundary of the outlet and give them the same flag as the wall (i.e., a no-slip boundary condition).
  
  If no additional information (to avoid ambigious boundary conditions) was provided by the user (via the mesh), the FEDDLib always uses the adjacent surface node with the larger flag as reference. 
  This trickery is supported, but it is not recommended to rely on it since the implementation may change. 
  Moreover, it is less clear and, thus, error prone. 
  Only make use of it if the condition is unambigious. 
  In the example above, we do not need to provide the edges at the inlet, since the prescribed velocity function should be zero at the boundary of the inlet, as is the no-slip boundary condition of the wall. 
  Thus, in this case, both conditions are identical at their intersection.

-----------------------------------
