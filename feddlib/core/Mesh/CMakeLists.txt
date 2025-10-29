set(Mesh_HEADERS
    Mesh/gmsh_io.hpp
	Mesh/Mesh.hpp
    Mesh/Mesh_decl.hpp
	Mesh/Mesh_def.hpp
	Mesh/MeshFileReader.hpp
	Mesh/MeshInterface.hpp
    Mesh/MeshInterface_decl.hpp
    Mesh/MeshInterface_def.hpp
    Mesh/MeshPartitioner.hpp
    Mesh/MeshPartitioner_decl.hpp
    Mesh/MeshPartitioner_def.hpp
    Mesh/MeshStructured.hpp
    Mesh/MeshStructured_decl.hpp
    Mesh/MeshStructured_def.hpp
	Mesh/MeshUnstructured.hpp
	Mesh/MeshUnstructured_decl.hpp
	Mesh/MeshUnstructured_def.hpp
    Mesh/AABBTree.hpp
    Mesh/AABBTree_decl.hpp
    Mesh/AABBTree_def.hpp
CACHE INTERNAL "")

set(Mesh_SOURCES
    Mesh/gmsh_io.cpp
    Mesh/Mesh.cpp
    Mesh/MeshPartitioner.cpp
    Mesh/MeshFileReader.cpp
    Mesh/MeshInterface.cpp
    Mesh/MeshStructured.cpp
    Mesh/MeshUnstructured.cpp
    Mesh/AABBTree.cpp
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${Mesh_HEADERS})
  string(REPLACE "Mesh/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
