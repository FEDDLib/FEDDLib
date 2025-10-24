### AABBTree
tribits_add_executable(
AABBTree
SOURCES AABBTree.cpp
)

tribits_add_test(
AABBTree
NAME AABBTree
ARGS
NUM_MPI_PROCS 1
)

tribits_copy_files_to_binary_dir(meshes_meshes_AABBTree
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES rectangle.mesh rectangle.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES rectangle.mesh rectangle.mesh
EXEDEPS AABBTree
)

### mesh_AABBTree
tribits_add_executable(
mesh_AABBTree
SOURCES mesh_AABBTree.cpp
)

tribits_add_test(
mesh_AABBTree
NAME mesh_AABBTree
ARGS
NUM_MPI_PROCS 1
)

tribits_copy_files_to_binary_dir(meshes_mesh_AABBTree
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES big_rectangle.mesh big_rectangle.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES big_rectangle.mesh big_rectangle.mesh
EXEDEPS mesh_AABBTree
)

### mesh_structured
tribits_add_executable(
mesh_structured
SOURCES mesh_structured.cpp
)

tribits_add_test(
mesh_structured
NAME mesh_structured
ARGS
NUM_MPI_PROCS 4
)

### mesh_structured_3D_5_element
tribits_add_executable(
mesh_structured_3D_5_element
SOURCES mesh_structured_5_element_subcube.cpp
)

tribits_add_test(
mesh_structured_3D_5_element
NAME mesh_structured_3D_5_element
ARGS
NUM_MPI_PROCS 8
)

### mesh_structured_3DCR
#tribits_add_executable(
#    mesh_structured_3DCR
#    SOURCES mesh_structured_3DCR.cpp
# )

#tribits_add_test(
#    mesh_structured_3DCR
#    NAME mesh_structured_3DCR
#    ARGS
#    NUM_MPI_PROCS 1
# )

### mesh_structured_3DCR_BFS
#tribits_add_executable(
#    mesh_structured_3DCR_BFS
#    SOURCES mesh_structured_3DCR_BFS.cpp
#)

#tribits_add_test(
#mesh_structured_3DCR_BFS
#    NAME mesh_structured_3DCR_BFS
#    ARGS
#    NUM_MPI_PROCS 3
#)

### mesh_elementFlags
tribits_add_executable(
mesh_elementFlags
SOURCES mesh_elementFlags.cpp
)

tribits_add_test(
mesh_elementFlags
NAME mesh_elementFlags
ARGS
NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_mesh_elementFlags
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES testFoam2.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES testFoam2.mesh
EXEDEPS mesh_elementFlags
)

### mesh_nodeFlagsP2
tribits_add_executable(
mesh_nodeFlagsP2
SOURCES mesh_nodeFlagsP2.cpp
)

tribits_add_test(
mesh_nodeFlagsP2
NAME mesh_nodeFlagsP2
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_mesh_nodeFlagsP2
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES meshNodeTestP23D_2.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES meshNodeTestP23D_2.mesh
EXEDEPS mesh_nodeFlagsP2
)

### testing mesh surfaces orientation
tribits_add_executable(
mesh_surfaceOrientation
SOURCES mesh_surfaceOrientation.cpp
)

tribits_add_test(
mesh_surfaceOrientation
NAME mesh_surfaceOrientation
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_mesh_surfaceOrientation
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES cube_tetr6_struct_h=0_2.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES cube_tetr6_struct_h=0_2.mesh
EXEDEPS mesh_surfaceOrientation
)

## mesh_unstructured_preProcessing
tribits_add_executable(
mesh_unstructured_preProcessing
SOURCES mesh_unstructured_preProcessing.cpp
)

tribits_add_test(
mesh_unstructured_preProcessing
NAME mesh_unstructured_preProcessing_3D_cube_h_1
ARGS "--FEType=P2 --dim=3 --file=cube_h_1.mesh"
NUM_MPI_PROCS 4
)

tribits_add_test(
mesh_unstructured_preProcessing
NAME mesh_unstructured_preProcessing_3D_fluid_mesh
ARGS "--FEType=P2 --dim=3 --file=fluidBenchmark2.mesh --ExportName=3D_fluid_mesh"
NUM_MPI_PROCS 4
)

tribits_add_test(
mesh_unstructured_preProcessing
NAME mesh_unstructured_preProcessing_2D_rectangle
ARGS "--FEType=P2 --dim=2 --file=rectangle.mesh --ExportName=2D_rectangle_mesh"
NUM_MPI_PROCS 4
)

tribits_add_test(
	mesh_unstructured_writeMesh
    NAME mesh_unstructured_writeMesh_3D_P2
    ARGS "--FEType=P2 --dim=3 --file=cube_h_1.mesh ---ExportName=2D_cube_mesh"
    NUM_MPI_PROCS 3
)

tribits_copy_files_to_binary_dir(data_mesh_unstructured_preProcessing
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES fluidBenchmark2.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES fluidBenchmark2.mesh
EXEDEPS mesh_unstructured_preProcessing
)

### mesh_unstructured
tribits_add_executable(
mesh_unstructured
SOURCES mesh_unstructured.cpp
)

tribits_add_test(
mesh_unstructured
NAME mesh_unstructured
ARGS
NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_mesh_unstructured
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES dfg_fsi_fluid_h004.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES dfg_fsi_fluid_h004.mesh
EXEDEPS mesh_unstructured
)

### Export of unstructured mesh
tribits_add_executable(
	mesh_unstructured_writeMesh
	SOURCES mesh_unstructured_writeMesh.cpp
)

tribits_add_test(
	mesh_unstructured_writeMesh
    NAME mesh_unstructured_writeMesh_2D_P2
    ARGS "--FEType=P2 --fileExport=export_2D_P2.mesh"
    NUM_MPI_PROCS 2
)

tribits_add_test(
	mesh_unstructured_writeMesh
    NAME mesh_unstructured_writeMesh_2D_P1
    ARGS "--fileExport=export_2D_P1.mesh"
    NUM_MPI_PROCS 4
)
tribits_add_test(
	mesh_unstructured_writeMesh
    NAME mesh_unstructured_writeMesh_3D_P1
    ARGS "--dim=3 --file=cube_h_1.mesh --fileExport=export_3D_P1.mesh"
    NUM_MPI_PROCS 5
)

tribits_copy_files_to_binary_dir(data_mesh_unstructured_writeMesh
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh cube_h_1.mesh cube_h_0_1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh cube_h_1.mesh cube_h_0_1.mesh
	EXEDEPS mesh_unstructured_writeMesh
)

### meshes_interface
tribits_add_executable(
meshes_interface
SOURCES meshes_interface.cpp
)

tribits_add_test(
meshes_interface
NAME meshes_interface
ARGS
NUM_MPI_PROCS 2
)

#tribits_add_test(
#    meshes_interface
#    NAME meshes_interface3d
#    ARGS --file=meshes_interface3D.xml
#    NUM_MPI_PROCS 2
#)

tribits_copy_files_to_binary_dir(data_meshes_interface
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
SOURCE_FILES meshes_interface.xml meshes_interface3D.xml
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES meshes_interface.xml meshes_interface3D.xml
EXEDEPS meshes_interface
)

tribits_copy_files_to_binary_dir(meshes_meshes_interface
SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
EXEDEPS meshes_interface
)

#tribits_copy_files_to_binary_dir(meshes_meshes_interface3d
#	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
#    SOURCE_FILES richter_3d_fluid_h0035.mesh richter_3d_solid_h0035.mesh
#    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#    DEST_FILES richter_3d_fluid_h0035.mesh richter_3d_solid_h0035.mesh
#    EXEDEPS meshes_interface
#)

### meshes_interface
tribits_add_executable(
meshes_custom_partition
SOURCES meshes_custom_partition.cpp
)

tribits_add_test(
meshes_custom_partition
NAME meshes_custom_partition
ARGS
NUM_MPI_PROCS 2
)

### Exporting mesh flags
tribits_add_executable(
mesh_flag_write
SOURCES mesh_flag_write.cpp
)

tribits_add_test(
mesh_flag_write
NAME mesh_flag_write
ARGS
NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_meshes_custom_partition
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
SOURCE_FILES meshes_custom_partition.xml
DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
 DEST_FILES meshes_custom_partition.xml
  EXEDEPS meshes_custom_partition
)

tribits_copy_files_to_binary_dir(meshes_meshes_custom_partition
 SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
  SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
   DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
  EXEDEPS meshes_custom_partition
)
