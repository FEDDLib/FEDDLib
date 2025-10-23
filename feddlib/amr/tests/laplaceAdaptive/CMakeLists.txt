tribits_add_executable(
	laplaceAdaptive
	SOURCES main.cpp
)

tribits_add_test(
	laplaceAdaptive
    NAME laplaceAdaptive2D
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

tribits_add_test(
	laplaceAdaptive
    NAME laplaceAdaptive3D
    ARGS 
	--problemfile=parametersProblem3D.xml
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_laplaceAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml 
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive2D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES lShape4.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape4.mesh
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive2D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES lShape1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape1.mesh
	EXEDEPS laplaceAdaptive
)


tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive2D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES square_paper.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square_paper.mesh
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive2D3
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES lShape3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape3D.mesh
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES cube_131.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube_131.mesh
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive3D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES cube_14.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube_14.mesh
	EXEDEPS laplaceAdaptive
)

tribits_copy_files_to_binary_dir(mesh_laplaceAdaptive3D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES cutCube.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cutCube.mesh
	EXEDEPS laplaceAdaptive
)

