set(MPIEXEC_EXECTUABLE "mpirun")

tribits_add_executable(
	steadyNavierStokesAssFE
	SOURCES main.cpp
)
tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_2D_unstructured_P2_P1
	ARGS "--problemfile=parametersProblem.xml --precfile=parametersPrec.xml"
	COMM mpi
	NUM_MPI_PROCS 4
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_2D_structured_P1_P1
	ARGS "--problemfile=parametersProblemP1_structured.xml --precfile=parametersPrec.xml"
	COMM mpi
	NUM_MPI_PROCS 9
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_2D_structured_P2_P1
	ARGS "--problemfile=parametersProblemP2_structured.xml --precfile=parametersPrec.xml"
	NUM_MPI_PROCS 9
)
tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_3D_structured_P1_P1
	ARGS "--problemfile=parametersProblemP1_structured_3D.xml --precfile=parametersPrec_3D.xml"
	NUM_MPI_PROCS 9
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_3D_structured_P2_P1
	ARGS "--problemfile=parametersProblemP2_structured_3D.xml --precfile=parametersPrec_3D.xml"
	NUM_MPI_PROCS 9
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_2D_structured_P2_P1_Newton
	ARGS "--problemfile=parametersProblemP2_structured_Newton.xml --precfile=parametersPrec.xml"
	NUM_MPI_PROCS 9
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE_3D_structured_P2_P1_Newton
	ARGS "--problemfile=parametersProblemP2_structured_3D_Newton.xml --precfile=parametersPrec_3D.xml"
	NUM_MPI_PROCS 9
)

tribits_copy_files_to_binary_dir(data_steadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersPrec_3D.xml parametersProblem.xml parametersProblem_3D.xml parametersProblemP1_structured.xml parametersProblemP1_structured_3D.xml parametersProblemP2_structured_3D.xml parametersProblemP2_structured.xml parametersProblemP2_structured_Newton.xml parametersProblemP2_structured_3D_Newton.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersPrec_3D.xml parametersProblem.xml parametersProblem_3D.xml parametersProblemP1_structured.xml parametersProblemP1_structured_3D.xml parametersProblemP2_structured_3D.xml parametersProblemP2_structured.xml parametersProblemP2_structured_Newton.xml parametersProblemP2_structured_3D_Newton.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyNavierStokesAssFE
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h008.mesh DFG3DCylinder_1k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h008.mesh DFG3DCylinder_1k.mesh 
	EXEDEPS steadyNavierStokesAssFE
)




