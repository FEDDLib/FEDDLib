tribits_add_executable(
	laplaceUnitTest
	SOURCES laplaceUnitTest.cpp
)

tribits_add_test(
	laplaceUnitTest
	NAME laplace_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	laplaceUnitTest
	NAME laplace_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	laplaceUnitTest
	NAME laplace_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
tribits_add_test(
	laplaceUnitTest
	NAME laplace_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)


tribits_copy_files_to_binary_dir(data_laplaceUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Laplace.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Laplace.xml  parametersSolver.xml
	EXEDEPS laplaceUnitTest
)

tribits_copy_files_to_binary_dir(data_laplacetests
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	EXEDEPS laplaceUnitTest
)
# -----------
# Linear Elasticity
# -----------
tribits_add_executable(
	linElasUnitTest
	SOURCES linElasUnitTest.cpp
)

tribits_add_test(
	linElasUnitTest
	NAME linElas_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	linElasUnitTest
	NAME linElas_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	linElasUnitTest
	NAME linElas_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
tribits_add_test(
	linElasUnitTest
	NAME linElas_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)

tribits_copy_files_to_binary_dir(data_linElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_linElas_2d_P1_4cores.h5 solution_linElas_2d_P2_4cores.h5 solution_linElas_3d_P1_8cores.h5 solution_linElas_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_linElas_2d_P1_4cores.h5 solution_linElas_2d_P2_4cores.h5 solution_linElas_3d_P1_8cores.h5 solution_linElas_3d_P2_8cores.h5
	EXEDEPS linElasUnitTest
)

# -----------
# Stokes
# -----------
tribits_add_executable(
	stokesUnitTest
	SOURCES stokesUnitTest.cpp
)
tribits_add_test(
	stokesUnitTest
	NAME stokes_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	stokesUnitTest
	NAME stokes_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	stokesUnitTest
	NAME stokes_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
tribits_add_test(
	stokesUnitTest
	NAME stokes_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)
tribits_copy_files_to_binary_dir(data_stokesUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Stokes_2D.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Stokes_2D.xml  parametersSolver.xml
	EXEDEPS stokesUnitTest
)
tribits_copy_files_to_binary_dir(data_stokesTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_stokes_velocity_2d_P1_4cores.h5 solution_stokes_velocity_2d_P2_4cores.h5 solution_stokes_velocity_3d_P1_8cores.h5 solution_stokes_velocity_3d_P2_8cores.h5 solution_stokes_pressure_2d_P1_4cores.h5 solution_stokes_pressure_2d_P2_4cores.h5 solution_stokes_pressure_3d_P1_8cores.h5 solution_stokes_pressure_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_stokes_velocity_2d_P1_4cores.h5 solution_stokes_velocity_2d_P2_4cores.h5 solution_stokes_velocity_3d_P1_8cores.h5 solution_stokes_velocity_3d_P2_8cores.h5 solution_stokes_pressure_2d_P1_4cores.h5 solution_stokes_pressure_2d_P2_4cores.h5 solution_stokes_pressure_3d_P1_8cores.h5 solution_stokes_pressure_3d_P2_8cores.h5
	EXEDEPS stokesUnitTest
)



# -----------
# Generalized Newtonian Power Law
# -----------
tribits_add_executable(
	generalizedNewtonianPowerLawTest
	SOURCES generalizedNewtonianPowerLawTest.cpp
)

tribits_add_test(
	generalizedNewtonianPowerLawTest
	NAME GNF_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_GNF_2D.xml parametersSolver.xml parametersProblem_GNF_2D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_GNF_2D.xml  parametersSolver.xml parametersProblem_GNF_2D.xml
	EXEDEPS generalizedNewtonianPowerLawTest
)

tribits_copy_files_to_binary_dir(data2_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_GNF_pressure_2d_P2P1_4cores.h5 solution_GNF_velocity_2d_P2P1_4cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_GNF_pressure_2d_P2P1_4cores.h5 solution_GNF_velocity_2d_P2P1_4cores.h5
	EXEDEPS generalizedNewtonianPowerLawTest
)

tribits_copy_files_to_binary_dir(mesh_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Meshes
	SOURCE_FILES rectangle_200.mesh rectangle_100.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/Meshes
	DEST_FILES rectangle_200.mesh rectangle_100.mesh
	EXEDEPS generalizedNewtonianPowerLawTest
)

# -----------
# Nonlinear Elasticity
# -----------
tribits_add_executable(
  nonLinElasticityUnitTest
	SOURCES nonLinElasticityUnitTest.cpp
)

tribits_add_test(
	nonLinElasticityUnitTest
	NAME nonLinElasticity_2D_P1
	ARGS "--dim=2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
  nonLinElasticityUnitTest
	NAME nonLinElasticity_2D_P1
	ARGS "--dim=3"
	NUM_MPI_PROCS 8
)

tribits_copy_files_to_binary_dir(param_nonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES parametersPrecNonLinElasticity.xml parametersSolver.xml parametersProblemNonLinElasticity.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  DEST_FILES parametersPrecNonLinElasticity.xml  parametersSolver.xml parametersProblemNonLinElasticity.xml
  EXEDEPS nonLinElasticityUnitTest
)

tribits_copy_files_to_binary_dir(data_nonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
  SOURCE_FILES solution_nonLinElasticity_2d_P1_4cores.h5 solution_nonLinElasticity_3d_P1_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
  DEST_FILES solution_nonLinElasticity_2d_P1_4cores.h5 solution_nonLinElasticity_3d_P1_8cores.h5
	EXEDEPS nonLinElasticityUnitTest
)


# -----------
# Nonlinear Laplace
# -----------
tribits_add_executable(
  nonLinLaplaceUnitTest
  SOURCES nonLinLaplaceUnitTest.cpp
)

tribits_add_test(
  nonLinLaplaceUnitTest
  NAME nonLinLaplace_2D_P1
	ARGS "--dim=2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
  nonLinLaplaceUnitTest
  NAME nonLinLaplace_2D_P1
	ARGS "--dim=3"
	NUM_MPI_PROCS 8
)

tribits_copy_files_to_binary_dir(param_nonLinLaplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES parametersPrecNonLinLaplace.xml parametersSolver.xml parametersProblemNonLinLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  DEST_FILES parametersPrecNonLinLaplace.xml  parametersSolver.xml parametersProblemNonLinLaplace.xml
	EXEDEPS nonLinLaplaceUnitTest
)

tribits_copy_files_to_binary_dir(data_nonLinLaplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
  SOURCE_FILES solution_nonLinLaplace_2d_P1_4cores.h5 solution_nonLinLaplace_3d_P1_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
  DEST_FILES solution_nonLinLaplace_2d_P1_4cores.h5 solution_nonLinLaplace_3d_P1_8cores.h5
  EXEDEPS nonLinLaplaceUnitTest
)



# -----------
# Navier-Stokes
# -----------
tribits_add_executable(
	navierStokesUnitTest
	SOURCES navierStokesUnitTest.cpp
)
# tribits_add_test(
# 	navierStokesUnitTest
# 	NAME navierStokes_2D_P1_NOX
# 	ARGS "--dim=2 --FEType=P1 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 4
# )
# tribits_add_test(
# 	navierStokesUnitTest
# 	NAME navierStokes_2D_P2_NOX
# 	ARGS "--dim=2 --FEType=P2 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 4
# )
# tribits_add_test(
# 	navierStokesUnitTest
# 	NAME navierStokes_3D_P1_NOX
# 	ARGS "--dim=3 --FEType=P1 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 8
# )
# tribits_add_test(
# 	navierStokesUnitTest
# 	NAME navierStokes_3D_P2_NOX
# 	ARGS "--dim=3 --FEType=P2 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 8
# )

tribits_add_test(
	navierStokesUnitTest
	NAME navierStokes_2D_P1_Newton
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	navierStokesUnitTest
	NAME navierStokes_2D_P2_Newton
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
tribits_add_test(
	navierStokesUnitTest
	NAME navierStokes_3D_P1_Newton
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
tribits_add_test(
	navierStokesUnitTest
	NAME navierStokes_3D_P2_Newton
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)


tribits_copy_files_to_binary_dir(param_navierStokesUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Stokes_2D.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Stokes_2D.xml  parametersSolver.xml
	EXEDEPS navierStokesUnitTest
)
tribits_copy_files_to_binary_dir(data_navierStokesUnitTest
 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
 	SOURCE_FILES solution_navier_stokes_velocity_2d_P1_4cores.h5 solution_navier_stokes_velocity_2d_P2_4cores.h5 solution_navier_stokes_velocity_3d_P1_8cores.h5 solution_navier_stokes_velocity_3d_P2_8cores.h5 solution_navier_stokes_pressure_2d_P1_4cores.h5 solution_navier_stokes_pressure_2d_P2_4cores.h5 solution_navier_stokes_pressure_3d_P1_8cores.h5 solution_navier_stokes_pressure_3d_P2_8cores.h5
 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
 	DEST_FILES solution_navier_stokes_velocity_2d_P1_4cores.h5 solution_navier_stokes_velocity_2d_P2_4cores.h5 solution_navier_stokes_velocity_3d_P1_8cores.h5 solution_navier_stokes_velocity_3d_P2_8cores.h5 solution_navier_stokes_pressure_2d_P1_4cores.h5 solution_navier_stokes_pressure_2d_P2_4cores.h5 solution_navier_stokes_pressure_3d_P1_8cores.h5 solution_navier_stokes_pressure_3d_P2_8cores.h5
 	EXEDEPS navierStokesUnitTest
)

# -----------
# Navier-Stokes
# -----------
tribits_add_executable(
	navierStokesUnitTestPrec
	SOURCES navierStokesUnitTestPrec.cpp
)

tribits_copy_files_to_binary_dir(param_navierStokesUnitTestPrec
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_NavierStokes_Block.xml parametersPrec_NavierStokes_Mono.xml parametersPrec_NavierStokes_Teko.xml parametersSolver_PrecTest.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_NavierStokes_Block.xml parametersPrec_NavierStokes_Mono.xml parametersPrec_NavierStokes_Teko.xml parametersSolver_PrecTest.xml
	EXEDEPS navierStokesUnitTestPrec
)

tribits_add_test(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P1_Monolithic
	ARGS "--dim=2 --FEType=P1 --precMethod=Monolithic"
	NUM_MPI_PROCS 9
)
tribits_add_test(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_Block
	ARGS "--dim=2 --FEType=P1 --precMethod=Diagonal"
	NUM_MPI_PROCS 9
)
tribits_add_test(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_PCD
	ARGS "--dim=2 --FEType=P1 --precMethod=PCD"
	NUM_MPI_PROCS 9
)
tribits_add_test(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_Teko
	ARGS "--dim=2 --FEType=P1 --precMethod=Teko"
	NUM_MPI_PROCS 9
)

tribits_copy_files_to_binary_dir(data_navierStokesUnitTestPrec
 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
 	SOURCE_FILES solution_navier_stokes_velocity_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_velocity_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_velocity_Teko_2d_P1_9cores.h5 solution_navier_stokes_velocity_PCD_2d_P1_9cores.h5 solution_navier_stokes_pressure_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_pressure_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_pressure_Teko_2d_P1_9cores.h5 solution_navier_stokes_pressure_PCD_2d_P1_9cores.h5  
 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
 	DEST_FILES solution_navier_stokes_velocity_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_velocity_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_velocity_Teko_2d_P1_9cores.h5 solution_navier_stokes_velocity_PCD_2d_P1_9cores.h5 solution_navier_stokes_pressure_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_pressure_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_pressure_Teko_2d_P1_9cores.h5 solution_navier_stokes_pressure_PCD_2d_P1_9cores.h5  
 	EXEDEPS navierStokesUnitTestPrec
)