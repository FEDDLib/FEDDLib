TRIBITS_ADD_EXECUTABLE(
	laplaceUnitTest
	SOURCES laplaceUnitTest.cpp
)

TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Laplace.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Laplace.xml  parametersSolver.xml
	EXEDEPS laplaceUnitTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplacetests
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	EXEDEPS laplaceUnitTest
)
# -----------
# Linear Elasticity
# -----------
TRIBITS_ADD_EXECUTABLE(
	linElasUnitTest
	SOURCES linElasUnitTest.cpp
)

TRIBITS_ADD_TEST(
	linElasUnitTest
	NAME linElas_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	linElasUnitTest
	NAME linElas_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	linElasUnitTest
	NAME linElas_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	linElasUnitTest
	NAME linElas_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_linElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_linElas_2d_P1_4cores.h5 solution_linElas_2d_P2_4cores.h5 solution_linElas_3d_P1_8cores.h5 solution_linElas_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_linElas_2d_P1_4cores.h5 solution_linElas_2d_P2_4cores.h5 solution_linElas_3d_P1_8cores.h5 solution_linElas_3d_P2_8cores.h5
	EXEDEPS linElasUnitTest
)

# -----------
# Stokes
# -----------
TRIBITS_ADD_EXECUTABLE(
	stokesUnitTest
	SOURCES stokesUnitTest.cpp
)
TRIBITS_ADD_TEST(
	stokesUnitTest
	NAME stokes_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	stokesUnitTest
	NAME stokes_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	stokesUnitTest
	NAME stokes_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	stokesUnitTest
	NAME stokes_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_stokesUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Stokes_2D.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Stokes_2D.xml  parametersSolver.xml
	EXEDEPS stokesUnitTest
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_stokesTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_stokes_velocity_2d_P1_4cores.h5 solution_stokes_velocity_2d_P2_4cores.h5 solution_stokes_velocity_3d_P1_8cores.h5 solution_stokes_velocity_3d_P2_8cores.h5 solution_stokes_pressure_2d_P1_4cores.h5 solution_stokes_pressure_2d_P2_4cores.h5 solution_stokes_pressure_3d_P1_8cores.h5 solution_stokes_pressure_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_stokes_velocity_2d_P1_4cores.h5 solution_stokes_velocity_2d_P2_4cores.h5 solution_stokes_velocity_3d_P1_8cores.h5 solution_stokes_velocity_3d_P2_8cores.h5 solution_stokes_pressure_2d_P1_4cores.h5 solution_stokes_pressure_2d_P2_4cores.h5 solution_stokes_pressure_3d_P1_8cores.h5 solution_stokes_pressure_3d_P2_8cores.h5
	EXEDEPS stokesUnitTest
)



# -----------
# Generalized Newtonian Power Law
# -----------
TRIBITS_ADD_EXECUTABLE(
	generalizedNewtonianPowerLawTest
	SOURCES generalizedNewtonianPowerLawTest.cpp
)

TRIBITS_ADD_TEST(
	generalizedNewtonianPowerLawTest
	NAME GNF_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_GNF_2D.xml parametersSolver.xml parametersProblem_GNF_2D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_GNF_2D.xml  parametersSolver.xml parametersProblem_GNF_2D.xml
	EXEDEPS generalizedNewtonianPowerLawTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data2_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_GNF_pressure_2d_P2P1_4cores.h5 solution_GNF_velocity_2d_P2P1_4cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_GNF_pressure_2d_P2P1_4cores.h5 solution_GNF_velocity_2d_P2P1_4cores.h5
	EXEDEPS generalizedNewtonianPowerLawTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_generalizedNewtonianPowerLawTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Meshes
	SOURCE_FILES rectangle_200.mesh rectangle_100.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/Meshes
	DEST_FILES rectangle_200.mesh rectangle_100.mesh
	EXEDEPS generalizedNewtonianPowerLawTest
)

# -----------
# Nonlinear Elasticity
# -----------
TRIBITS_ADD_EXECUTABLE(
  nonLinElasticityUnitTest
	SOURCES nonLinElasticityUnitTest.cpp
)

TRIBITS_ADD_TEST(
	nonLinElasticityUnitTest
	NAME nonLinElasticity_2D_P1
	ARGS "--dim=2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
  nonLinElasticityUnitTest
	NAME nonLinElasticity_2D_P1
	ARGS "--dim=3"
	NUM_MPI_PROCS 8
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(param_nonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES parametersPrecNonLinElasticity.xml parametersSolver.xml parametersProblemNonLinElasticity.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  DEST_FILES parametersPrecNonLinElasticity.xml  parametersSolver.xml parametersProblemNonLinElasticity.xml
  EXEDEPS nonLinElasticityUnitTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_nonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
  SOURCE_FILES solution_nonLinElasticity_2d_P1_4cores.h5 solution_nonLinElasticity_3d_P1_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
  DEST_FILES solution_nonLinElasticity_2d_P1_4cores.h5 solution_nonLinElasticity_3d_P1_8cores.h5
	EXEDEPS nonLinElasticityUnitTest
)


# -----------
# Nonlinear Laplace
# -----------
TRIBITS_ADD_EXECUTABLE(
  nonLinLaplaceUnitTest
  SOURCES nonLinLaplaceUnitTest.cpp
)

TRIBITS_ADD_TEST(
  nonLinLaplaceUnitTest
  NAME nonLinLaplace_2D_P1
	ARGS "--dim=2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
  nonLinLaplaceUnitTest
  NAME nonLinLaplace_2D_P1
	ARGS "--dim=3"
	NUM_MPI_PROCS 8
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(param_nonLinLaplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES parametersPrecNonLinLaplace.xml parametersSolver.xml parametersProblemNonLinLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  DEST_FILES parametersPrecNonLinLaplace.xml  parametersSolver.xml parametersProblemNonLinLaplace.xml
	EXEDEPS nonLinLaplaceUnitTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_nonLinLaplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
  SOURCE_FILES solution_nonLinLaplace_2d_P1_4cores.h5 solution_nonLinLaplace_3d_P1_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
  DEST_FILES solution_nonLinLaplace_2d_P1_4cores.h5 solution_nonLinLaplace_3d_P1_8cores.h5
  EXEDEPS nonLinLaplaceUnitTest
)



# -----------
# Navier-Stokes
# -----------
TRIBITS_ADD_EXECUTABLE(
	navierStokesUnitTest
	SOURCES navierStokesUnitTest.cpp
)
# TRIBITS_ADD_TEST(
# 	navierStokesUnitTest
# 	NAME navierStokes_2D_P1_NOX
# 	ARGS "--dim=2 --FEType=P1 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 4
# )
# TRIBITS_ADD_TEST(
# 	navierStokesUnitTest
# 	NAME navierStokes_2D_P2_NOX
# 	ARGS "--dim=2 --FEType=P2 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 4
# )
# TRIBITS_ADD_TEST(
# 	navierStokesUnitTest
# 	NAME navierStokes_3D_P1_NOX
# 	ARGS "--dim=3 --FEType=P1 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 8
# )
# TRIBITS_ADD_TEST(
# 	navierStokesUnitTest
# 	NAME navierStokes_3D_P2_NOX
# 	ARGS "--dim=3 --FEType=P2 --nonlinearSolver=NOX"
# 	NUM_MPI_PROCS 8
# )

TRIBITS_ADD_TEST(
	navierStokesUnitTest
	NAME navierStokes_2D_P1_Newton
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	navierStokesUnitTest
	NAME navierStokes_2D_P2_Newton
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	navierStokesUnitTest
	NAME navierStokes_3D_P1_Newton
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	navierStokesUnitTest
	NAME navierStokes_3D_P2_Newton
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(param_navierStokesUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Stokes_2D.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Stokes_2D.xml  parametersSolver.xml
	EXEDEPS navierStokesUnitTest
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_navierStokesUnitTest
 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
 	SOURCE_FILES solution_navier_stokes_velocity_2d_P1_4cores.h5 solution_navier_stokes_velocity_2d_P2_4cores.h5 solution_navier_stokes_velocity_3d_P1_8cores.h5 solution_navier_stokes_velocity_3d_P2_8cores.h5 solution_navier_stokes_pressure_2d_P1_4cores.h5 solution_navier_stokes_pressure_2d_P2_4cores.h5 solution_navier_stokes_pressure_3d_P1_8cores.h5 solution_navier_stokes_pressure_3d_P2_8cores.h5
 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
 	DEST_FILES solution_navier_stokes_velocity_2d_P1_4cores.h5 solution_navier_stokes_velocity_2d_P2_4cores.h5 solution_navier_stokes_velocity_3d_P1_8cores.h5 solution_navier_stokes_velocity_3d_P2_8cores.h5 solution_navier_stokes_pressure_2d_P1_4cores.h5 solution_navier_stokes_pressure_2d_P2_4cores.h5 solution_navier_stokes_pressure_3d_P1_8cores.h5 solution_navier_stokes_pressure_3d_P2_8cores.h5
 	EXEDEPS navierStokesUnitTest
)

# -----------
# Navier-Stokes
# -----------
TRIBITS_ADD_EXECUTABLE(
	navierStokesUnitTestPrec
	SOURCES navierStokesUnitTestPrec.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(param_navierStokesUnitTestPrec
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_NavierStokes_Block.xml parametersPrec_NavierStokes_Mono.xml parametersPrec_NavierStokes_Teko.xml parametersSolver_PrecTest.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_NavierStokes_Block.xml parametersPrec_NavierStokes_Mono.xml parametersPrec_NavierStokes_Teko.xml parametersSolver_PrecTest.xml
	EXEDEPS navierStokesUnitTestPrec
)

TRIBITS_ADD_TEST(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P1_Monolithic
	ARGS "--dim=2 --FEType=P1 --precMethod=Monolithic"
	NUM_MPI_PROCS 9
)
TRIBITS_ADD_TEST(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_Block
	ARGS "--dim=2 --FEType=P1 --precMethod=Diagonal"
	NUM_MPI_PROCS 9
)
TRIBITS_ADD_TEST(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_PCD
	ARGS "--dim=2 --FEType=P1 --precMethod=PCD"
	NUM_MPI_PROCS 9
)
TRIBITS_ADD_TEST(
	navierStokesUnitTestPrec
	NAME navierStokes_2D_P2_Teko
	ARGS "--dim=2 --FEType=P1 --precMethod=Teko"
	NUM_MPI_PROCS 9
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_navierStokesUnitTestPrec
 	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
 	SOURCE_FILES solution_navier_stokes_velocity_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_velocity_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_velocity_Teko_2d_P1_9cores.h5 solution_navier_stokes_velocity_PCD_2d_P1_9cores.h5 solution_navier_stokes_pressure_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_pressure_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_pressure_Teko_2d_P1_9cores.h5 solution_navier_stokes_pressure_PCD_2d_P1_9cores.h5  
 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
 	DEST_FILES solution_navier_stokes_velocity_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_velocity_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_velocity_Teko_2d_P1_9cores.h5 solution_navier_stokes_velocity_PCD_2d_P1_9cores.h5 solution_navier_stokes_pressure_Diagonal_2d_P1_9cores.h5 solution_navier_stokes_pressure_Monolithic_2d_P1_9cores.h5 solution_navier_stokes_pressure_Teko_2d_P1_9cores.h5 solution_navier_stokes_pressure_PCD_2d_P1_9cores.h5  
 	EXEDEPS navierStokesUnitTestPrec
)