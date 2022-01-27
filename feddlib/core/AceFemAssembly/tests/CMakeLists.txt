INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### FE
TRIBITS_ADD_EXECUTABLE(
fe_test
SOURCES fe_test.cpp
)

TRIBITS_ADD_TEST(
fe_test
NAME fe_test
ARGS
NUM_MPI_PROCS 4
)


