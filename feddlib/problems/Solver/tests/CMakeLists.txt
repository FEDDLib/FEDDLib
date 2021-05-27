INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### Domain
TRIBITS_ADD_EXECUTABLE(
facsi
SOURCES facsi.cpp
)

TRIBITS_ADD_TEST(
facsi
NAME facsi
ARGS
NUM_MPI_PROCS 2
)
