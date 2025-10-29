### Domain
tribits_add_executable(
facsi
SOURCES facsi.cpp
)

tribits_add_test(
facsi
NAME facsi
ARGS
NUM_MPI_PROCS 2
)
