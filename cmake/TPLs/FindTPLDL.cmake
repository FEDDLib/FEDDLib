# We only want to search for the shared library.
# Enable the search for shared libraries and reset it to the user setting afterwards.
set(_SAVE_TPL_FIND_SHARED_LIBS "${TPL_FIND_SHARED_LIBS}")
set(TPL_FIND_SHARED_LIBS ON)


TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( DL
  REQUIRED_LIBS_NAMES "dl"
)

SET(CMAKE_REQUIRED_INCLUDES ${TPL_DL_INCLUDE_DIRS})
SET(CMAKE_REQUIRED_LIBRARIES ${TPL_DL_LIBRARY_DIRS})

# Reset user setting for finding shared libraries.
set(TPL_FIND_SHARED_LIBS "${_SAVE_TPL_FIND_SHARED_LIBS}")
unset(_SAVE_TPL_FIND_SHARED_LIBS)

message (" dl library found! ")
