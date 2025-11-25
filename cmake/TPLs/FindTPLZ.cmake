# We only want to search for the shared library.
# Enable the search for shared libraries and reset it to the user setting afterwards.
set(_SAVE_TPL_FIND_SHARED_LIBS "${TPL_FIND_SHARED_LIBS}")
set(TPL_FIND_SHARED_LIBS ON)

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Z
  REQUIRED_LIBS_NAMES "z"
)

SET(CMAKE_REQUIRED_INCLUDES ${TPL_Z_INCLUDE_DIRS})
SET(CMAKE_REQUIRED_LIBRARIES ${TPL_Z_LIBRARY_DIRS})

# Reset user setting for finding shared libraries.
set(TPL_FIND_SHARED_LIBS "${_SAVE_TPL_FIND_SHARED_LIBS}")
unset(_SAVE_TPL_FIND_SHARED_LIBS)

message (" z library found! ")
