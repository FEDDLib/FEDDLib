set(IO_HEADERS
	IO/IO_ParameterList_decl.hpp
CACHE INTERNAL "")

set(LOCAL_HEADERS)
foreach(INC ${IO_HEADERS})
  string(REPLACE "IO/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
