set(Utils_HEADERS
	Utils/FEDDUtils.hpp
CACHE INTERNAL "")

set(Utils_SOURCES
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${FE_HEADERS})
  string(REPLACE "Utils/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
