SET(LinearAlgebra_HEADERS
LinearAlgebra/BlockMap_decl.hpp
LinearAlgebra/BlockMap_def.hpp
LinearAlgebra/BlockMap.hpp
LinearAlgebra/BlockMatrix.hpp
LinearAlgebra/BlockMatrix_decl.hpp
LinearAlgebra/BlockMatrix_def.hpp
LinearAlgebra/BlockMultiVector.hpp
LinearAlgebra/BlockMultiVector_decl.hpp
LinearAlgebra/BlockMultiVector_def.hpp
LinearAlgebra/Map_decl.hpp
LinearAlgebra/Map_def.hpp
LinearAlgebra/Map.hpp
LinearAlgebra/Map_Xpetra_decl.hpp
LinearAlgebra/Map_Xpetra_def.hpp
LinearAlgebra/Map_Xpetra.hpp
LinearAlgebra/Matrix_decl.hpp
LinearAlgebra/Matrix_def.hpp
LinearAlgebra/Matrix.hpp
LinearAlgebra/MultiVector_decl.hpp
LinearAlgebra/MultiVector_def.hpp
LinearAlgebra/MultiVector.hpp
CACHE INTERNAL "")

SET(LinearAlgebra_SOURCES
LinearAlgebra/BlockMap.cpp
LinearAlgebra/BlockMatrix.cpp
LinearAlgebra/BlockMultiVector.cpp
LinearAlgebra/Map.cpp
LinearAlgebra/Map_Xpetra.cpp
LinearAlgebra/Matrix.cpp
LinearAlgebra/MultiVector.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${LinearAlgebra_HEADERS})
  STRING(REPLACE "LinearAlgebra/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/LinearAlgebra"
  COMPONENT ${PACKAGE_NAME}
)