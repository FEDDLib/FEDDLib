set(LinearAlgebra_HEADERS
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
# LinearAlgebra/Map_Xpetra_decl.hpp
# LinearAlgebra/Map_Xpetra_def.hpp
# LinearAlgebra/Map_Xpetra.hpp
LinearAlgebra/Matrix_decl.hpp
LinearAlgebra/Matrix_def.hpp
LinearAlgebra/Matrix.hpp
LinearAlgebra/MultiVector_decl.hpp
LinearAlgebra/MultiVector_def.hpp
LinearAlgebra/MultiVector.hpp
CACHE INTERNAL "")

set(LinearAlgebra_SOURCES
LinearAlgebra/BlockMap.cpp
LinearAlgebra/BlockMatrix.cpp
LinearAlgebra/BlockMultiVector.cpp
LinearAlgebra/Map.cpp
# LinearAlgebra/Map_Xpetra.cpp
LinearAlgebra/Matrix.cpp
LinearAlgebra/MultiVector.cpp
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${LinearAlgebra_HEADERS})
  string(REPLACE "LinearAlgebra/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
