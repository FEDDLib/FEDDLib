#ifndef HDF5IMPORT_DECL_hpp
#define HDF5IMPORT_DECL_hpp

#include <fstream>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
// Trilinos
#include <Teuchos_Array.hpp>

#include <hdf5.h>
#include "HDF5Toolbox_decl.hpp"


/*!
 Importing a HDF5 file

 @brief  HDF5Import
 
 Now, we use the HDF5Toolbox as HDF5 export/import toolbox. It is popssible to import a HDF5 file with the command 'read'.
 The user must provide the corresponding map, to distribute the vector and the correct file and variable name.
 
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class HDF5Import {
public:

    typedef HDF5Toolbox<SC,LO,GO,NO> HDF5_Type;
    typedef Teuchos::RCP<HDF5_Type> HDF5Ptr_Type;

    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;


    /// @brief Constructor of HDF import. Here the general setting are defined. An tpetra map build based on the read map
    /// @param readMap Map for reading file. Parallel distribution for the to be imported multivector. 
    /// @param inputFilename Name of input file
    HDF5Import(MapConstPtr_Type readMap, std::string inputFilename);

    /// @brief Reading a variable 'varName' from the inputFile with inputFilename of file type HDF5
    /// @param varName Name of variable contained in file
    /// @return Tpetra formatted multivector distributed as defined with readMap
    MultiVectorPtr_Type readVariablesHDF5(std::string varName);

    // Closing Importer
    void closeImporter();

  protected:
    
    /// @brief HDF5 importer based on HDF5Toolbox
    HDF5Ptr_Type hdf5importer_; 
    CommConstPtr_Type comm_;
    
    // ------------------------
    // READ 
    // ------------------------
    /// @brief Name of input file
    std::string inputFilename_; 
    /// @brief Name of Map of import multivector
    MapConstPtr_Type  readMap_;
    /// @brief Imported file in Tpetra format
    MultiVectorPtr_Type u_import_Tpetra_;

};

}

#endif
