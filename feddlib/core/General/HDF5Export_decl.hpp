#ifndef HDF5EXPORT_DECL_hpp
#define HDF5EXPORT_DECL_hpp

#include <fstream>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
// Trilinos
#include <Teuchos_Array.hpp>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_LongLongVector.h>
#include <Epetra_IntVector.h>

#include <EpetraExt_HDF5.h>
#include <hdf5.h>

/*!
 Exporting a HDF5 file

 @brief  HDF5Export
 
 Based on Epetra_Ext it is popssible to export a MultiVector to a HDF5 file with the command 'write'.
 The use must provide the corresponding writemap, to correctly store the (parallely distributed) vector and the file and variable name to stroe it in.
 
 The Structure is as follows:
    The HDF5 file which stores data is set up via Create(filename) -> We have a HDF5 file
    The file can now contain one or more variables with different MultiVectors, but they all correspond to the same map
    When we use more than one variable, the different variables can correspond ie. to different time steps
    -> This can be useful for restarts, where we can then import the different time steps
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class HDF5Export{
public:
    typedef Teuchos::RCP<Epetra_Map> EpetraMapPtr_Type;
    typedef Teuchos::RCP<Epetra_MultiVector> EpetraMVPtr_Type;

    typedef EpetraExt::HDF5 HDF5_Type;
    typedef Teuchos::RCP<HDF5_Type> HDF5Ptr_Type;
    
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    
    /// @brief Constructor for HDF5 Exporter
    /// @param writeMap Map for writing file. Parallel distribution for of the exported multivector. 
    /// @param outputFilename Name for output file
    HDF5Export(MapConstPtr_Type writeMap, std::string outputFilename);

    /// @brief Exporting MultiVector writeVector as HDF5 File with the variable name varName
    /// @param varName Variable name of MultiVector
    /// @param writeVector Vector to be exported, corresponding to writeMap_ 
    void writeVariablesHDF5(std::string varName,MultiVectorConstPtr_Type writeVector);

     /// @brief Closing Exporter
    void closeExporter();

protected:
    
    HDF5Ptr_Type hdf5exporter_;
    CommConstPtr_Type comm_;
    Teuchos::RCP<Epetra_MpiComm> commEpetra_;
    
    // ------------------------
    // READ 
    // ------------------------
    std::string outputFilename_;
    std::vector<std::string>   		varNamesRead_;
    EpetraMapPtr_Type               writeMap_;

};

}

#endif
