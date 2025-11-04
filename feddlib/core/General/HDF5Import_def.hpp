#ifndef HDF5IMPORT_DEF_hpp
#define HDF5IMPORT_DEF_hpp

/*!
 Importing a HDF5 file

 @brief  HDF5Import
 
 Now, we use the HDF5Toolbox as HDF5 export/import toolbox. It is possible to import a HDF5 file with the command 'read'.
 The use must provide the corresponding map, to distribute the vector and the correct file and variable name.
 
 */

namespace FEDD {

template<class SC,class LO,class GO,class NO>
HDF5Import<SC,LO,GO,NO>::HDF5Import(MapConstPtr_Type readMap, std::string inputFilename):
hdf5importer_(),
comm_()
{
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( readMap->getComm() );
    comm_ = mpiComm;
 
    readMap_ = readMap; // Defining the read map. All different variables that might be contained in the .h5 file must have the same map

    hdf5importer_.reset( new HDF5_Type(comm_) ); // Build HDF5 importer as HDF5Toolbox Type
  
    inputFilename_ = inputFilename; //  Name of input file
    hdf5importer_->open(inputFilename_+".h5"); // We 'open' the file and connect to HDF5 importer 

    u_import_Tpetra_.reset(new MultiVector_Type(readMap)); // The general import vector is defined via readMap

}

template<class SC,class LO,class GO,class NO>
typename HDF5Import<SC, LO, GO, NO>::MultiVectorPtr_Type HDF5Import<SC,LO,GO,NO>::readVariablesHDF5(std::string varName){
    // Testing if requested vairable is contained in file
    TEUCHOS_TEST_FOR_EXCEPTION( !hdf5importer_->isContained(varName), std::logic_error, "Requested varName: " << varName << " not contained in hdf file "<< inputFilename_ << ".h5.");
    hdf5importer_->read(varName,readMap_,u_import_Tpetra_); // Reading the variable 'varName' from the file and distribute according to readMap_ to u_import_

    return u_import_Tpetra_;
    
}
template<class SC,class LO,class GO,class NO>
void HDF5Import<SC,LO,GO,NO>::closeImporter(){
    hdf5importer_->close();
}

}
#endif
