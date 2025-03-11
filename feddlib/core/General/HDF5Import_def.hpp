#ifndef HDF5IMPORT_DEF_hpp
#define HDF5IMPORT_DEF_hpp

#include "HDF5Import_decl.hpp"

/*!
 Importing a HDF5 file

 @brief  HDF5Import
 
 Based on Epetra_Ext it is popssible to import a HDF5 file with the command 'read'.
 The use must provide the corresponding map, to distribute the vector and the correct file and variable name.
 
 */

using namespace std;
namespace FEDD {

template<class SC,class LO,class GO,class NO>
HDF5Import<SC,LO,GO,NO>::HDF5Import(MapConstPtr_Type readMap, std::string inputFilename):
hdf5importer_(),
comm_(),
commEpetra_()
{
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( readMap->getComm() );
    commEpetra_.reset( new Epetra_MpiComm( *mpiComm->getRawMpiComm() ) ); // Communicator for epetra

    // We convert the current map to an eptra map
    Teuchos::ArrayView< const GO > indices = readMap->getNodeElementList(); // Global Ids of map on this processor
    int* intGlobIDs = new int[indices.size()];
    for (int i=0; i<indices.size(); i++) {
        intGlobIDs[i] = (int) indices[i];
    }

    int nmbPointsGlob = readMap->getGlobalNumElements(); // Number of global points

    EpetraMapPtr_Type mapEpetra = Teuchos::rcp(new Epetra_Map((int)nmbPointsGlob,indices.size(),intGlobIDs,0,*commEpetra_)); // Building epetra map

    readMap_ = mapEpetra; // Defining the read map. All different variables that might be contained in the .h5 file must have the same map

    hdf5importer_.reset( new HDF5_Type(*commEpetra_) ); // Build HDF5 importer as Epetra_Ext HDF5 Type
  
    inputFilename_ = inputFilename; //  Name of input file
    hdf5importer_->Open(inputFilename_+".h5"); // We 'open' the file and connect to HDF5 iXpetramporter 

    u_import_Tpetra_.reset(new MultiVector_Type(readMap)); // The general import vector is defined via readMap


}

template<class SC,class LO,class GO,class NO>
typename HDF5Import<SC, LO, GO, NO>::MultiVectorPtr_Type HDF5Import<SC,LO,GO,NO>::readVariablesHDF5(string varName){
    // Testing if requested vairable is contained in file
    TEUCHOS_TEST_FOR_EXCEPTION( !hdf5importer_->IsContained(varName), std::logic_error, "Requested varName: " << varName << " not contained in hdf file "<< inputFilename_ << ".h5.");
    hdf5importer_->Read(varName,*readMap_,u_import_Epetra_); // Reading the variable 'varName' from the file and distribute according to readMap_ to u_import_

    //  Now we write the contents of the eptera multivector u_import to u_import_mv
    Teuchos::ArrayRCP<SC> tmpData = u_import_Tpetra_->getDataNonConst(0);
    for (int i=0; i<u_import_Tpetra_->getLocalLength(); i++) {
        tmpData[i] = u_import_Epetra_->Values()[i]; // this is how you access epetra values of multi vectors
    }

    hdf5importer_->Flush();

    return u_import_Tpetra_;
    
}
template<class SC,class LO,class GO,class NO>
void HDF5Import<SC,LO,GO,NO>::closeImporter(){
    hdf5importer_->Close();
}

}
#endif
