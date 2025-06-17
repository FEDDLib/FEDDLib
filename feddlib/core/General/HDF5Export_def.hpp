#ifndef HDF5EXPORT_DEF_hpp
#define HDF5EXPORT_DEF_hpp

#include "HDF5Export_decl.hpp"


namespace FEDD {
 
template<class SC,class LO,class GO,class NO>
HDF5Export<SC,LO,GO,NO>::HDF5Export(MapConstPtr_Type writeMap, std::string outputFilename):
hdf5exporter_(),
comm_(),
commEpetra_()
{

    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( writeMap->getComm() );
    commEpetra_.reset( new Epetra_MpiComm( *mpiComm->getRawMpiComm() ) );

    // We convert the current map to an eptra map
    Teuchos::ArrayView< const GO > indices = writeMap->getNodeElementList();
    int* intGlobIDs = new int[indices.size()];
    for (int i=0; i<indices.size(); i++) {
        intGlobIDs[i] = (int) indices[i];
    }

    int nmbPointsGlob = writeMap->getGlobalNumElements(); // Number of global points

    EpetraMapPtr_Type mapEpetra = Teuchos::rcp(new Epetra_Map((int)nmbPointsGlob,indices.size(),intGlobIDs,0,*commEpetra_));

    writeMap_ = mapEpetra; // Defining the read map. All different variables that might be written to the .h5 file have the same map

    hdf5exporter_.reset( new HDF5_Type(*commEpetra_) ); // Building HDF5 Exporter

    hdf5exporter_->Create(outputFilename+".h5"); // Creating output file with the 'outoutFilename'

    outputFilename_ = outputFilename; 
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::writeVariablesHDF5(std::string varName,MultiVectorConstPtr_Type writeVector){

    EpetraMVPtr_Type u_export(new Epetra_MultiVector(*(writeMap_),1)); // Epetra export vector

    TEUCHOS_TEST_FOR_EXCEPTION( std::abs(writeMap_->NumMyElements() - writeVector->getLocalLength()) > 1e-12, std::logic_error, " The local length of map does not match the local mv length. Map and MultiVector are not compatible");

    // We need to write the contents of the writeVector into the Epetra export vector: Convert Xpetra -> Epetra
    Teuchos::ArrayRCP<const SC> tmpData = writeVector->getData(0);
    for (int i=0; i<writeVector->getLocalLength(); i++) {
        u_export->ReplaceMyValue( i, 0, tmpData[i] );
    }

    hdf5exporter_->Write(varName,*u_export); // Writing u_export as variable 'varName' in file
    
    if(writeVector->getMap()->getComm()->getRank() == 0 )
        std::cout << " HDF5_Export:: Exporting to file " << outputFilename_ << " with variable name " << varName << std::endl;

    hdf5exporter_->Flush();
    
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::closeExporter(){
    hdf5exporter_->Close();
}

}
#endif
