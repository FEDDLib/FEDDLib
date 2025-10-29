#ifndef HDF5EXPORT_DEF_hpp
#define HDF5EXPORT_DEF_hpp

#include "HDF5Export_decl.hpp"


namespace FEDD {
 
template<class SC,class LO,class GO,class NO>
HDF5Export<SC,LO,GO,NO>::HDF5Export(MapConstPtr_Type writeMap, std::string outputFilename):
hdf5exporter_(),
comm_()
{

    Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( writeMap->getComm() );

    hdf5exporter_.reset( new HDF5_Type(mpiComm) ); // Building HDF5 Exporter

    hdf5exporter_->create(outputFilename+".h5"); // Creating output file with the 'outoutFilename'

    outputFilename_ = outputFilename; 
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::writeVariablesHDF5(std::string varName, const MultiVectorConstPtr_Type writeVector){

    hdf5exporter_->write(varName,writeVector); // Writing u_export as variable 'varName' in file
    
    if(writeVector->getMap()->getComm()->getRank() == 0 )
        std::cout << " HDF5_Export:: Exporting to file " << outputFilename_ << " with variable name " << varName << std::endl;

    hdf5exporter_->flush();
    
}

template<class SC,class LO,class GO,class NO>
void HDF5Export<SC,LO,GO,NO>::closeExporter(){
    hdf5exporter_->close();
}


}
#endif
