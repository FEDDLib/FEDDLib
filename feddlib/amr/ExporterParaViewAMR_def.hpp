#ifndef ExporterParaViewAMR_def_hpp
#define ExporterParaViewAMR_def_hpp

#include "ExporterParaViewAMR_decl.hpp"

/*!
 Definition of ExporterParaView

 @brief  ExporterParaView
 @author Lea Saßmannshausen

 */

namespace FEDD {

template<class SC,class LO,class GO,class NO>
ExporterParaViewAMR<SC,LO,GO,NO>::ExporterParaViewAMR():
ExporterParaView<SC,LO,GO,NO>()
{

}

template<class SC,class LO,class GO,class NO>
void ExporterParaViewAMR<SC,LO,GO,NO>::reSetup(MeshPtr_Type mesh){
    
	this->mesh_ = mesh;
    this->nmbElementsGlob_ = this->mesh_->getNumElementsGlobal();
    this->pointsUnique_ = this->mesh_->getPointsUnique();
    this->nmbPointsGlob_ = this->mesh_->getMapUnique()->getGlobalNumElements();
   
    MapConstPtr_Type elementMap = this->mesh_->getElementMap();
    Teuchos::ArrayView<const GO> nodeElementList = elementMap->getNodeElementList();
    vec_int_Type nodeElementListInteger( this->nmbPointsPerElement_ * nodeElementList.size() );
    int counter=0;

    for (int i=0; i<nodeElementList.size(); i++) {
        for (int j=0; j< this->nmbPointsPerElement_; j++){
            nodeElementListInteger[counter] = (int) this->nmbPointsPerElement_*nodeElementList[i] + j;
            counter++;
        }
    }
    Teuchos::RCP<Epetra_BlockMap> mapElements;

    if (nodeElementListInteger.size()>0)
        mapElements.reset(new Epetra_BlockMap( (int) (this->nmbPointsPerElement_*this->nmbElementsGlob_), (int) nodeElementListInteger.size(), &nodeElementListInteger[0],1, 0, *this->commEpetra_));
    else
        mapElements.reset(new Epetra_BlockMap( (int) (this->nmbPointsPerElement_*this->nmbElementsGlob_), (int) nodeElementListInteger.size(), NULL,1, 0, *this->commEpetra_));
    
    this->elementsHDF_.reset(new Epetra_IntVector(*mapElements));
    
    ElementsPtr_Type elements = this->mesh_->getElementsC();

    counter = 0;
    for (int i=0; i<elements->numberElements(); i++) {
        for (int j=0; j<this->nmbPointsPerElement_; j++) {
            int globalIndex = (int) mesh->getMapRepeated()->getGlobalElement( elements->getElement(i).getNode(j) );
            (*this->elementsHDF_)[counter] = globalIndex;
            counter++;
        }
    }

    Teuchos::ArrayView< const GO > indices = this->mesh_->getMapUnique()->getNodeElementList();
    int* intGlobIDs = new int[indices.size()];
    for (int i=0; i<indices.size(); i++) {
        intGlobIDs[i] = (int) indices[i];
    }
    
    EpetraMapPtr_Type mapEpetra = Teuchos::rcp(new Epetra_Map((int)this->nmbPointsGlob_,indices.size(),intGlobIDs,0,*this->commEpetra_));
    delete [] intGlobIDs;
    
    this->pointsHDF_.reset(new Epetra_MultiVector(*mapEpetra,this->dim_));

    this->updatePoints(); 

	this->redo_=true;
    std::string nameConn = "Connections" + std::to_string(this->timeIndex_);
   	this->writeMeshElements(nameConn);

	
}

template<class SC,class LO,class GO,class NO>
void ExporterParaViewAMR<SC,LO,GO,NO>::updateVariables(MultiVectorConstPtr_Type &u, std::string varName){

    for (int i=0; i<this->variables_.size(); i++) {
		if(this->varNames_[i] == varName){
			this->variables_[i] = u;
	  		if (this->FEType_ == "P0") {
				this->mapUniqueVariables_ = this->mesh_->getElementMap();
			}
			else 
				this->mapUniqueVariables_= this->mesh_->getMapUnique();

			this->nmbExportValuesGlob_ = this->mapUniqueVariables_->getGlobalNumElements();
			Teuchos::ArrayView< const GO > indices = this->mapUniqueVariables_->getNodeElementList();
			int* intGlobIDs = new int[indices.size()];
			for (int j=0; j<indices.size(); j++) {
				intGlobIDs[j] = (int) indices[j];
			}

			EpetraMapPtr_Type mapToStore = Teuchos::rcp(new Epetra_Map( (int) this->mapUniqueVariables_->getGlobalNumElements(), indices.size(), intGlobIDs,0, *this->commEpetra_ ) );

			this->uniqueMaps_[i] =mapToStore;
			delete [] intGlobIDs;
		}
	}

}




}
#endif
