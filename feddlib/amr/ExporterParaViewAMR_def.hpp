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

    // Something different happens to the element List and the elements
    // Probably have the following form:
    // ElementMap :   0     1       2       3       4       5  
    // ->            0 1 2  3 4 5   6 7 8   ...
    // Element GIDs:  0 1 2 3 ...
    // Where the map tells us which IDs belong to which element
    MapConstPtr_Type elementMap = this->mesh_->getElementMap();
    Teuchos::ArrayView<const GO> nodeElementList = elementMap->getNodeElementList(); // Global Ids on this processor
    vec_GO_Type nodeElementListInteger( this->nmbPointsPerElement_ * nodeElementList.size() );
    int counter=0;
    for (int i=0; i<nodeElementList.size(); i++) { // Number of elements
        for (int j=0; j<this->nmbPointsPerElement_; j++){ // number of points per element
            nodeElementListInteger[counter] = (int) this->nmbPointsPerElement_*nodeElementList[i] + j;
            counter++; // from 0 ... to nmbPointsPerElement_*localElements
        }
    }
    Teuchos::ArrayView<GO> globalMapIDs = Teuchos::arrayViewFromVector( nodeElementListInteger);
    MapPtr_Type	mapElements = Teuchos::rcp( new Map_Type((int) (this->nmbPointsPerElement_*this->nmbElementsGlob_), globalMapIDs,elementMap->getIndexBase()*this->nmbPointsPerElement_, this->comm_));
  
    // They contain global IDs of nodes corresponding to 'elements'
    this->elementsHDF_.reset(new MultiVector_Type(mapElements,1));
    
    ElementsPtr_Type elements = this->mesh_->getElementsC();
    counter = 0;
    for (int i=0; i<elements->numberElements(); i++) {
        for (int j=0; j<this->nmbPointsPerElement_; j++) {
            int globalIndex = (int) this->mesh_->getMapRepeated()->getGlobalElement( elements->getElement(i).getNode(j) );
            (this->elementsHDF_->getDataNonConst(0))[counter] = globalIndex;
            counter++;
        }
    }

    this->pointsHDF_.reset(new MultiVector_Type(this->mesh_->getMapUnique(),this->dim_));

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

			this->uniqueMaps_[i] =this->mapUniqueVariables_;
		}
	}

}




}
#endif
