#include "TriangleElements.hpp"
/*!
 Definition of Elements
 
 @brief  Elements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
    

SurfaceElements::SurfaceElements():
Elements(),
elementsOfSurfaceGlobal_(0),
elementsOfSurfaceLocal_(0),
surfacesOfElements_(0)
{
//    elements_.reset(new FE_vec_Type());
//    globalIDs_.reset( new vec_GO_Type( 0 ) );
};

SurfaceElements::SurfaceElements( SurfaceElements& elements  ):
Elements( *(Teuchos::rcp_dynamic_cast<Elements_Type> ( Teuchos::rcpFromRef( elements ) ) ) ),
elementsOfSurfaceGlobal_( elements.elementsOfSurfaceGlobal_ ),
elementsOfSurfaceLocal_( elements.elementsOfSurfaceLocal_ ),
surfacesOfElements_(0)
{

};

void SurfaceElements::addSurface( FiniteElement& fe, GO globalID ){
    this->addElement( fe );
    elementsOfSurfaceGlobal_.push_back( vec_GO_Type( 1, globalID ) );
};
void SurfaceElements::setElementsSurface( vec2D_GO_Type& elementsOfSurface ){
    // We assume that elementsOfSurfaceGlobal_ has 1 element per redundant surface.
    // SortUnique was already called for the surfaces, so they are unique now.
    // Here we want to set all elements for an surfaces. We have the information in elementsOfSurfaces, which holds the information of the from unique to redundant surfaces.
    // and elementsOfSurfaceGlobal_, which holds the global element ID.
    
    vec2D_GO_Type elementsOfRedundantSurfaceGlobal =  elementsOfSurfaceGlobal_;

    vec2D_GO_Type newElementsOfSurfaceGlobal( this->numberElements(), vec_GO_Type(0) );
    
    for (int i=0; i<newElementsOfSurfaceGlobal.size(); i++) {
        for (int j=0; j<elementsOfSurface[i].size(); j++) {
            newElementsOfSurfaceGlobal[i].push_back( elementsOfRedundantSurfaceGlobal[ elementsOfSurface[i][j] ][0] );
        }
    }
    
    elementsOfSurfaceGlobal_ = newElementsOfSurfaceGlobal;
    
};

void SurfaceElements::partitionSurfaces( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated ){
   
    typedef Teuchos::OrdinalTraits<LO> OTLO;
    FE_vec_ptr_Type elementsTmp = Teuchos::rcp( new FE_vec_Type ( *elements_ ) );
    
    // Here it is assumed that elementsOfSurfaceGlobal_ is still the list for redundant Surfaces. The correct partitioned and combine list is only set here.
    // We might want to make shift the setup of combined element information to the function sortUniqueAndSetGlobalIDs() in the future
    vec2D_GO_Type elementsOfSurfaceGlobalTmp = elementsOfSurfaceGlobal_;
    
    elementsOfSurfaceGlobal_.resize( 0 );
    elementsOfSurfaceLocal_.resize( 0 );

    this->elements_.reset( new FE_vec_Type ( ) );
    vec_GO_Type globaIDs = *(this->globalIDs_);
    this->globalIDs_.reset( new vec_GO_Type(0) );
    bool setSurface = false;
    for (int i=0; i<elementsTmp->size(); i++) {
        vec_GO_Type elementsOfThisSurfaceGlobal(0);
        vec_LO_Type elementsOfThisSurfaceLocal(0);
        vec_LO_Type localElementIDs(0);
        for (int j=0; j<elementsOfSurfaceGlobalTmp[i].size(); j++) {
            LO idLocal = elementMap->getLocalElement( elementsOfSurfaceGlobalTmp[i][j] );
            // we need to determine which ancestor elements are owned by this rank
            // the rank which owns the first element of all ancestor elements gets the active Surface i
            elementsOfThisSurfaceLocal.push_back( idLocal );
            elementsOfThisSurfaceGlobal.push_back( elementsOfSurfaceGlobalTmp[i][j] );
            if (idLocal != OTLO::invalid())
                setSurface = true;
        }

        if ( setSurface ) {
            
            elementsOfSurfaceGlobal_.push_back( elementsOfThisSurfaceGlobal );
            elementsOfSurfaceLocal_.push_back( elementsOfThisSurfaceLocal );
            
            FiniteElement surface = (*elementsTmp)[i];
            vec_int_Type surfaceVec = { nodeMapRepeated->getLocalElement( surface.getNode(0) ),
                                     nodeMapRepeated->getLocalElement( surface.getNode(1) ),
					nodeMapRepeated->getLocalElement( surface.getNode(2) )
                                    };
            
            this->globalIDs_->push_back( globaIDs[i] );
            FiniteElement surfaceLocal( surfaceVec );
            elements_->push_back( surfaceLocal );
            setSurface= false;
        }
    }
    
};
    
void SurfaceElements::sortUniqueAndSetGlobalIDs( vec2D_GO_Type& combinedElements ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUniqueAndSetGlobalIDs( ) not possible.");
    vec2D_int_Type elementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            elementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
   // elements_.reset( new FE_vec_Type( ) );
    vec_GO_Type elementsOfSurfaceGlobalTmp( elementsOfSurfaceGlobal_.size() );
    for (int i=0; i<elementsOfSurfaceGlobalTmp.size(); i++)
        elementsOfSurfaceGlobalTmp[i] = elementsOfSurfaceGlobal_[i][0];
    
    //we also need to call sort but not unique on the elements belonging to the redundant edges
    makeUniqueWithCombines( elements_, combinedElements, elementsOfSurfaceGlobal_ );
    
//    for (int i=0; i<elementsOfEdgeGlobal_.size(); i++) {
//        std::cout << i << " elementsOfEdgeGlobal_[i].size(): " << elementsOfEdgeGlobal_[i].size() << std::endl;
//    }
    globalIDs_->resize(0);
    for (int i=0; i<numberElements(); i++)
        globalIDs_->push_back(i);
	
    

  
}

const vec_LO_Type& SurfaceElements::getElementsOfSurfaceLocal( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfSurfaceLocal_.size()-1 < i, std::logic_error, "No local elements for this Surface." );

    return elementsOfSurfaceLocal_[i];
}
    
const vec_GO_Type& SurfaceElements::getElementsOfSurfaceGlobal( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfSurfaceGlobal_.size()-1 < i, std::logic_error, "No global elements for this Surface." );
    return elementsOfSurfaceGlobal_[i];
}
        
void SurfaceElements::makeUniqueWithCombines( FE_vec_ptr_Type& elements, vec2D_GO_Type& combinedElements, vec2D_GO_Type& globaIDs )
{
    // We assume that each inner vector of globalIDs has only one element which means that globalIDs can be represented by a vec_GO_Type. However, we want to use vec2D_GO_Type here because we would have to copy the values otherwise in the method sortUniqueAndSetGlobalIDs(), where this method is used
    {
        std::vector<int> index(elements->size());
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        std::sort(index.begin(), index.end(),
                  [&](const int& a, const int& b) {
                      return  (*elements)[a].getVectorNodeList() < (*elements)[b].getVectorNodeList();
                  }
                  );
        elements = sort_from_ref( elements, index );
        globaIDs = sort_from_ref( globaIDs, index );
    }
    {
        std::vector<int> index(elements->size());
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        combinedElements.resize( elements->size() );
        
        auto it = uniqueWithCombines( elements->begin(), elements->end(), combinedElements );
        //std::cout <<"size pre:"<< elements->size()<< std::endl;
        elements->resize( distance( elements->begin(), it ) );
        //std::cout <<"size after:"<< elements->size()<< std::endl;
        combinedElements.resize( elements->size() );
        /*for (int i=0; i<combinedElements.size(); i++) {
            std::cout <<  i << " combinedElementssize:" << combinedElements[i].size() << std::endl;
        }*/
    }
};

SurfaceElements::FE_vec_ptr_Type SurfaceElements::sort_from_ref( FE_vec_ptr_Type& elements,
                                                           std::vector<int> const& reference ) {
    FE_vec_ptr_Type ret = Teuchos::rcp(new FE_vec_Type(elements->size()));
    int const size = elements->size();
    for (long long i = 0; i < size; ++i)
        (*ret)[i] = (*elements)[reference[i]];
    
    return ret;
};


// We assume that vec2D_GO_Type in has only one inner element per outer entry
vec2D_GO_Type SurfaceElements::sort_from_ref( vec2D_GO_Type const& in,
                             std::vector<int> const& reference ) {
    vec2D_GO_Type ret(in.size(), vec_GO_Type(1));
    
    int const size = in.size();
    for (int i = 0; i < size; ++i)
        ret[i][0] = in[reference[i]][0];
    
    return ret;
};




// returns the Surface of Element i
const vec_int_Type SurfaceElements::getSurfacesOfElement( int i ){
	TEUCHOS_TEST_FOR_EXCEPTION( surfacesOfElements_.size()-1 < i, std::logic_error, "No Surface for this Element." );
	return surfacesOfElements_.at(i);
};

// function that matches the edges to the right elements, local indexing
void SurfaceElements::matchSurfacesToElements(MapConstPtr_Type elementMap){

	int numberOfElements=elementMap->getMaxLocalIndex() +1;
	
    vec2D_int_Type newSurfacesOfElements(numberOfElements , vec_int_Type(0) );
	int j1;
	//cout << " Number of Surface Elements " << numberElements()<< " und elementen " << numberOfElements << endl;
	for(int i=0; i< numberElements() ; i ++){
		for(int j=0;j< elementsOfSurfaceLocal_.at(i).size(); j++){
			if(elementsOfSurfaceLocal_.at(i).at(j) != -1){
				j1 = elementsOfSurfaceLocal_.at(i).at(j);
				newSurfacesOfElements.at(j1).push_back(i);
			}
		}
		
		
	}
	/*for(int i=0; i< numberOfElements ; i ++){
	std::cout << " Surfaces of Elements [" << i <<"] :"  ;
		for(int j=0; j<newSurfacesOfElements[i].size() ; j++){
			std::cout << newSurfacesOfElements[i][j] << " " ;
		
		}
	std::cout << std::endl;
	}
	cout << " Matched surface to elements " << endl;*/
    surfacesOfElements_ = newSurfacesOfElements;

};

// Functions to set elementsOfEdgeLocal und elementsOfEdgeGlobal entries remotely 
// (used by meshUnstructured Refinement to complete vector updates)
void SurfaceElements::setElementsOfSurfaceLocalEntry(int index, int entry){
	elementsOfSurfaceLocal_[index].push_back(entry);

}
void SurfaceElements::setElementsOfSurfaceGlobalEntry(int index, int entry){
	elementsOfSurfaceGlobal_[index].push_back(entry);

}
// Parallel Functions for building edgelists, elementsOfEdgeGlobal and elementsOfEdgesLocal
// The difference is, that we dont build the elementsVec and use the elementMap to set the right global IDs to elementsOfEdgeGlobal
void SurfaceElements::sortUniqueAndSetGlobalIDsParallel(MapConstPtr_Type elementMap, vec2D_GO_Type& combinedElements ){

	LO ind;
    for (int i=0; i<elementsOfSurfaceGlobal_.size(); i++){
		ind =  elementMap->getGlobalElement( elementsOfSurfaceGlobal_[i][0]);
        elementsOfSurfaceGlobal_[i][0] = ind ;
	}

	makeUniqueWithCombines( elements_, combinedElements, elementsOfSurfaceGlobal_ );

    globalIDs_->resize(0);
    for (int i=0; i<numberElements(); i++)
        globalIDs_->push_back(i);
	
}
void SurfaceElements::setUpElementsOfSurface( MapConstPtr_Type elementMap, MapConstPtr_Type edgeMap, EdgeElementsPtr_Type edgeElements){

    // Here it is assumed that elementsOfEdgeGlobal_ is still the list for redundant edges. The correct partitioned and combine list is only set here.
    // We might want to make shift the setup of combined element information to the function sortUniqueAndSetGlobalIDs() in the future
    vec2D_GO_Type elementsOfSurfaceGlobalTmp = elementsOfSurfaceGlobal_;

	vec_LO_Type interfaceSurfaceIDs(0);    

    elementsOfSurfaceGlobal_.resize( 0 );
    elementsOfSurfaceLocal_.resize( 0 );

    vec_GO_Type globaIDs = *(this->globalIDs_);
    this->globalIDs_.reset( new vec_GO_Type(0) );
    for (int i=0; i<this->numberElements(); i++) {

        vec_GO_Type elementsOfThisSurfaceGlobal(0);
        vec_LO_Type elementsOfThisSurfaceLocal(0);
        vec_LO_Type localElementIDs(0);
        for (int j=0; j<elementsOfSurfaceGlobalTmp[i].size(); j++) {

            LO idLocal = elementMap->getLocalElement( elementsOfSurfaceGlobalTmp[i][j] );
            // we need to determine which ancestor elements are owned by this rank
            // the rank which owns the first element of all ancestor elements gets the active edge i
            elementsOfThisSurfaceLocal.push_back( idLocal );
            elementsOfThisSurfaceGlobal.push_back( elementsOfSurfaceGlobalTmp[i][j] );

       		 }   
            elementsOfSurfaceGlobal_.push_back( elementsOfThisSurfaceGlobal );
            elementsOfSurfaceLocal_.push_back( elementsOfThisSurfaceLocal );

			if(this->getElement(i).isInterfaceElement()){
				elementsOfSurfaceLocal_[i].push_back(-1);
				interfaceSurfaceIDs.push_back(i);
			}

           // this->globalIDs_->push_back( edgeMap->getGlobalElement((LO)globaIDs[i]) );    
    }
	/*for(int i=0; i<elementsOfSurfaceGlobal_.size(); i++){
		cout << " ElementsOfSurfaceGlobal_["<<i<<"] : " ;
		for(int j=0; j< elementsOfSurfaceGlobal_[i].size(); j++){
			cout<< elementsOfSurfaceGlobal_[i][j] << " " ;
		}
		cout << endl; 
	}*/



	/*
	GO myElement;
	vec2D_int_Type edgesTmp(3,vec_int_Type(2));
	vec_LO_Type triangle;
	vec2D_LO_Type edgeList = edgeElements->getElementList();
	int entry(3);

	for(int i=0; i < interfaceSurfaceIDs.size() ; i++){
		myElement = elementsOfSurfaceGlobal_[interfaceSurfaceIDs[i]][0];
		triangle = this->getElement(interfaceSurfaceIDs[i]).getVectorNodeList();
		edgesTmp[0]={triangle[0],triangle[1]};
		edgesTmp[1]={triangle[0],triangle[2]};
		edgesTmp[2]={triangle[1],triangle[2]};

		vec_GO_Type elementsOfSurfaceTmp(0);
		for(int j=0; j<3; j++){
			auto it1 = find( edgeList.begin(), edgeList.end() ,edgesTmp[j] );
		    entry= distance( edgeList.begin() , it1 );	
			for(int k=0; k<edgeElements->getElementsOfEdgeGlobal( entry ).size(); k++)
				elementsOfSurfaceTmp.push_back(edgeElements->getElementsOfEdgeGlobal( entry )[k]);
		}
		sort(elementsOfSurfaceTmp.begin(),elementsOfSurfaceTmp.end());
		std::cout << " elementsOfSurfaceTmp " ;
		for(int j=0; j< elementsOfSurfaceTmp.size(); j++)
			std::cout << " " << elementsOfSurfaceTmp[j];
		std::cout << " " << std::endl;
	}*/



    
}

}
