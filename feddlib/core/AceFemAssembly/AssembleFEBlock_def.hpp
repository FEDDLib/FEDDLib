#ifndef AssembleFEBLOCK_DEF_hpp
#define AssembleFEBLOCK_DEF_hpp

namespace FEDD {


template <class SC, class LO, class GO, class NO>
AssembleFEBlock<SC,LO,GO,NO>::AssembleFEBlock(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	dofsChem_ = std::get<2>(this->diskTuple_->at(1)); // Degrees of freedom per node

	numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
	numNodesChem_ = std::get<3>(this->diskTuple_->at(1)); // Number of nodes of element

	dofsElement_ = dofsSolid_*numNodesSolid_ + dofsChem_*numNodesChem_; // "Dimension of return matrix"
	
	dimSystem_ = 2;
/// Element Numbering for triangular elements:

}

/*!

 \brief Assembly Jacobian

*/ 

template <class SC, class LO, class GO, class NO>
void  AssembleFEBlock<SC,LO,GO,NO>::assembleJacobian() {


	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( this->dofsElement_)); // Matrix we fill with entries.

	assembleMonolithicSystem(elementMatrix); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix ; // We init the jacobian matrix with the matrix we just build.
}

/*!

 \brief Assembly Jacobian

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void  AssembleFEBlock<SC,LO,GO,NO>::assembleJacobianBlock(LO i) {


	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( this->dofsElement_)); // Matrix we fill with entries.

	assembleMonolithicSystem(elementMatrix); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix ; // We init the jacobian matrix with the matrix we just build.
}

template <class SC, class LO, class GO, class NO>
void AssembleFEBlock<SC,LO,GO,NO>::assembleMonolithicSystem(SmallMatrixPtr_Type elementMatrix){


};



}
#endif
