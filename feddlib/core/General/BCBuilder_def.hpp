#ifndef BCBuilder_def_hpp
#define BCBuilder_def_hpp

#include "BCBuilder_decl.hpp"
#include "feddlib/core/Utils/FEDDUtils.hpp"
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_ScalarTraitsDecl.hpp>
#include <Xpetra_MatrixFactory.hpp>

/*!
 Definition of BCBuilder
 
 @brief  BCBuilder
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template<class SC,class LO,class GO,class NO>
BCBuilder<SC,LO,GO,NO>::BCBuilder():
    vecBC_func_(),
    vecFlag_(),
    vecBlockID_(),
    vecDomain_(),
    vecBCType_(),
    vecDofs_(),
    vecBC_Parameters_(),
    vecExternalSol_(0),
    resultPtr_(),
    pointPtr_()
#ifdef BCBuilder_TIMER
,SetSystemRowTimer_(Teuchos::TimeMonitor::getNewCounter("BCBuilder: SetSystemRow")),
BlockRowHasDirichletTimer_(Teuchos::TimeMonitor::getNewCounter("BCBuilder: SetSystemRow: BlockHasDiri")),
FindFlagTimer_(Teuchos::TimeMonitor::getNewCounter("BCBuilder: SetSystemRow: BlockHasDiri: FindFlag"))
#endif
{
    
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs){

    vecBC_func_.push_back(funcBC);
    vecFlag_.push_back(flag);
    vecBlockID_.push_back(block);
    vecDomain_.push_back(domain);
    vecBCType_.push_back(type);
    vecDofs_.push_back(dofs);
    vec_dbl_Type dummy(1,0.);
    vecBC_Parameters_.push_back(dummy);
    MultiVectorConstPtr_Type dummyMV;
    vecExternalSol_.push_back(dummyMV);
    
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type &parameter_vec){
    
    vecBC_func_.push_back(funcBC);
    vecFlag_.push_back(flag);
    vecBlockID_.push_back(block);
    vecDomain_.push_back(domain);
    vecBCType_.push_back(type);
    vecDofs_.push_back(dofs);
    vecBC_Parameters_.push_back(parameter_vec);
    MultiVectorConstPtr_Type dummyMV;
    vecExternalSol_.push_back(dummyMV);

}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type &parameter_vec, MultiVectorConstPtr_Type& externalSol){
    
    vecBC_func_.push_back(funcBC);
    vecFlag_.push_back(flag);
    vecBlockID_.push_back(block);
    vecDomain_.push_back(domain);
    vecBCType_.push_back(type);
    vecDofs_.push_back(dofs);
    vecBC_Parameters_.push_back(parameter_vec);
    vecExternalSol_.push_back(externalSol);

}


template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::set(const BlockMatrixPtr_Type &blockMatrix, const BlockMultiVectorPtr_Type &blockMV, double t) const{
    
    setSystem(blockMatrix);
    
    setRHS(blockMV, t);
    
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setRHS(const BlockMultiVectorPtr_Type &blockMV, double t) const{
    vec_dbl_Type result;
    vec_dbl_Type point;
    
    TEUCHOS_TEST_FOR_EXCEPTION( blockMV->getNumVectors()>1, std::runtime_error, "BCBuilder::setRHS() only for getNumVectors == 1.");
    
    int loc = 0;
    for (int block = 0; block < blockMV->size(); block++) { // blocks of RHS vector
        if(blockHasDirichletBC(block, loc)){
                        
            vec_int_ptr_Type bcFlags = vecDomain_.at(loc)->getBCFlagUnique();
            vec2D_dbl_ptr_Type vecPoints = vecDomain_.at(loc)->getPointsUnique();
            int dim = vecDomain_.at(loc)->getDimension();
            int dofs = vecDofs_.at(loc);

            result.resize(dofs,0.);
            point.resize(dim,0.);

            resultPtr_ = Teuchos::rcp(new vec_dbl_Type(dofs,0.));
            pointPtr_ = Teuchos::rcp(new vec_dbl_Type(dim,0.));
            
            Teuchos::ArrayRCP<SC> valuesRHS = blockMV->getBlock(block)->getDataNonConst(0);
            for (int i = 0 ; i < bcFlags->size(); i++) {
                int flag = bcFlags->at(i);
                if(findFlag(flag, block, loc)){ // loc changed to postion corresponding to bcFlag
                    if (!vecBCType_.at(loc).compare("Dirichlet") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Z")||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y_Z")) {
                        if (!vecExternalSol_[loc].is_null()) { // we use the external vector here
                            std::string type = "standard";
                            this->setDirichletBoundaryFromExternal(valuesRHS, i, loc, t, type);
                        }
                        else {
                            for (int d=0; d < dim; d++) {
                                point.at(d) = vecPoints->at(i).at(d);
                            }
                            for (int d=0; d < dofs; d++) {
                                result.at(d) = vecPoints->at(i).at(d);
                            }
                            vecBC_func_.at(loc)(&point.at(0), &result.at(0), t, &(vecBC_Parameters_.at(loc).at(0)));

                            for (UN j=0; j<blockMV->getNumVectors(); j++) {
                                Teuchos::ArrayRCP<SC> values = blockMV->getBlock(block)->getDataNonConst(j);
                                if ( !vecBCType_.at(loc).compare("Dirichlet") ) {
                                    for (int dd=0; dd < dofs; dd++)
                                        values[ dofs * i + dd ] = result.at(dd);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X") ) {
                                    values[ dofs * i + 0 ] = result.at(0);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y") ) {
                                    values[ dofs * i + 1 ] = result.at(1);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Z") ) {
                                    values[ dofs * i + 2 ] = result.at(2);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Y") ) {
                                    values[ dofs * i + 0 ] = result.at(0);
                                    values[ dofs * i + 1 ] = result.at(1);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Z") ) {
                                    values[ dofs * i + 0 ] = result.at(0);
                                    values[ dofs * i + 2 ] = result.at(2);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y_Z") ) {
                                    values[ dofs * i + 1 ] = result.at(1);
                                    values[ dofs * i + 2 ] = result.at(2);
                                }
                            }
                        }
                    }
                }
            }
        }
        //Code below is experimental
        for (int i=0; i<vecBCType_.size(); i++) {
            if( vecBCType_.at(i) == "Neumann" &&  vecBlockID_.at(i)==block){
                DomainPtr_Type domain = vecDomain_.at(i);
                FEFacPtr_Type feFactory = Teuchos::rcp( new FEFac_Type() );
                feFactory->addFE(domain);
                
                std::vector<SC> funcParameter(3 , 0.);//0: order, 1:time, 2:surface flag (no specific flag is used in function above)
                funcParameter[1] = t;
                funcParameter[2] = vecFlag_.at(i);
                int dim = domain->getDimension();
                int dofs = vecDofs_.at(i);
                MultiVectorPtr_Type a;
                MultiVectorPtr_Type aUnique;
                if (dofs>1){
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Vector", vecBC_func_.at(i), funcParameter );
                }
                else{
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Scalar", vecBC_func_.at(i), funcParameter );
                }
                aUnique->exportFromVector( a, false, "Add" );
                MultiVectorPtr_Type mv = blockMV->getBlockNonConst(block);
                mv->update(1.,*aUnique,1.);
            }
        }
    }
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setDirichletBoundaryFromExternal(Teuchos::ArrayRCP<SC>& values/*values will be set to this vector*/, LO index, int loc, double time, std::string type, Teuchos::ArrayRCP<SC> valuesSubstract) const{
    
    Teuchos::ArrayRCP<SC> valuesEx = vecExternalSol_[loc]->getDataNonConst(0);//We assume that our Multivector only has one column
        
    (*pointPtr_)[0] = valuesEx[index]; // laplace solution in current point
                
    vecBC_func_.at(loc)( &((*pointPtr_)[0]), &((*resultPtr_)[0]), time, &(vecBC_Parameters_.at(loc).at(0)));
        
    int dofs = resultPtr_->size();
    if ( !vecBCType_.at(loc).compare("Dirichlet") ) {
        if (type == "standard") {
            for (int dd=0; dd < dofs; dd++){
                values[ dofs * index + dd ] = (*resultPtr_)[dd];
			}
        }
        else if (type == "BCMinusVec"){
            for (int dd=0; dd < dofs; dd++)
                values[ dofs * index + dd ] = (*resultPtr_)[dd] - valuesSubstract[ dofs * index + dd ];
        }
        else if (type == "VecMinusBC"){
            for (int dd=0; dd < dofs; dd++)
                values[ dofs * index + dd ] = valuesSubstract[ dofs * index + dd ] - (*resultPtr_)[dd];
        }
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "BCBuilder::setDirichletBoundaryFromExternal only implemented for full Dirichlet.");
    }
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setBCMinusVector(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t ) const{
    
    TEUCHOS_TEST_FOR_EXCEPTION( outBlockMV->getNumVectors()>1, std::runtime_error, "BCBuilder::setRHS() only for getNumVectors == 1.");
    
    vec_dbl_Type result;
    vec_dbl_Type point;
    int loc = 0;

    for (int block = 0; block < outBlockMV->size(); block++) { // blocks of RHS vector
        if(blockHasDirichletBC(block, loc)){
            vec_int_ptr_Type     bcFlags = vecDomain_.at(loc)->getBCFlagUnique();
            vec2D_dbl_ptr_Type     vecPoints = vecDomain_.at(loc)->getPointsUnique();
            int dim = vecDomain_.at(loc)->getDimension();
            int dofs = vecDofs_.at(loc);
            
            result.resize(dofs,0.);
            point.resize(dim,0.);
            
            resultPtr_ = Teuchos::rcp(new vec_dbl_Type(dofs,0.));
            pointPtr_ = Teuchos::rcp(new vec_dbl_Type(dim,0.));
            
            Teuchos::ArrayRCP<SC> valuesRHS = outBlockMV->getBlock(block)->getDataNonConst(0);
            Teuchos::ArrayRCP<SC> valuesSubstract = substractBlockMV->getBlock(block)->getDataNonConst(0);
            
            for (int i = 0 ; i < bcFlags->size(); i++) {
                int flag = bcFlags->at(i);
                if(findFlag(flag,block, loc)){ // loc changed to postion corresponding to bcFlag
                    if (!vecBCType_.at(loc).compare("Dirichlet") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Z")||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y_Z")) {
                        if (!vecExternalSol_[loc].is_null()) { // we use the external vector here
                            std::string type = "BCMinusVec";
                            this->setDirichletBoundaryFromExternal(valuesRHS, i, loc, t, type, valuesSubstract);
                        }
                        else{
                            for (int d=0; d < dim; d++) {
                                point.at(d) = vecPoints->at(i).at(d);
                            }
                            for (int d=0; d < dofs; d++) {
                                result.at(d) = vecPoints->at(i).at(d);
                            }
                            vecBC_func_.at(loc)(&point.at(0), &result.at(0), t, &(vecBC_Parameters_.at(loc).at(0)));
                            
                            for (UN j=0; j<outBlockMV->getNumVectors(); j++) {
                                Teuchos::ArrayRCP<SC> values = outBlockMV->getBlock(block)->getDataNonConst(j);
                                Teuchos::ArrayRCP<const SC> substractValues = substractBlockMV->getBlock(block)->getData(j);
                                if ( !vecBCType_.at(loc).compare("Dirichlet") ) {
                                    for (int dd=0; dd < dofs; dd++)
                                        values[ dofs * i + dd ] = result.at(dd) - substractValues[ dofs * i + dd ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X") ) {
                                    values[ dofs * i + 0 ] = result.at(0) - substractValues[ dofs * i + 0 ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y") ) {
                                    values[ dofs * i + 1 ] = result.at(1) - substractValues[ dofs * i + 1 ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Z") ) {
                                    values[ dofs * i + 2 ] = result.at(2) - substractValues[ dofs * i + 2 ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Y") ) {
                                    values[ dofs * i + 0 ] = result.at(0) - substractValues[ dofs * i + 0 ];
                                    values[ dofs * i + 1 ] = result.at(1) - substractValues[ dofs * i + 1 ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Z") ) {
                                    values[ dofs * i + 0 ] = result.at(0) - substractValues[ dofs * i + 0 ];
                                    values[ dofs * i + 2 ] = result.at(2) - substractValues[ dofs * i + 2 ];
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y_Z") ) {
                                    values[ dofs * i + 1 ] = result.at(1) - substractValues[ dofs * i + 1 ];
                                    values[ dofs * i + 2 ] = result.at(2) - substractValues[ dofs * i + 2 ];
                                }
                            }
                        }
                    }
                }
            }
        }
        //Code below is experimental
        for (int i=0; i<vecBCType_.size(); i++) {
            if( vecBCType_.at(i) == "Neumann" &&  vecBlockID_.at(i)==block){
                DomainPtr_Type domain = vecDomain_.at(i);
                FEFacPtr_Type feFactory = Teuchos::rcp( new FEFac_Type() );
                feFactory->addFE(domain);
                
                std::vector<SC> funcParameter(3 , 0.);//0: order, 1:time, 2:surface flag (no specific flag is used in function above)
                funcParameter[1] = t;
                funcParameter[2] = vecFlag_.at(i);
                int dim = domain->getDimension();
                int dofs = vecDofs_.at(i);
                MultiVectorPtr_Type a;
                MultiVectorPtr_Type aUnique;
                if (dofs>1){
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Vector", vecBC_func_.at(i), funcParameter );
                }
                else{
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Scalar", vecBC_func_.at(i), funcParameter );
                }
                aUnique->exportFromVector( a, false, "Add" );
                MultiVectorPtr_Type mv = outBlockMV->getBlockNonConst(block);
                mv->update(1.,*aUnique,1.);
                MultiVectorPtr_Type mvMinus = substractBlockMV->getBlockNonConst(block);
                mv->update(-1.,*mvMinus,1.);
            }
        }
    }
}
    
template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setVectorMinusBC(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t) const{
    vec_dbl_Type result;
    vec_dbl_Type point;
    int loc = 0;

    TEUCHOS_TEST_FOR_EXCEPTION( outBlockMV->getNumVectors()>1, std::runtime_error, "BCBuilder::setRHS() only for getNumVectors == 1.");
    
    for (int block = 0; block < outBlockMV->size(); block++) { // blocks of RHS vector
        if(blockHasDirichletBC(block, loc)){
            vec_int_ptr_Type     bcFlags = vecDomain_.at(loc)->getBCFlagUnique();
            vec2D_dbl_ptr_Type     vecPoints = vecDomain_.at(loc)->getPointsUnique();
            int dim = vecDomain_.at(loc)->getDimension();
            int dofs = vecDofs_.at(loc);
            
            result.resize(dofs,0.);
            point.resize(dim,0.);

            resultPtr_ = Teuchos::rcp(new vec_dbl_Type(dofs,0.));
            pointPtr_ = Teuchos::rcp(new vec_dbl_Type(dim,0.));
            
            Teuchos::ArrayRCP<SC> valuesRHS = outBlockMV->getBlock(block)->getDataNonConst(0);
            Teuchos::ArrayRCP<SC> valuesSubstract = substractBlockMV->getBlock(block)->getDataNonConst(0);
            for (int i = 0 ; i < bcFlags->size(); i++) {
                int flag = bcFlags->at(i);
                if(findFlag(flag,block, loc)){ // loc changed to postion corresponding to bcFlag
                    if (!vecBCType_.at(loc).compare("Dirichlet") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y_Z")) {
                        if (!vecExternalSol_[loc].is_null()) { // we use the external vector here
                            std::string type = "VecMinusBC";
                            this->setDirichletBoundaryFromExternal(valuesRHS, i, loc, t, type, valuesSubstract);
                        }
                        else{
                            for (int d=0; d < dim; d++) {
                                point.at(d) = vecPoints->at(i).at(d);
                            }
                            for (int d=0; d < dofs; d++) {
                                result.at(d) = vecPoints->at(i).at(d);
                            }
                            vecBC_func_.at(loc)(&point.at(0), &result.at(0), t, &(vecBC_Parameters_.at(loc).at(0)));
                            
                            for (UN j=0; j<outBlockMV->getNumVectors(); j++) {
                                Teuchos::ArrayRCP<SC> values = outBlockMV->getBlock(block)->getDataNonConst(j);
                                Teuchos::ArrayRCP<const SC> substractValues = substractBlockMV->getBlock(block)->getData(j);
                                if ( !vecBCType_.at(loc).compare("Dirichlet") ) {
                                    for (int dd=0; dd < dofs; dd++)
                                        values[ dofs * i + dd ] = substractValues[ dofs * i + dd ] - result.at(dd);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X") ) {
                                    values[ dofs * i + 0 ] = substractValues[ dofs * i + 0 ] - result.at(0);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y") ) {
                                    values[ dofs * i + 1 ] = substractValues[ dofs * i + 1 ] - result.at(1);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Z") ) {
                                    values[ dofs * i + 2 ] = substractValues[ dofs * i + 2 ] - result.at(2);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Y") ) {
                                    values[ dofs * i + 0 ] = substractValues[ dofs * i + 0 ] - result.at(0);
                                    values[ dofs * i + 1 ] = substractValues[ dofs * i + 1 ] - result.at(1);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Z") ) {
                                    values[ dofs * i + 0 ] = substractValues[ dofs * i + 0 ] - result.at(0);
                                    values[ dofs * i + 2 ] = substractValues[ dofs * i + 2 ] - result.at(2);
                                }
                                else if ( !vecBCType_.at(loc).compare("Dirichlet_Y_Z") ) {
                                    values[ dofs * i + 1 ] = substractValues[ dofs * i + 1 ] - result.at(1);
                                    values[ dofs * i + 2 ] = substractValues[ dofs * i + 2 ] - result.at(2);
                                }
                            }
                        }
                    }
                }
            }
        }
        //Code below is experimental
        for (int i=0; i<vecBCType_.size(); i++) {
            if( vecBCType_.at(i) == "Neumann" &&  vecBlockID_.at(i)==block){
                DomainPtr_Type domain = vecDomain_.at(i);
                FEFacPtr_Type feFactory = Teuchos::rcp( new FEFac_Type() );
                feFactory->addFE(domain);
                
                std::vector<SC> funcParameter(3 , 0.);//0: order, 1:time, 2:surface flag (no specific flag is used in function above)
                funcParameter[1] = t;
                funcParameter[2] = vecFlag_.at(i);
                int dim = domain->getDimension();
                int dofs = vecDofs_.at(i);
                MultiVectorPtr_Type a;
                MultiVectorPtr_Type aUnique;
                if (dofs>1){
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapVecFieldUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Vector", vecBC_func_.at(i), funcParameter );
                }
                else{
                    a = Teuchos::rcp(new MultiVector_Type ( domain->getMapRepeated() ) );
                    aUnique = Teuchos::rcp(new MultiVector_Type ( domain->getMapUnique() ) );
                    feFactory->assemblySurfaceIntegralFlag( dim, domain->getFEType(), a, "Scalar", vecBC_func_.at(i), funcParameter );
                }
                aUnique->exportFromVector( a, false, "Add" );
                MultiVectorPtr_Type mv = outBlockMV->getBlockNonConst(block);
                mv->update(1.,*aUnique,1.);
                MultiVectorPtr_Type mvMinus = substractBlockMV->getBlockNonConst(block);
                mv->update(1.,*mvMinus,-1.);
            }
        }
    }
}
    
template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setAllDirichletZero(const BlockMultiVectorPtr_Type &blockMV) const{
    vec_dbl_Type result;
    int loc = 0;
    for (int block = 0; block < blockMV->size(); block++) { // blocks of RHS vector
        if(blockHasDirichletBC(block, loc)){
            vec_int_ptr_Type 	bcFlags = vecDomain_.at(loc)->getBCFlagUnique();
            vec2D_dbl_ptr_Type 	vecPoints = vecDomain_.at(loc)->getPointsUnique();
            int dim = vecDomain_.at(loc)->getDimension();
            int dofs = vecDofs_.at(loc);
            
            result.resize(dofs,0.);
            
            for (int i = 0 ; i < bcFlags->size(); i++) {
                int flag = bcFlags->at(i);
                if(findFlag(flag,block, loc)){ // loc changed to postion corresponding to bcFlag
                    if (    !vecBCType_.at(loc).compare("Dirichlet") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Y") ||
                        !vecBCType_.at(loc).compare("Dirichlet_X_Z") ||
                        !vecBCType_.at(loc).compare("Dirichlet_Y_Z")
                        ){
                        
                        for (UN j=0; j<blockMV->getNumVectors(); j++) {
                            Teuchos::ArrayRCP<SC> values = blockMV->getBlock(block)->getDataNonConst(j);
                            if ( !vecBCType_.at(loc).compare("Dirichlet") ) {
                                for (int dd=0; dd < dofs; dd++)
                                    values[ dofs * i + dd ] = result[dd];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_X") ) {
                                values[ dofs * i + 0 ] = result[0];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_Y") ) {
                                values[ dofs * i + 1 ] = result[1];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_Z") ) {
                                values[ dofs * i + 2 ] = result[2];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Y") ) {
                                values[ dofs * i + 0 ] = result[0];
                                values[ dofs * i + 1 ] = result[1];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_X_Z") ) {
                                values[ dofs * i + 0 ] = result[0];
                                values[ dofs * i + 2 ] = result[2];
                            }
                            else if ( !vecBCType_.at(loc).compare("Dirichlet_Y_Z") ) {
                                values[ dofs * i + 1 ] = result[1];
                                values[ dofs * i + 2 ] = result[2];
                            }
                        }
                    }
                }
            }
        }
    }
}

    
    
template<class SC,class LO,class GO,class NO>
bool BCBuilder<SC,LO,GO,NO>::blockHasDirichletBC(int block) const{
    int dummyloc;

    return blockHasDirichletBC(block, dummyloc);
}

template<class SC,class LO,class GO,class NO>
bool BCBuilder<SC,LO,GO,NO>::blockHasDirichletBC(int block, int &loc) const{

    bool hasBC = false;
    bool found = false;
    loc = -1;
    std::vector<int>::const_iterator it = vecBlockID_.begin();
    while (it!=vecBlockID_.end() && !found) {
        it = std::find(it,vecBlockID_.end(),block);
        if (it!=vecBlockID_.end()) {
            loc = distance(vecBlockID_.begin(),it);
            it++;
            if (!vecBCType_.at(loc).compare("Dirichlet") ||
                !vecBCType_.at(loc).compare("Dirichlet_X") ||
                !vecBCType_.at(loc).compare("Dirichlet_Y") ||
                !vecBCType_.at(loc).compare("Dirichlet_Z") ||
                !vecBCType_.at(loc).compare("Dirichlet_X_Y") ||
                !vecBCType_.at(loc).compare("Dirichlet_X_Z") ||
                !vecBCType_.at(loc).compare("Dirichlet_Y_Z")
                ) {
                found = true;
                hasBC = true;
            }
        }
    }
    
    return hasBC;
}

template<class SC,class LO,class GO,class NO>
bool BCBuilder<SC,LO,GO,NO>::findFlag(LO flag, int block, int &loc) const{

#ifdef BCBuilder_TIMER
    Teuchos::TimeMonitor FindFlagMonitor(*FindFlagTimer_);
#endif

    bool hasFlag = false;
    bool found = false;
    loc = -1;
    std::vector<int>::const_iterator it = vecFlag_.begin();
    while (it!=vecFlag_.end() && !found) {
        it = std::find(it,vecFlag_.end(),flag);
        if (it!=vecFlag_.end()) {
            loc = distance(vecFlag_.begin(),it);
            if (vecBlockID_.at(loc)==block) {
                hasFlag = true;
                found 	= true;
            }
            else{
                it++;
            }
        }
    }
    
    return hasFlag;
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setSystem(const BlockMatrixPtr_Type &blockMatrix) const{

    UN numBlocks = blockMatrix->size();
    int loc;
    for (UN blockRow = 0; blockRow < numBlocks; blockRow++) {

#ifdef BCBuilder_TIMER
        Teuchos::TimeMonitor SetSystemRowMonitor(*SetSystemRowTimer_);
#endif
        bool boolBlockHasDirichlet = false;
        
        {
#ifdef BCBuilder_TIMER
            Teuchos::TimeMonitor BlockRowHasDirichletMonitor(*BlockRowHasDirichletTimer_);
#endif

            boolBlockHasDirichlet = blockHasDirichletBC(blockRow,loc);
        }

        if (boolBlockHasDirichlet){
            for (UN blockCol = 0; blockCol < numBlocks ; blockCol++) {
                if (blockMatrix->blockExists(blockRow, blockCol)) {
                    MatrixPtr_Type matrix = blockMatrix->getBlock(blockRow, blockCol);
                    setDirichletBC(matrix, loc, blockRow, blockRow == blockCol);
                } else if (blockRow == blockCol) {
                    // If the current block is a diagonal block but does not exist we need to create it and insert the
                    // value 1 in all diagonal entries of the block corresponding to a Dirichlet boundary condition. We
                    // only need max. one entry per row.

                    // To get the correct domain from vecDomain_ we need to know the index of any (in this case the
                    // first) entry that was done with addBC for this block
                    const auto it = std::find(this->vecBlockID_.begin(), this->vecBlockID_.end(), blockRow);
                    cout << " loc " << loc << " it - something " << it - this->vecBlockID_.begin() << endl;
                    auto matrixMap = this->vecDomain_.at(loc)->getMapUnique();
                    // // Use Xpetra::MatrixFactory to build a matrix with known row and column maps
                    auto tpetraMatrix = Teuchos::rcp(new Tpetra::CrsMatrix<SC,LO,GO,NO>(matrixMap->getTpetraMap(), matrixMap->getTpetraMap(), 1)); // Tpetra::MatrixFactory<SC, LO, GO, NO>::Build(matrixMap->getTpetraMap(), matrixMap->getTpetraMap(), 1);

                    MatrixPtr_Type matrix = Teuchos::rcp(new Matrix_Type(tpetraMatrix));

                    // Fill the matrix with zeros on the diagonal
                    Teuchos::Array<LO> colIndex(1, 0);
                    Teuchos::Array<SC> val(1, Teuchos::ScalarTraits<SC>::zero());
                    for (auto i = 0; i < matrixMap->getNodeNumElements(); i++) {
                        colIndex[0] = i;
                        // Since the matrix has a column map we can use insertLocalValues(). This is more efficient!
                        matrix->insertLocalValues(i, colIndex, val);
                    }
                    matrix->fillComplete(matrixMap, matrixMap);
                    setDirichletBC(matrix, loc, blockRow, true);
                    blockMatrix->addBlock(matrix, blockRow, blockCol);
                }
            }
        }
    }
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setDirichletBC(const MatrixPtr_Type &matrix, int loc, int blockRow, bool isDiagonalBlock) const{
    
    matrix->resumeFill();
    bool isDirichlet;
    UN dofsPerNode = vecDofs_.at(loc);
    vec_int_ptr_Type bcFlags = vecDomain_.at(loc)->getBCFlagUnique();

    for (LO i=0; i<bcFlags->size(); i++) { // flags are for nodes.
        int flag = bcFlags->at(i);
        isDirichlet = false;
        int locThisFlag;

        if( findFlag(flag, blockRow, locThisFlag) ){

            if( !vecBCType_.at(locThisFlag).compare("Dirichlet") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_X") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_Y") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_Z") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_X_Y") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_X_Z") ||
                !vecBCType_.at(locThisFlag).compare("Dirichlet_Y_Z") )
                isDirichlet = true;
            if (isDirichlet) {

                if (isDiagonalBlock)
                    setLocalRowOne(matrix, i, dofsPerNode, locThisFlag );
                else
                    setLocalRowZero(matrix, i, dofsPerNode, locThisFlag );
            }
        }
    }
    matrix->fillComplete( matrix->getMap("domain"), matrix->getMap("range") );
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setLocalRowOne(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const{
    
    Teuchos::ArrayView<const SC> valuesOld;
    Teuchos::ArrayView<const LO> indices;

    MapConstPtr_Type colMap = matrix->getMap("col");
    LO localDof = (LO) (dofsPerNode * localNode);
    for (UN dof=0; dof<dofsPerNode; dof++) {
        if ( vecBCType_.at(loc) == "Dirichlet" ||
            (vecBCType_.at(loc) == "Dirichlet_X" && dof==0 ) ||
            (vecBCType_.at(loc) == "Dirichlet_Y" && dof==1 ) ||
            (vecBCType_.at(loc) == "Dirichlet_Z" && dof==2 ) ||
            (vecBCType_.at(loc) == "Dirichlet_X_Y" && dof!=2 )||
            (vecBCType_.at(loc) == "Dirichlet_X_Z" && dof!=1 )||
            (vecBCType_.at(loc) == "Dirichlet_Y_Z" && dof!=0 )
            ) {
           // cout << " Setting Dirichlet Row 1 for node " << localDof << " of type " << vecBCType_.at(loc)  <<endl;
            GO globalDof = matrix->getMap()->getGlobalElement( localDof );
            matrix->getLocalRowView(localDof, indices, valuesOld);
            Teuchos::Array<SC> values( valuesOld.size(), Teuchos::ScalarTraits<SC>::zero() );
            bool setOne = false;
            for (UN j=0; j<indices.size() && !setOne; j++) {
                if ( colMap->getGlobalElement( indices[j] )  == globalDof ){
                    values[j] = Teuchos::ScalarTraits<SC>::one();
                    setOne = true;
                }
            }
            matrix->replaceLocalValues(localDof, indices(), values());
        }
        localDof++;
    }
    
}

template<class SC,class LO,class GO,class NO>
void BCBuilder<SC,LO,GO,NO>::setLocalRowZero(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const{
    
    Teuchos::ArrayView<const SC> valuesOld;
    Teuchos::ArrayView<const LO> indices;
    LO localDof = (LO) (dofsPerNode * localNode);
    for (UN dof=0; dof<dofsPerNode; dof++) {
        if ( vecBCType_.at(loc) == "Dirichlet" ||
        (vecBCType_.at(loc) == "Dirichlet_X" && dof==0 ) ||
        (vecBCType_.at(loc) == "Dirichlet_Y" && dof==1 ) ||
        (vecBCType_.at(loc) == "Dirichlet_Z" && dof==2 ) ||
        (vecBCType_.at(loc) == "Dirichlet_X_Y" && dof!=2 )||
        (vecBCType_.at(loc) == "Dirichlet_X_Z" && dof!=1 )||
        (vecBCType_.at(loc) == "Dirichlet_Y_Z" && dof!=0 )
        ) {
            matrix->getLocalRowView(localDof, indices, valuesOld);
            Teuchos::Array<SC> values( valuesOld.size(), Teuchos::ScalarTraits<SC>::zero() );
            matrix->replaceLocalValues(localDof, indices(), values());
        }
        localDof++;
    }
}



template<class SC,class LO,class GO,class NO>
int BCBuilder<SC,LO,GO,NO>::dofsPerNodeAtBlock(int block) {
    std::vector<int>::iterator it;
    it = find(vecBlockID_.begin(),vecBlockID_.end(),block);
#ifdef ASSERTS_WARNINGS
    MYASSERT(it!=vecBlockID_.end(),"Requested dof information is not known to BCBuilder class.");
#endif
    
    return vecDofs_.at(distance(vecBlockID_.begin(),it));
}
}

#endif
