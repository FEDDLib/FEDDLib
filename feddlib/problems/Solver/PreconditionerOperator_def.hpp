#ifndef PreconditionerOperator_DEF_hpp
#define PreconditionerOperator_DEF_hpp


/*!
 Definition of PreconditionerOperator
 
 @brief  PreconditionerOperator, copied from
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
        
// Constructors


template<class SC, class LO, class GO, class NO>
PreconditionerOperator<SC, LO, GO, NO>::PreconditionerOperator()
:numRowBlocks_(0), numColBlocks_(0), blockFillIsActive_(false)
{}


// Overridden from PhysicallyBlockedLinearOpBase


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::beginBlockFill()
{
    uninitialize();
    resetStorage(0,0);
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::beginBlockFill(
                                                    const int numRowBlocks, const int numColBlocks
                                                    )
{
    assertBlockFillIsActive(false);
    uninitialize();
    resetStorage(numRowBlocks,numColBlocks);
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::beginBlockFill(
                                                    const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > &new_productRange
                                                    ,const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > &new_productDomain
                                                    )
{
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"beginBlockFill is used but not implemented!");
}


template<class SC, class LO, class GO, class NO>
bool PreconditionerOperator<SC, LO, GO, NO>::blockFillIsActive() const
{
    return blockFillIsActive_;
}


template<class SC, class LO, class GO, class NO>
bool PreconditionerOperator<SC, LO, GO, NO>::acceptsBlock(
                                                  const int i, const int j
                                                  ) const
{
    assertBlockFillIsActive(true);
    assertBlockRowCol(i,j);
    return true;
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::setNonconstBlock(
                                                      const int i, const int j
                                                      ,const Teuchos::RCP<Thyra::LinearOpBase<SC> > &block
                                                      )
{
    setBlockImpl(i, j, block);
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::setBlock(
                                              const int i, const int j
                                              ,const Teuchos::RCP<const Thyra::LinearOpBase<SC> > &block
                                              )
{
    setBlockImpl(i, j, block);
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::endBlockFill()
{
    
    using Teuchos::as;
    
    assertBlockFillIsActive(true);
    
    // 2009/05/06: rabartl: ToDo: When doing a flexible block fill
    // (Ops_stack_.size() > 0), we need to assert that all of the block rows and
    // columns have been filled in.  I don't think we do that here.
    
    // Get the number of block rows and columns
    if (nonnull(productRange_)) {
        numRowBlocks_ = productRange_->numBlocks();
        numColBlocks_ = productDomain_->numBlocks();
    }
    else {
        numRowBlocks_ = rangeBlocks_.size();
        numColBlocks_ = domainBlocks_.size();
        // NOTE: Above, whether doing a flexible fill or not, all of the blocks
        // must be set in order to have a valid filled operator so this
        // calculation should be correct.
    }
    
    // Assert that all of the block rows and columns have at least one entry if
    // the spaces were not given up front.
//#ifdef TEUCHOS_DEBUG
//    if (is_null(productRange_)) {
//        for (int i = 0; i < numRowBlocks_; ++i) {
//            TEUCHOS_TEST_FOR_EXCEPTION(
//                                       !rangeBlocks_[i].get(), std::logic_error
//                                       ,"PreconditionerOperator<SC>::endBlockFill():"
//                                       " Error, no linear operator block for the i="<<i<<" block row was added"
//                                       " and we can not complete the block fill!"
//                                       );
//        }
//        for(int j = 0; j < numColBlocks_; ++j) {
//            TEUCHOS_TEST_FOR_EXCEPTION(
//                                       !domainBlocks_[j].get(), std::logic_error
//                                       ,"PreconditionerOperator<SC>::endBlockFill():"
//                                       " Error, no linear operator block for the j="
//                                       <<j<<" block column was added"
//                                       " and we can not complete the block fill!"
//                                       );
//        }
//    }
//#endif
    
    // Insert the block LOB objects if doing a flexible fill.
//    if (Ops_stack_.size()) {
//        Ops_.resize(numRowBlocks_*numColBlocks_);
//        for ( int k = 0; k < as<int>(Ops_stack_.size()); ++k ) {
//            const BlockEntry<SC> &block_i_j = Ops_stack_[k];
//            Ops_[numRowBlocks_*block_i_j.j + block_i_j.i] = block_i_j.block;
//        }
//        Ops_stack_.resize(0);
//    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"endBlockFill is used!");

    // Set the product range and domain spaces if not already set
//    if (is_null(productRange_)) {
//        adjustBlockSpaces();
//        defaultProductRange_ = productVectorSpace<SC>(rangeBlocks_());
//        defaultProductDomain_ = productVectorSpace<SC>(domainBlocks_());
//        productRange_ = defaultProductRange_;
//        productDomain_ = defaultProductDomain_;
//    }
//    
//    rangeBlocks_.resize(0);
//    domainBlocks_.resize(0);
    
    blockFillIsActive_ = false;
    
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::uninitialize()
{
    productRange_ = Teuchos::null;
    productDomain_ = Teuchos::null;
    numRowBlocks_ = 0;
    numColBlocks_ = 0;
    Ops_.resize(0);
    Ops_stack_.resize(0);
    rangeBlocks_.resize(0);
    domainBlocks_.resize(0);
    blockFillIsActive_ = false;
}


// Overridden from BlockedLinearOpBase


template<class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::productRange() const
{
    return productRange_;
}


template<class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::productDomain() const
{
    return productDomain_;
}


template<class SC, class LO, class GO, class NO>
bool PreconditionerOperator<SC, LO, GO, NO>::blockExists(
                                                 const int i, const int j
                                                 ) const
{
    assertBlockFillIsActive(false);
    assertBlockRowCol(i,j);
    return true;
}


template<class SC, class LO, class GO, class NO>
bool PreconditionerOperator<SC, LO, GO, NO>::blockIsConst(
                                                  const int i, const int j
                                                  ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
    assertBlockFillIsActive(false);
    assertBlockRowCol(i,j);
    return Ops_[numRowBlocks_*j+i].isConst();
}


template<class SC, class LO, class GO, class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::getNonconstBlock(const int i, const int j)
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
    assertBlockFillIsActive(false);
    assertBlockRowCol(i,j);
    return Ops_[numRowBlocks_*j+i].getNonconstObj();
}


template<class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::LinearOpBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::getBlock(const int i, const int j) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(!blockExists(i,j));
#endif
    assertBlockFillIsActive(false);
    assertBlockRowCol(i,j);
    return Ops_[numRowBlocks_*j+i];
}


// Overridden from LinearOpBase


template<class SC, class LO, class GO, class NO>
Teuchos::RCP< const Thyra::VectorSpaceBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::range() const
{
    return productRange_;
}


template<class SC, class LO, class GO, class NO>
Teuchos::RCP< const Thyra::VectorSpaceBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::domain() const
{
    return productDomain_;
}


template<class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::LinearOpBase<SC> >
PreconditionerOperator<SC, LO, GO, NO>::clone() const
{
    return Teuchos::null; // ToDo: Implement this when needed!
}


// Overridden from Teuchos::Describable


template<class SC, class LO, class GO, class NO>
std::string PreconditionerOperator<SC, LO, GO, NO>::description() const
{
    assertBlockFillIsActive(false);
    std::ostringstream oss;
    oss
    << Teuchos::Describable::description() << "{"
    << "numRowBlocks="<<numRowBlocks_
    << ",numColBlocks="<<numColBlocks_
    << "}";
    return oss.str();
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::describe(
                                              Teuchos::FancyOStream &out_arg
                                              ,const Teuchos::EVerbosityLevel verbLevel
                                              ) const
{
    using Teuchos::rcpFromRef;
    using Teuchos::FancyOStream;
    using Teuchos::OSTab;
    assertBlockFillIsActive(false);
    Teuchos::RCP<FancyOStream> out = rcpFromRef(out_arg);
    OSTab tab1(out);
    switch(verbLevel) {
        case Teuchos::VERB_DEFAULT:
        case Teuchos::VERB_LOW:
            *out << this->description() << std::endl;
            break;
        case Teuchos::VERB_MEDIUM:
        case Teuchos::VERB_HIGH:
        case Teuchos::VERB_EXTREME:
        {
            *out
            << Teuchos::Describable::description() << "{"
            << "rangeDim=" << this->range()->dim()
            << ",domainDim=" << this->domain()->dim()
            << ",numRowBlocks=" << numRowBlocks_
            << ",numColBlocks=" << numColBlocks_
            << "}\n";
            OSTab tab2(out);
            *out
            << "Constituent LinearOpBase objects for M = [ Op[0,0] ..."
            << " ; ... ; ... Op[numRowBlocks-1,numColBlocks-1] ]:\n";
            tab2.incrTab();
            for( int i = 0; i < numRowBlocks_; ++i ) {
                for( int j = 0; j < numColBlocks_; ++j ) {
                    *out << "Op["<<i<<","<<j<<"] = ";
                    Teuchos::RCP<const Thyra::LinearOpBase<SC> >
                    block_i_j = getBlock(i,j);
                    if(block_i_j.get())
                        *out << Teuchos::describe(*getBlock(i,j),verbLevel);
                    else
                        *out << "NULL\n";
                }
            }
            break;
        }
        default:
            TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
    }
}


// protected


// Overridden from LinearOpBase


template<class SC, class LO, class GO, class NO>
bool PreconditionerOperator<SC, LO, GO, NO>::opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    bool supported = true;
    for( int i = 0; i < numRowBlocks_; ++i ) {
        for( int j = 0; j < numColBlocks_; ++j ) {
            Teuchos::RCP<const Thyra::LinearOpBase<SC> >
            block_i_j = getBlock(i,j);
            if( block_i_j.get() && !Thyra::opSupported(*block_i_j,M_trans) )
                supported = false;
        }
    }
    return supported;
}

// private


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::resetStorage(
                                                  const int numRowBlocks, const int numColBlocks
                                                  )
{
    numRowBlocks_ = numRowBlocks;
    numColBlocks_ = numColBlocks;
    Ops_.resize(numRowBlocks_*numColBlocks_);
    if (is_null(productRange_)) {
        rangeBlocks_.resize(numRowBlocks);
        domainBlocks_.resize(numColBlocks);
    }
    blockFillIsActive_ = true;
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::assertBlockFillIsActive(
                                                             bool wantedValue
                                                             ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(!(blockFillIsActive_==wantedValue));
#else
    (void)wantedValue;
#endif
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::assertBlockRowCol(
                                                       const int i, const int j
                                                       ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
                               !( 0 <= i ), std::logic_error
                               ,"Error, i="<<i<<" is invalid!"
                               );
    TEUCHOS_TEST_FOR_EXCEPTION(
                               !( 0 <= j ), std::logic_error
                               ,"Error, j="<<j<<" is invalid!"
                               );
    // Only validate upper range if the number of row and column blocks is
    // fixed!
    if(Ops_.size()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
                                   !( 0 <= i && i < numRowBlocks_ ), std::logic_error
                                   ,"Error, i="<<i<<" does not fall in the range [0,"<<numRowBlocks_-1<<"]!"
                                   );
        TEUCHOS_TEST_FOR_EXCEPTION(
                                   !( 0 <= j && j < numColBlocks_ ), std::logic_error
                                   ,"Error, j="<<j<<" does not fall in the range [0,"<<numColBlocks_-1<<"]!"
                                   );
    }
#else
    (void)i;
    (void)j;
#endif
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::setBlockSpaces(
                                                    const int i, const int j, const Thyra::LinearOpBase<SC> &block
                                                    )
{
    using Teuchos::toString;
    assertBlockFillIsActive(true);
    assertBlockRowCol(i,j);
    
    // Validate that if the vector space block is already set that it is
    // compatible with the block that is being set.
    if( i < numRowBlocks_ && j < numColBlocks_ ) {
#ifdef TEUCHOS_DEBUG
        Teuchos::RCP<const Thyra::VectorSpaceBase<SC> >
        rangeBlock = (
                      productRange_.get()
                      ? productRange_->getBlock(i)
                      : rangeBlocks_[i]
                      ),
        domainBlock = (
                       productDomain_.get()
                       ? productDomain_->getBlock(j)
                       : domainBlocks_[j]
                       );
//        if(rangeBlock.get()) {
//            THYRA_ASSERT_VEC_SPACES_NAMES(
//                                          "PreconditionerOperator<SC>::setBlockSpaces(i,j,block):\n\n"
//                                          "Adding block: " + block.description(),
//                                          *rangeBlock,("(*productRange->getBlock("+toString(i)+"))"),
//                                          *block.range(),("(*block["+toString(i)+","+toString(j)+"].range())")
//                                          );
//        }
//        if(domainBlock.get()) {
//            THYRA_ASSERT_VEC_SPACES_NAMES(
//                                          "PreconditionerOperator<SC>::setBlockSpaces(i,j,block):\n\n"
//                                          "Adding block: " + block.description(),
//                                          *domainBlock,("(*productDomain->getBlock("+toString(j)+"))"),
//                                          *block.domain(),("(*block["+toString(i)+","+toString(j)+"].domain())")
//                                          );
//        }
#endif // TEUCHOS_DEBUG
    }
    
    // Add spaces missing range and domain space blocks if we are doing a
    // flexible fill (otherwise these loops will not be executed)
    for( int k = numRowBlocks_; k <= i; ++k )
        rangeBlocks_.push_back(Teuchos::null);
    for( int k = numColBlocks_; k <= j; ++k )
        domainBlocks_.push_back(Teuchos::null);
    
    // Set the incoming range and domain blocks if not already set
    if(!productRange_.get()) {
        if(!rangeBlocks_[i].get())
            rangeBlocks_[i] = block.range().assert_not_null();
        if(!domainBlocks_[j].get()) {
            domainBlocks_[j] = block.domain().assert_not_null();
        }
    }
    
    // Update the current number of row and columns blocks if doing a flexible
    // fill.
    if(!Ops_.size()) {
        numRowBlocks_ = rangeBlocks_.size();
        numColBlocks_ = domainBlocks_.size();
    }
    
}


template<class SC, class LO, class GO, class NO>
template<class LinearOpType>
void PreconditionerOperator<SC, LO, GO, NO>::setBlockImpl(
                                                  const int i, const int j,
                                                  const Teuchos::RCP<LinearOpType> &block
                                                  )
{
    setBlockSpaces(i, j, *block);
    if (Ops_.size()) {
        // We are doing a fill with a fixed number of row and column blocks so we
        // can just set this.
        Ops_[numRowBlocks_*j+i] = block;
    }
    else {
        // We are doing a flexible fill so add the block to the stack of blocks or
        // replace a block that already exists.
        bool foundBlock = false;
        for( unsigned int k = 0; k < Ops_stack_.size(); ++k ) {
            BlockEntry<SC> &block_i_j = Ops_stack_[k];
            if( block_i_j.i == i && block_i_j.j == j ) {
                block_i_j.block = block;
                foundBlock = true;
                break;
            }
        }
        if(!foundBlock)
            Ops_stack_.push_back(BlockEntry<SC>(i,j,block));
    }
}


template<class SC, class LO, class GO, class NO>
void PreconditionerOperator<SC, LO, GO, NO>::adjustBlockSpaces()
{
    
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_INEQUALITY(Ops_.size(), !=, 0);
#endif
    
    //
    // Loop through the rows and columns looking for rows with mixed
    // single-space range and/or domain spaces on operators and set the single
    // spaces as the block space if it exists.
    //
    // NOTE: Once we get here, we can safely assume that all of the operators
    // are compatible w.r.t. their spaces so if there are rows and/or columns
    // with mixed product and single vector spaces that we can just pick the
    // single vector space for the whole row and/or column.
    //
    
    // Adjust blocks in the range space
    for (int i = 0; i < numRowBlocks_; ++i) {
        for (int j = 0; j < numColBlocks_; ++j) {
            const Teuchos::RCP<const Thyra::LinearOpBase<SC> >
            op_i_j = Ops_[numRowBlocks_*j+i];
            if (is_null(op_i_j))
                continue;
            const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > range_i_j = op_i_j->range();
            if (is_null(Thyra::productVectorSpaceBase<SC>(range_i_j, false))) {
                rangeBlocks_[i] = range_i_j;
                break;
            }
        }
    }
    
    // Adjust blocks in the domain space
    for (int j = 0; j < numColBlocks_; ++j) {
        for (int i = 0; i < numRowBlocks_; ++i) {
            const Teuchos::RCP<const Thyra::LinearOpBase<SC> >
            op_i_j = Ops_[numRowBlocks_*j+i];
            if (is_null(op_i_j))
                continue;
            const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> >
            domain_i_j = op_i_j->domain();
            if (is_null(Thyra::productVectorSpaceBase<SC>(domain_i_j, false))) {
                domainBlocks_[j] = domain_i_j;
                break;
            }
        }
    }
    
}
    
};


#endif
