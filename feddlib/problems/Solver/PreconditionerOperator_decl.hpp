#ifndef PreconditionerOperator_DECL_hpp
#define PreconditionerOperator_DECL_hpp
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/problems/problems_config.h"
#include <Thyra_PhysicallyBlockedLinearOpBase.hpp>
#include <Thyra_DefaultProductVectorSpace_decl.hpp>
#include <Teuchos_ConstNonconstObjectContainer.hpp>
#include <Thyra_DefaultMultiVectorProductVectorSpace_decl.hpp>

/*!
 Declaration of PreconditionerOperator
 
 @brief  PreconditionerOperator
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class PreconditionerOperator : virtual public Thyra::PhysicallyBlockedLinearOpBase<SC> {

public:

    PreconditionerOperator();

    /** @name Overridden from PhysicallyBlockedLinearOpBase */
    //@{
    
    /** \brief . */
    virtual void beginBlockFill();
    /** \brief . */
    virtual void beginBlockFill(
                        const int numRowBlocks, const int numColBlocks
                        );
    /** \brief . */
    virtual void beginBlockFill(
                        const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > &productRange,
                        const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > &productDomain
                        );
    /** \brief . */
    virtual bool blockFillIsActive() const;
    /** \brief . */
    virtual bool acceptsBlock(const int i, const int j) const;
    /** \brief . */
    virtual void setNonconstBlock(
                          const int i, const int j,
                          const Teuchos::RCP<Thyra::LinearOpBase<SC> > &block
                          );
    /** \brief . */
    virtual void setBlock(
                  const int i, const int j
                  ,const Teuchos::RCP<const Thyra::LinearOpBase<SC> > &block
                  );
    /** \brief . */
    virtual void endBlockFill();
    /** \brief . */
    virtual void uninitialize();
    
    //@}
    
    /** @name Overridden from BlockedLinearOpBase */
    //@{
    
    /** \brief . */
    Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> >
    productRange() const;
    /** \brief . */
    Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> >
    productDomain() const;
    /** \brief . */
    bool blockExists(const int i, const int j) const;
    /** \brief . */
    bool blockIsConst(const int i, const int j) const;
    /** \brief . */
    Teuchos::RCP<Thyra::LinearOpBase<SC> >
    getNonconstBlock(const int i, const int j);
    /** \brief . */
    Teuchos::RCP<const Thyra::LinearOpBase<SC> >
    getBlock(const int i, const int j) const;
    
    Teuchos::RCP< const Thyra::VectorSpaceBase<SC> > range() const;

    Teuchos::RCP< const Thyra::VectorSpaceBase<SC> > domain() const;

    Teuchos::RCP<const Thyra::LinearOpBase<SC> > clone() const;
    
    std::string description() const;
    
    void describe(
                  Teuchos::FancyOStream &out,
                  const Teuchos::EVerbosityLevel verbLevel
                  ) const;
    
protected:
    
    /** @name Overridden from LinearOpBase */
    //@{
    
    /** \brief Returns <tt>true</tt> only if all constituent operators support
     * <tt>M_trans</tt>.
     */
    bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
    
    /** \brief . */
    virtual void applyImpl(
                   const Thyra::EOpTransp M_trans,
                   const Thyra::MultiVectorBase<SC> &X,
                    const Teuchos::Ptr<Thyra::MultiVectorBase<SC> > &Y,
                   const SC alpha,
                   const SC beta
                   ) const = 0;
    
    //@}           
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductRange_;
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductDomain_;

//    Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace< SC > > defaultProductRange_;
//    Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace< SC > > defaultProductDomain_;

private:
    
    // ///////////////////
    // Private types
    
    typedef Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<SC> > CNCLO;
    typedef Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > > vec_array_t;
    
    template<class SC2>
    struct BlockEntry {
        BlockEntry() : i(-1), j(-1) {}
        BlockEntry( const int i_in, const int j_in, const CNCLO &block_in )
        :i(i_in),j(j_in),block(block_in)
        {}
        int i;
        int j;
        CNCLO block;
    };
    
    // /////////////////////////
    // Private data members
    
    Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > productRange_;
    Teuchos::RCP<const Thyra::ProductVectorSpaceBase<SC> > productDomain_;

    int numRowBlocks_; // M
    int numColBlocks_; // N
    
    std::vector<CNCLO> Ops_; // Final M x N storage
    
    vec_array_t rangeBlocks_;
    vec_array_t domainBlocks_;
    std::vector<BlockEntry<SC> > Ops_stack_; // Temp stack of ops begin filled (if Ops_.size()==0).
    bool blockFillIsActive_;
    
    // ///////////////////////////
    // Private member functions
    
    void resetStorage( const int numRowBlocks, const int numColBlocks );
    void assertBlockFillIsActive(bool) const;
    void assertBlockRowCol(const int i, const int j) const;
    void setBlockSpaces(
                        const int i, const int j, const Thyra::LinearOpBase<SC> &block
                        );
    template<class LinearOpType>
    void setBlockImpl(
                      const int i, const int j,
                      const Teuchos::RCP<LinearOpType> &block
                      );
    void adjustBlockSpaces();
        
};
}
#endif
