#ifndef MAP_DECL_hpp
#define MAP_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Xpetra_ThyraUtils.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Import.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVectorSpace.hpp>
/*!
 Declaration of Map
 
 @brief  Map
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template < class LO = default_lo, class GO = default_go, class NO = default_no>
class Map {
    
public:
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

    typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef Teuchos::RCP<TpetraMap_Type> TpetraMapPtr_Type;
    typedef Teuchos::RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;
    typedef const TpetraMapConstPtr_Type TpetraMapConstPtrConst_Type;
    
    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;

    typedef Thyra::VectorSpaceBase<default_sc> ThyraVSB_Type;
    typedef Teuchos::RCP<ThyraVSB_Type> ThyraVSBPtr_Type;
    typedef Teuchos::RCP<const ThyraVSB_Type> ThyraVSBConstPtr_Type;

    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<Comm_Type> CommPtr_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    Map();
    
    Map( const TpetraMapConstPtr_Type& tpetraMatPtrIn );
    
    Map( const Map_Type& mapIn );
    
    Map(GO numGlobalElements,
        const Teuchos::ArrayView<const GO> &elementList,
        GO indexBase,
        const CommConstPtr_Type &comm);

    Map(GO numGlobalElements,
        LO numLocalElements,
        GO indexBase,
        const CommConstPtr_Type &comm);

    
    ~Map();

    LO getNodeNumElements() const;
    
    GO getGlobalNumElements() const;
    
    GO getGlobalElement(LO id) const;
    
    LO getLocalElement(GO id) const;

    CommConstPtr_Type getComm() const;

    Teuchos::ArrayView<const GO> getNodeElementList() const;

    GO getIndexBase() const; 

    std::string getUnderlyingLib( ) const;       
    
    MapPtr_Type buildVecFieldMap(UN numDofs, std::string ordering="NodeWise") const;
   
    TpetraMapConstPtr_Type getTpetraMap() const;

    XpetraMapConstPtr_Type getXpetraMap();

    GO getMaxAllGlobalIndex() const;
    
    LO getMaxLocalIndex() const;
    
    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME) const;

    ThyraVSBConstPtr_Type getThyraVectorSpaceBase() const;

    CommPtr_Type getCommNonConst();

    /*!
     @param[in] numFreeProcs: Do not use the last numFreeProcs of MPI communicator in the building process
     */
    Teuchos::RCP<Map<LO,GO,NO> > buildUniqueMap( int numFreeProcs=0 ) const;
    
    Teuchos::RCP<Map<LO,GO,NO> > buildUniqueMap( tuple_intint_Type rankRange ) const;
   
    
private:
    
    TpetraMapConstPtr_Type map_;
    XpetraMapConstPtr_Type mapX_;

};
}

#endif
