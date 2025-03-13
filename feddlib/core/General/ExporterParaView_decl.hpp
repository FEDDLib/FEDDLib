#ifndef ExporterParaView_DECL_hpp
#define ExporterParaView_DECL_hpp

#include <fstream>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/Mesh.hpp"
// Trilinos
#include <Teuchos_Array.hpp>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_LongLongVector.h>
#include <Epetra_IntVector.h>

#include <EpetraExt_HDF5.h>
#include <hdf5.h>

/*!
 Declaration of ExporterParaView
 
 @brief  ExporterParaView
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 
 Moving mesh are possible, here connections are constant but coordinates can be shifted between exports.
 Use set the following paramter in your parameter list in order to export updated mesh points sublist("Exporter")->get("Write new mesh") = true
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class ExporterParaView {
public:
    typedef std::vector<double>										vec_dbl;
    typedef std::vector<std::vector<double> >						vec2D_dbl;
    typedef std::vector<std::vector<int> >							vec2D_int;
    typedef std::vector<std::vector<long long> >					vec2D_longlong;
    typedef Teuchos::RCP<std::vector<int> >							vec_int_ptr;
    typedef Teuchos::RCP<std::vector<long long> >					vec_longlong_ptr;
    typedef Teuchos::RCP<vec_dbl>									vec_dbl_ptr;
    typedef Teuchos::RCP<std::vector<std::vector<double> > >     	vec2D_dbl_ptr;
    typedef Teuchos::RCP<std::vector<std::vector<int> > >        	vec2D_int_ptr;
    typedef Teuchos::RCP<vec2D_longlong >				        	vec2D_longlong_ptr;
    typedef Teuchos::RCP<Epetra_Vector> 							EpetraVec_ptr;
    typedef Teuchos::RCP<Epetra_MpiComm>		 					EpetraComm_ptr;
    typedef Teuchos::RCP<Epetra_IntVector>	 						EpetraVecInt_ptr;
    typedef Teuchos::RCP<Epetra_LongLongVector>	 					EpetraVecLongLong_ptr;
    typedef Teuchos::RCP<Epetra_MultiVector>	 					EpetraMVPtr_Type;
    typedef Teuchos::RCP<Epetra_Map>                               	EpetraMapPtr_Type;

    
    typedef EpetraExt::HDF5 HDF5_Type;
    typedef Teuchos::RCP<HDF5_Type> HDF5Ptr_Type;
    
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;
    typedef const Teuchos::RCP<const Comm_Type> CommConstPtrConst_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef const MapConstPtr_Type MapConstPtrConst_Type;
    
    typedef MultiVector<SC,LO,GO,NO> MultiVec_Type;
    typedef Teuchos::RCP<const MultiVec_Type> MultiVecConstPtr_Type;
    typedef const MultiVecConstPtr_Type MultiVecConstPtrConst_Type;
    
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<Mesh_Type> MeshPtr_Type;
    
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;
    
    ExporterParaView();

    void setup(std::string filename,
               MeshPtr_Type mesh,
               std::string FEType,
               ParameterListPtr_Type parameterList=Teuchos::null);
    
    void setup(std::string filename,
               MeshPtr_Type mesh,
               std::string FEType,               
               int saveTimestep,
               ParameterListPtr_Type parameterList=Teuchos::null);

//    void setup(int dim,
//               GO nmbElementsGlob,
//               vec2D_int_ptr elements,
//               vec2D_dbl_ptr pointsUni,
//               MapConstPtrConst_Type& mapUnique,
//               MapConstPtrConst_Type& mapRepeated,
//               std::string FEType,
//               std::string filename,
//               CommConstPtrConst_Type &comm,
//               ParameterListPtr_Type parameterList=Teuchos::null);
//    
//    void setup(int dim,
//               GO nmbElementsGlob,
//               vec2D_int_ptr elements,
//               vec2D_dbl_ptr pointsUni,
//               MapConstPtrConst_Type& mapUnique,
//               MapConstPtrConst_Type& mapRepeated,
//               std::string FEType,
//               std::string filename,
//               int saveTimestep,
//               CommConstPtrConst_Type &comm,
//               ParameterListPtr_Type parameterList=Teuchos::null);

    
    void addVariable(MultiVecConstPtr_Type &u,
                     std::string varName,
                     std::string varType,
                     int dofPerNode,
                     MapConstPtrConst_Type& mapUnique=Teuchos::null);
    
    void save(double time);
    
    void save(double time, double dt);
    
    void closeExporter();
    
    void writeVariablesHDF5();
    
    void initHDF5();
    
    void writeMeshElements( std::string nameConn );
    
    void writeMeshPointsHDF5();
    
    void writeMeshPoints(std::string nameP_X,
                    std::string nameP_Y,
                    std::string nameP_Z );
    
    void updatePoints();
    
    
    void initXmf();
    
    void initXmfTimes();
    
    void writeXmf(double time);
    
    void writeXmfElements(std::string nameConn, double time);
    
    void writeXmfPoints(std::string nameP_X,
                             std::string nameP_Y,
                             std::string nameP_Z);
    
    void writeXmfVariables();
    
    void writeXmfTime(double time, double dt);
    
    void prepareVectorField(MultiVecConstPtr_Type &u,
                            EpetraMVPtr_Type &u_export,
                            int dof) const;
    
    void prepareScalar(MultiVecConstPtr_Type &u,
                       EpetraMVPtr_Type &u_export) const;
    
    void makePostfix();
    
protected:
    
    HDF5Ptr_Type hdf5exporter_;
    CommConstPtr_Type comm_;
    Teuchos::RCP<Epetra_MpiComm> commEpetra_;
    
    std::streampos				closingLinesPosition_;
    std::streampos              closingLinesPositionTimes_;
    std::string          		closingLines_;
    std::ofstream 				xmf_out_;
    std::ofstream 				xmf_times_out_;
    std::string 				filename_;
    std::string			 		outputFilename_;
    std::string                 postfix_;
    std::string                 FEType_;
    
    
    std::vector<MultiVecConstPtr_Type> variables_;
    std::vector<EpetraMapPtr_Type >   uniqueMaps_;
    std::vector<std::string>   		varNames_;
    std::vector<std::string>   		varTypes_;
    std::vector<int> 				varDofPerNode_;
    EpetraMVPtr_Type				pointsHDF_;
    EpetraVecInt_ptr				elementsHDF_;
    
    UN     			dim_;
    GO             	nmbElementsGlob_;
    GO 				nmbPointsGlob_;
    GO			nmbExportValuesGlob_;
    int 				nmbPointsPerElement_;
    int 				timeIndex_;
    bool				writeDt_;
    int 				saveTimestep_;
    bool verbose_;
    ParameterListPtr_Type parameterList_;
    vec2D_dbl_ptr pointsUnique_;

	bool redo_ = false;
	MeshPtr_Type mesh_;
    MapConstPtr_Type mapUniqueVariables_;
    
    };
}

#endif
