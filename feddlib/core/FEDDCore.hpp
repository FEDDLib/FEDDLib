#ifndef FEDDCORE_hpp
#define FEDDCORE_hpp
#define UNDERLYING_LIB_TPETRA

#define FEDD_TIMER
#define FEDD_DETAIL_TIMER

#include <stdint.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iosfwd>
#include <string>
#include <limits>
#include <chrono> 

#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"

#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPBoostSharedPtrConversions.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <boost/function.hpp>
namespace FEDD {
//    template<class SC, class LO, class GO, class NO>
//    class InterfaceEntity;
//    enum EntityType {DefaultType,VertexType,EdgeType,FaceType,InteriorType,InterfaceType};
//    enum EntityFlag {DefaultFlag,StraightFlag,ShortFlag,NodeFlag};
//    enum DistanceFunction {ConstantDistanceFunction,InverseEuclideanDistanceFunction};
//
//    template <class SC,class LO,class GO,class NO>
//    bool compareInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
//                                  Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb)
//    {
//        return iEa->getUniqueID()<iEb->getUniqueID();
//    }
//
//    template <class SC,class LO,class GO,class NO>
//    bool equalInterfaceEntities(Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEa,
//                                Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > iEb)
//    {
//        return iEa->getUniqueID()==iEb->getUniqueID();
//    }

//    
//#ifndef FEDD_TIMER_START
//#define FEDD_TIMER_START(A,S) Teuchos::RCP<TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("FEDD:") + std::string(S))));
//#endif
//    
//#ifndef FEDD_TIMER_STOP
//#define FEDD_TIMER_STOP(A) A.reset();
//#endif

#define FEDDLIB_WARNING(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cerr << std::setw(5) << " " << "[WARNING] " << CLASS << ": " << OUTPUT << std::endl;

#define FEDDLIB_NOTIFICATION(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cerr << std::setw(5) << " " << "[NOTIFICATION] " << CLASS << ": " << OUTPUT << std::endl;

    
typedef unsigned UN;

typedef std::tuple<int,int> tuple_intint_Type;
typedef std::tuple<std::string,std::string,int,int> tuple_ssii_Type;
typedef std::tuple<std::string,double> tuple_sd_Type;
typedef std::vector<tuple_sd_Type> tuple_sd_vec_Type;
typedef Teuchos::RCP<tuple_sd_vec_Type> tuple_sd_vec_ptr_Type;

typedef std::vector<tuple_ssii_Type> tuple_disk_vec_Type;
typedef Teuchos::RCP<tuple_disk_vec_Type> tuple_disk_vec_ptr_Type;

typedef std::vector<std::string>                                    string_vec_Type;
typedef Teuchos::RCP<string_vec_Type>								string_vec_ptr_Type;

typedef std::vector<std::string>                                    vec_string_Type;
typedef std::vector<double>											vec_dbl_Type;
typedef std::vector<int>											vec_int_Type;
typedef std::vector<long long>										vec_long_Type;
typedef std::vector<default_lo>										vec_LO_Type;
typedef std::vector<default_go>										vec_GO_Type;
typedef std::vector<bool>											vec_bool_Type;
typedef std::vector<std::vector<double> >							vec2D_dbl_Type;
typedef std::vector<std::vector<int> >                              vec2D_int_Type;
typedef std::vector<vec_long_Type >                                 vec2D_long_Type;
typedef std::vector<vec_LO_Type >                                   vec2D_LO_Type;
typedef std::vector<vec_GO_Type >                                   vec2D_GO_Type;
typedef std::vector<vec2D_dbl_Type>                                 vec3D_dbl_Type;
typedef std::vector<vec2D_long_Type >                               vec3D_long_Type;

typedef Teuchos::RCP<vec_dbl_Type>									vec_dbl_ptr_Type;
typedef Teuchos::RCP<vec_int_Type>									vec_int_ptr_Type;
typedef Teuchos::RCP<vec_long_Type>									vec_long_ptr_Type;
typedef Teuchos::RCP<vec_GO_Type>									vec_GO_ptr_Type;
typedef Teuchos::RCP<vec2D_dbl_Type >				     			vec2D_dbl_ptr_Type;
typedef Teuchos::RCP<vec2D_int_Type >				     			vec2D_int_ptr_Type;
typedef Teuchos::RCP<vec2D_long_Type >				     			vec2D_long_ptr_Type;
typedef Teuchos::RCP<vec3D_dbl_Type>				     			vec3D_dbl_ptr_Type;
typedef std::vector<vec_long_ptr_Type>                              vec_long_ptr_vec_Type;
typedef Teuchos::RCP<vec_long_ptr_vec_Type>                         vec_long_ptr_vec_ptr_Type;
typedef std::vector<vec2D_int_Type>                                 vec2D_int_vec_Type;
typedef std::vector<vec2D_dbl_Type>                                 vec2D_dbl_vec_Type;
typedef std::vector<vec2D_int_ptr_Type>                             vec2D_int_ptr_vec_Type;
typedef std::vector<vec2D_dbl_ptr_Type>                             vec2D_dbl_ptr_vec_Type;

typedef Teuchos::RCP<vec2D_int_vec_Type>							vec2D_int_vec_ptr_Type;
typedef Teuchos::RCP<vec2D_dbl_vec_Type>							vec2D_dbl_vec_ptr_Type;
typedef Teuchos::RCP<vec2D_int_ptr_vec_Type>						vec2D_int_ptr_vec_ptr_Type;
typedef Teuchos::RCP<vec2D_dbl_ptr_vec_Type>						vec2D_dbl_ptr_vec_ptr_Type;

typedef Teuchos::RCP<vec3D_long_Type>				     			vec3D_long_ptr_Type;

typedef Teuchos::RCP<std::vector<vec_int_ptr_Type> >				vec_int_ptr_vec_ptr_Type;

typedef Teuchos::Time                                               Time_Type;
typedef Teuchos::RCP<Time_Type>                                     TimePtr_Type;
typedef Teuchos::TimeMonitor                                        TimeMonitor_Type;

typedef Teuchos::ParameterList                                      ParameterList_Type;
typedef Teuchos::RCP<ParameterList_Type>                            ParameterListPtr_Type;
typedef Teuchos::RCP<const ParameterList_Type>                      ParameterListConstPtr_Type;

typedef boost::function<double(double* x, int* parameters)>                             CoeffFunc_Type;
typedef boost::function<double(double* x, double* parameters)>                          CoeffFuncDbl_Type;
typedef boost::function<void(double* x, double* res, double* parameters)>               RhsFunc_Type;
typedef boost::function<void(double* x, double* res, double t, double* parameters)>     GeneralFunc_Type;        
typedef boost::function<void(double* x, double* res)>     								Func_Type;          
    
}
#endif
