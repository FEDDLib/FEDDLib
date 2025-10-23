#ifndef ExporterTxt_hpp
#define ExporterTxt_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/core_config.h"
#include <fstream>

/*!
 Declaration of ExporterTxt
 
 @brief  ExporterTxt
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
class ExporterTxt {
public:
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;
    
    bool verbose_;

    std::ofstream txt_out_;
    
    ExporterTxt();

    void setup(std::string filename, CommConstPtr_Type comm, int targetRank=0);
    
    void exportData(double data);

    void writeTxt(double data);

    void exportData(double data1, double data2);

    void exportData(std::string data1, double data2);

    void exportData(std::string data1, std::string data2);

    void writeTxt(double data1, double data2);
    
    void writeTxt(std::string data1, double data2);

    void writeTxt(std::string data1, std::string data2);

    void closeExporter();

    
private:
    };
}

#endif
