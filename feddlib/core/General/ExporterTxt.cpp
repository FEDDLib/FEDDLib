#include "ExporterTxt.hpp"
/*!
 Definition of ExporterTxt
 
 @brief  ExporterTxt
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
ExporterTxt::ExporterTxt():
    verbose_(false),
    txt_out_()
{
    
}


void ExporterTxt::setup(std::string filename, CommConstPtr_Type comm, int targetRank){
    
    if (comm->getRank() == targetRank)
        verbose_ = true;

    if (verbose_)
        txt_out_.open((filename + ".txt").c_str(),std::ios_base::out);
}
    
void ExporterTxt::exportData(double data){
   
    writeTxt(data);
    
}

void ExporterTxt::writeTxt(double data){
    if (verbose_) {
        txt_out_ << data << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::closeExporter(){
    if (verbose_) {
        txt_out_.close();
    }

}
}


