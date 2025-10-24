#include "ExporterTxt.hpp"
/*!
 Definition of ExporterTxt
 
 @brief  ExporterTxt
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
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
        txt_out_.open((filename + ".txt").c_str(),ios_base::out);
}
    
void ExporterTxt::exportData(double data){
   
    writeTxt(data);
    
}
void ExporterTxt::exportData(std::string data1, double data2){
   
    writeTxt(data1,data2);
    
}

void ExporterTxt::exportData(std::string data1, std::string data2){
   
    writeTxt(data1,data2);
    
}
void ExporterTxt::exportData(double data1, double data2){
   
    writeTxt(data1,data2);
    
}

void ExporterTxt::writeTxt(double data){
    if (verbose_) {
        txt_out_ << data << "\n";
        txt_out_.flush();
    }

}


void ExporterTxt::writeTxt(double data1, double data2){
    if (verbose_) {
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::writeTxt(std::string data1, double data2){
    if (verbose_) {
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::writeTxt(std::string data1, std::string data2){
    if (verbose_) {
        txt_out_ << data1 << " " << data2 << "\n";
        txt_out_.flush();
    }

}

void ExporterTxt::closeExporter(){
    if (verbose_) {
        txt_out_.close();
    }

}
}


