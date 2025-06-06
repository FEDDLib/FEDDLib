#include "TimeSteppingTools.hpp"
/*!
 Definition of TimeSteppingTools
 
 @brief  TimeSteppingTools
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
TimeSteppingTools::TimeSteppingTools():
comm_(),
parameterList_(),
butcherTableNmb_(0),
BDFNmb_(0),
tsType_(),
tEnd_(0.1),
dt_(0.1),
t_(0.0),
/* adaptive variables*/
dt_prev_(0.),
dt_adaptive_(0.),
rho_(1),
tolAdaptive_(0.01),
dtmin_(0.001),
dtmax_(1.),
convOrder_(1),
adaptiveError_(),
adaptiveCalculation_(),
error_(0.),
error_prev_(0.),
stages_(1),
stifflyAcc_(false),
stifflyAccEmbedded_(false),
butcherTable_(),
BDFInformation_(),
b_embedded_(),
gamma_vec_(),
verbose_(false),
exporterTxtTime_(),
exporterTxtDt_(),
exporterTxtError_(),
beta_(0.25),
gamma_(0.5)
{

}

TimeSteppingTools::TimeSteppingTools(ParameterListPtr_Type parameterList, CommConstPtr_Type comm):
comm_(comm),
parameterList_(parameterList),
butcherTableNmb_(0),
BDFNmb_(0),
tsType_(),
tEnd_(0.1),
dt_(0.1),
t_(0.0),
/* adaptive variables*/
dt_prev_(0.1),
dt_adaptive_(0.),
rho_(1),
tolAdaptive_(0.01),
dtmin_(0.001),
dtmax_(1.),
convOrder_(1),
adaptiveError_(),
adaptiveCalculation_(),
error_(0.),
error_prev_(0.),
stages_(1),
stifflyAcc_(false),
stifflyAccEmbedded_(false),
butcherTable_(),
BDFInformation_(),
b_embedded_(),
gamma_vec_(),
verbose_(comm->getRank() == 0),
exporterTxtTime_(),
exporterTxtDt_(),
exporterTxtError_(),
beta_(0.25),
gamma_(0.5)
{
    setParameter();
}

void TimeSteppingTools::setParameter(){
    
    tEnd_ 	= parameterList_->get("Final time",1.);
    dt_		= parameterList_->get("dt",0.01);
    dt_prev_ = dt_;

    beta_ = parameterList_->get("beta",0.25);
    gamma_ = parameterList_->get("gamma",0.5);

    if (!parameterList_->get("Class","Singlestep").compare("Singlestep")) {
        butcherTableNmb_ = parameterList_->get("Butcher table",0);
        if ( !parameterList_->get("Timestepping type","non-adaptive").compare("adaptive") ) {
            tsType_ = ADAPTIVE;
            setupTxtExporter();
        }
        else{
            tsType_ = NON_ADAPTIVE;
        }
        /* adaptive variables*/
        rho_	= parameterList_->get("Safety factor",1.);
        tolAdaptive_ = parameterList_->get("Adaptive tolerance",0.01);
        dtmin_	= parameterList_->get("dtmin",0.001);
        dtmax_	= parameterList_->get("dtmax",0.1);
        if (parameterList_->get("Adaptive error",0)==0) {
            adaptiveError_ = EUCLIDIAN;
        }
        else if(parameterList_->get("Adaptive error",0)==1){
            adaptiveError_ = L2;
        }

        if (parameterList_->get("Adaptive calculation",0)==0) {
            adaptiveCalculation_ = 0;
        }
        else if(parameterList_->get("Adaptive calculation",0)==1){
            adaptiveCalculation_ = 1;
        }


        setTableInformationRK();

    }
    else if(!parameterList_->get("Class","Singlestep").compare("Multistep")){
        tsType_ = NON_ADAPTIVE;
        BDFNmb_ = parameterList_->get("BDF",1);
        setInformationBDF();
    }
}

double TimeSteppingTools::currentTime(){
    return t_;
}

bool TimeSteppingTools::continueTimeStepping(){
    if (t_+1.e-10<tEnd_)
        return true;
    else
        return false;
}

double TimeSteppingTools::getButcherTableCoefficient(int row , int col){

    return butcherTable_->at(row).at(col+1);

}

double TimeSteppingTools::getButcherTableC(int row){

    return butcherTable_->at(row).at(0);
}

double TimeSteppingTools::get_dt(){

    return dt_;

}

double TimeSteppingTools::get_dt_prev(){

    return dt_prev_;

}

void TimeSteppingTools::calculateSolution(BlockMultiVectorPtr_Type &sol, BlockMultiVectorPtrArray_Type &rkSolVec, BlockMatrixPtr_Type massSystem, BlockMultiVectorPtr_Type solShort){

    if (stifflyAcc_){

    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should not happen yet. RK method not stiffly accurate.");
    
    
    if (tsType_ == ADAPTIVE)
        adaptiveTimestep(sol , rkSolVec, massSystem, solShort);
}

void TimeSteppingTools::adaptiveTimestep(BlockMultiVectorPtr_Type &sol, BlockMultiVectorPtrArray_Type &rkSolVec, BlockMatrixPtr_Type massSystem, BlockMultiVectorPtr_Type solShortIN){


    BlockMultiVectorPtr_Type solShort = Teuchos::rcp(new BlockMultiVector_Type ( sol ) ); //We do not need a copy here!

    if (stifflyAccEmbedded_) {
        solShort = rkSolVec.at(stages_-2);
    }
    else{
        if (gamma_vec_->at(0)!=-1.) {
            solShort->putScalar(0.);
            for (int i=0; i<gamma_vec_->size(); i++) {
                if (gamma_vec_->at(i)!=0.0) {
                    BlockMultiVectorPtr_Type RU = Teuchos::rcp(new BlockMultiVector_Type( rkSolVec[i] ) );
                    RU->update(-1.,*sol,1.);//
                    solShort->update( gamma_vec_->at(i), *RU, 1. );
                }
            }
        }
    }
    if(!solShortIN.is_null())
        solShortIN->update( 1., *solShort, 0. );
        
    solShort->update( -1., *sol, 1. );

    calculateNewDt(solShort, massSystem);

}

void TimeSteppingTools::calculateNewDt(BlockMultiVectorPtr_Type &solDiff, BlockMatrixPtr_Type massSystem){

    TEUCHOS_TEST_FOR_EXCEPTION(solDiff->getNumVectors()>1, std::logic_error, "Implement for MV.size()>1.");
    error_prev_ = error_;
    Teuchos::Array<SC> errorTmp(1);
    
    if ( adaptiveError_ == EUCLIDIAN ) {
        solDiff->norm2(errorTmp());
        error_ = errorTmp[0];
    }

    else if ( adaptiveError_ == L2 ) {
        BlockMultiVectorPtr_Type result = Teuchos::rcp(new BlockMultiVector_Type ( solDiff) ); //We do not need a copy here!
        massSystem->apply( *solDiff, *result );
        solDiff->dot( result, errorTmp );
        error_ = std::sqrt( errorTmp[0] );
    }

    if ( adaptiveCalculation_ == 0 ) {
        if (t_==0) {
            error_prev_ = error_;
        }
        dt_adaptive_ = rho_ * dt_*dt_ / dt_prev_ * std::pow(tolAdaptive_*(error_prev_)/((error_)*(error_)),1./convOrder_);
    }
    else if ( adaptiveCalculation_ == 1 ) {
        dt_adaptive_ = rho_ * dt_ * std::pow(tolAdaptive_/(error_),1./convOrder_);
    }

    dt_prev_ = dt_;
    if(dt_adaptive_ > dtmax_){
        dt_ = dtmax_;
    }
    else if(dt_adaptive_ < dtmin_){
        dt_ = dtmin_;
    }
    else{
        dt_ = dt_adaptive_;
    }
}

void TimeSteppingTools::correctPressure(MultiVectorPtr_Type &newP/*should be the lastest (false) pressure solution*/, MultiVectorConstPtr_Type lastP){

    TEUCHOS_TEST_FOR_EXCEPTION(butcherTableNmb_!=1, std::logic_error, "Only for CN.");
    
    if (butcherTableNmb_ == 1) {
        newP->scale( 1./butcherTable_->at(1).at(2) );
        newP->update( -butcherTable_->at(1).at(1)/butcherTable_->at(1).at(2), *lastP, 1.);
    }
}


void TimeSteppingTools::advanceTime(bool printInfo){

    if (tsType_ == ADAPTIVE) {
        exporterTxtTime_->exportData(t_);
        exporterTxtDt_->exportData(dt_prev_);
        exporterTxtError_->exportData(t_);
    }

    t_+=dt_;

    if (printInfo)
        this->printInfo();

}

void TimeSteppingTools::printInfo(){

    if (verbose_) {
        std::cout <<  std::endl;
        std::cout << " --------------------------------------------------------------" << std::endl;
        std::cout << " -------- Time step " << currentTime() << std::endl;
        std::cout << "        - dt     = " << get_dt() << "\t calculated dt = " << dt_adaptive_ << std::endl;
        std::cout << "        - r_m = " << error_ << "  r_m_prev_ = " << error_prev_ << std::endl;
//        std::cout << "        - rho*tm^2/tm_1 = " <<rho * dt*dt / dt_prev_ << std::endl;
//        std::cout << "        - TOL*rm-1 / (rm^2) = "<<TOLRK*(r_m_prev_)/((r_m_)*(r_m_)) << "  (...)^1/p = " << std::pow(TOLRK*(r_m_prev_)/((r_m_)*(r_m_)),1./p) << std::endl;
        std::cout << " --------------------------------------------------------------" << std::endl;
        std::cout <<  std::endl;
        }
}

int TimeSteppingTools::getNmbStages(){
    return stages_;
}


void TimeSteppingTools::setupTxtExporter(){

    exporterTxtTime_.reset(new ExporterTxt());
    exporterTxtDt_.reset(new ExporterTxt());
    exporterTxtError_.reset(new ExporterTxt());

    exporterTxtTime_->setup("time",comm_);
    exporterTxtDt_->setup("dt",comm_);
    exporterTxtError_->setup("error",comm_);
}

void TimeSteppingTools::setTableInformationRK(){

    // Butcher tables :
    //                      c | A
    //                      --------
    //                        | b^T
    double Theta = 1. - std::sqrt(2.) / 2.;
    double Theta_t = 1. - 2.*Theta;
    double alpha = Theta_t / (1. - Theta);
    double beta = 1. - alpha;
    double gamma;
    double sigma;

    b_embedded_.reset(new vec_dbl_Type(1,-1.));
    gamma_vec_.reset(new vec_dbl_Type(1,-1.));
    RKType_ = butcherTableNmb_;
    switch (butcherTableNmb_) {
        case 0: //implicit;
            
            butcherTable_.reset(new vec2D_dbl_Type(3,vec_dbl_Type(3,0.)));

            butcherTable_->at(1).at(1) = 0.; //a21
            butcherTable_->at(1).at(2) = 1.;

            butcherTable_->at(1).at(0) = 1.;

            stifflyAcc_ = true;

            stifflyAccEmbedded_ = false;
            break;
        case 1: //CN
            
            butcherTable_.reset(new vec2D_dbl_Type(3,vec_dbl_Type(3,0.)));

            butcherTable_->at(1).at(1) = .5; //a21
            butcherTable_->at(1).at(2) = .5;

            butcherTable_->at(0).at(0) = 0.; //c1
            butcherTable_->at(1).at(0) = 1.; //c2

            stifflyAcc_ = true;

            stifflyAccEmbedded_ = false;
            break;
        case 2: //pressure-corrected fractional-setp theta-scheme
            
            butcherTable_.reset(new vec2D_dbl_Type(5,vec_dbl_Type(5)));
            butcherTable_->at(0).at(0) = 0.; //c1
            butcherTable_->at(0).at(1) = 0.; //a11
            butcherTable_->at(0).at(2) = 0.;
            butcherTable_->at(0).at(3) = 0.; //a13
            butcherTable_->at(0).at(4) = 0.; //a14
            butcherTable_->at(1).at(0) = Theta; //c2
            butcherTable_->at(1).at(1) = Theta * beta; //a21
            butcherTable_->at(1).at(2) = Theta * alpha;
            butcherTable_->at(1).at(3) = 0.;
            butcherTable_->at(1).at(4) = 0.;
            butcherTable_->at(2).at(0) = Theta + Theta_t;
            butcherTable_->at(2).at(1) = Theta*beta;
            butcherTable_->at(2).at(2) = (Theta + Theta_t) * alpha;
            butcherTable_->at(2).at(3) = Theta_t * alpha;
            butcherTable_->at(2).at(4) = 0.;
            butcherTable_->at(3).at(0) = 1.;
            butcherTable_->at(3).at(1) = Theta * beta;
            butcherTable_->at(3).at(2) = (Theta + Theta_t) * alpha;
            butcherTable_->at(3).at(3) = (Theta + Theta_t) * beta;
            butcherTable_->at(3).at(4) = Theta * alpha;
            butcherTable_->at(4).at(0) = 0.;
            butcherTable_->at(4).at(1) = Theta * beta; //b1
            butcherTable_->at(4).at(2) = (Theta + Theta_t) * alpha;
            butcherTable_->at(4).at(3) = (Theta + Theta_t) * beta;
            butcherTable_->at(4).at(4) = Theta * alpha;

            stifflyAcc_ = true;//true;
            // do we need this?
            b_embedded_.reset(new vec_dbl_Type(4,0.));
            b_embedded_->at(0)=0.11785113033497070959;
            b_embedded_->at(1)=0.49509379160690495120;
            b_embedded_->at(2)=0.29636243203812433921;
            b_embedded_->at(3)=0.09069264621404818692;

            gamma_vec_.reset(new vec_dbl_Type(4,0.));

            gamma_vec_->at(1) = -0.382148867894378;
            gamma_vec_->at(2) = 0.824957911172935;
            gamma_vec_->at(3) = 0.528595479020469;

            stifflyAccEmbedded_ = false;
            convOrder_ = 2;
            break;
        case 3: //DIRK3L; p=2
            
            alpha = 1.- std::sqrt(2.)/2.;
            butcherTable_.reset(new vec2D_dbl_Type(4,vec_dbl_Type(4,0.)));

            butcherTable_->at(1).at(0) = 2.*alpha;

            butcherTable_->at(1).at(1) = alpha;
            butcherTable_->at(1).at(2) = alpha;

            butcherTable_->at(2).at(0) = 1.;

            butcherTable_->at(2).at(1) = 1.- (1. / (4.*(1.-alpha))) - alpha;
            butcherTable_->at(2).at(2) = 1. / (4.*(1.-alpha));
            butcherTable_->at(2).at(3) = alpha;

            gamma_vec_.reset(new vec_dbl_Type(3,0.));

            gamma_vec_->at(1) = -0.353553390593274;
            gamma_vec_->at(2) = 1.207106781186547;

            stifflyAccEmbedded_ = false;
            stifflyAcc_ = true;
            convOrder_ = 2;
            break;
        case 4: //DIRK34; p=3
            
            alpha = 0.1558983899988677;
            beta = 1.072486270734370;
            gamma = 0.7685298292769537;
            sigma = 0.09666483609791597;
            butcherTable_.reset(new vec2D_dbl_Type(5,vec_dbl_Type(5,0.)));


            butcherTable_->at(1).at(0) = 2.*alpha; //c2
            butcherTable_->at(2).at(0) = 1.; //c3
            butcherTable_->at(3).at(0) = 1.; //c4

            butcherTable_->at(1).at(1) = alpha; //a21
            butcherTable_->at(1).at(2) = alpha;
            butcherTable_->at(2).at(3) = alpha;
            butcherTable_->at(3).at(4) = alpha;

            butcherTable_->at(2).at(2) = beta;

            butcherTable_->at(2).at(1) = 1.-butcherTable_->at(2).at(2)-butcherTable_->at(1).at(2);

            butcherTable_->at(3).at(2) = gamma;
            butcherTable_->at(3).at(3) = sigma;

            butcherTable_->at(3).at(1) = 1.-gamma-sigma-alpha;

            stifflyAcc_ = true;
            stifflyAccEmbedded_ = true;
            convOrder_ = 3;
            break;

//        case 5: //DIRK3; p=3
//            alpha = 0.5 + std::sqrt(3)/6.;
//            butcherTable_.reset(new vec2D_dbl_Type(4,vec_dbl_Type(4,0.)));
//
//            butcherTable_->at(1).at(0) = 2.*alpha;
//            butcherTable_->at(2).at(0) = 1.;
//
//
//            butcherTable_->at(1).at(1) = alpha;
//            butcherTable_->at(1).at(2) = alpha;
//
//            butcherTable_->at(2).at(2) = 1. / (12.*alpha*(2.*alpha-1.));
//            butcherTable_->at(2).at(1) = 1.- butcherTable_->at(2).at(3) - butcherTable_->at(2).at(2);
//
//
//            butcherTable_->at(2).at(3) = alpha;
//
//            gamma_vec_.reset(new vec_dbl_Type(3,0.));
//
//            gamma_vec_->at(1) = 0.;
//            gamma_vec_->at(2) = 0.;
//
//            stifflyAccEmbedded_ = false;
//            stifflyAcc_ = true;
//            convOrder_ = 3;
//            break;

        default:
            break;
    }

    stages_ = butcherTable_->size()-1;
}
double TimeSteppingTools::getInformationBDF(int i){
    TEUCHOS_TEST_FOR_EXCEPTION(i+1>BDFInformation_->size() || i<0, std::logic_error, "Wrong bdf table access!");
    return BDFInformation_->at(i);
}

void TimeSteppingTools::setInformationBDF(){
    BDFInformation_.reset(new vec_dbl_Type(BDFNmb_+2));
    switch (BDFNmb_) {
        case 1: //implicit;
            BDFInformation_->at(0) = 1.; //Mun+1/dt
            BDFInformation_->at(1) = 1.; //Aun+1
            BDFInformation_->at(2) = 1.; //Mun/dt
            break;
        case 2:
            // BDFInformation_->at(0) = 3.; //Mun+2/dt
            // BDFInformation_->at(1) = 2.; //Aun+2
            // BDFInformation_->at(2) = 4.; //Mun+1/dt
            // BDFInformation_->at(3) = -1.; //Mun
            BDFInformation_->at(0) = 1.5; //Mun+2/dt
            BDFInformation_->at(1) = 1.0; //Aun+2
            BDFInformation_->at(2) = 2.0; //Mun+1/dt
            BDFInformation_->at(3) = -0.5; //Mun/dt
            break;
        default:
            break;
    }
}

int TimeSteppingTools::getBDFNumber(){

    return BDFNmb_;
}

double TimeSteppingTools::get_beta()
{
    return beta_;
}

double TimeSteppingTools::get_gamma()
{
    return gamma_;
}

}
