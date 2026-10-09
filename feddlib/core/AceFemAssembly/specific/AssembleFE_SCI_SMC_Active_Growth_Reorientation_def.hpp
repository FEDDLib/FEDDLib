#ifndef AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
#define AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp

#include <string>
#include <set>
#include <vector>

#include "AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::AssembleFE_SCI_SMC_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple),
                                                                                                                                                                                                                       segmentsActive_(0),
                                                                                                                                                                                                                       segmentsGrowth_(0),
                                                                                                                                                                                                                       segmentsReorientation_(0) {
#ifndef FEDD_HAVE_ACEGENINTERFACE
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "AssembleFE_SCI_SMC_Active_Growth_Reorientation needs FEDDLib built with the AceGen interface (Interface2): configure with -D TPL_ENABLE_AceGENInterface=ON.");
#endif
    activeAcceleratedEndTime_ = this->params_->sublist("Parameter").get("Accelerated Active Until", -1.);
    TEUCHOS_TEST_FOR_EXCEPTION(activeAcceleratedEndTime_ < 0, std::logic_error, "!!! Warning: Accelerated Active Until not set correctly. Please Check Parameterlist !!!");
										
    activeAcceleratedMultiplier_ = this->params_->sublist("Parameter").get("Active Acceleration Multiplier", -1.);
    TEUCHOS_TEST_FOR_EXCEPTION(activeAcceleratedMultiplier_ < 0, std::logic_error, "!!! Warning: Active Acceleration Multiplier not set correctly. Please Check Parameterlist !!!");

    int numMaterials = this->params_->sublist("Parameter Solid").get("Number of Materials", 1);
    int materialID = 0;
    for (int i = 1; i <= numMaterials; i++)
        if (this->params_->sublist("Parameter Solid").sublist(std::to_string(i)).get("Volume Flag", 15) == this->flag_)
            materialID = i;
	
    TEUCHOS_TEST_FOR_EXCEPTION(materialID == 0, std::logic_error, "!!! Warning: No corresponding parameterslist for the element flag = " << this->flag_ << ". Please Check volume flags of elements and Mesh Data !!!");

    this->iCode_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Integration Code", 18);

#ifdef FEDD_HAVE_ACEGENINTERFACE

    AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 tempElem = AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10(this->iCode_);
    this->historyLength_ = tempElem.getHistoryLength();
    this->numberOfIntegrationPoints_ = tempElem.getNumberOfGaussPoints();
    this->postDataLength_ = tempElem.getNumberOfPostData();
    this->domainDataLength_ = tempElem.getNumberOfDomainData();
    char** domainDataNames = tempElem.getDomainDataNames();
    char** postDataNames = tempElem.getPostDataNames();

    std::vector<std::string> subString(this->domainDataLength_);

    this->domainDataNames_.resize(this->domainDataLength_);
    this->postDataNames_.resize(this->postDataLength_);
    this->domainData_.resize(this->domainDataLength_, 0.0);

    for (int i = 0; i < this->domainDataLength_; i++) {
        this->domainDataNames_[i] = std::string(domainDataNames[i]);
        int pos1 = domainDataNames_[i].find("-");
        int pos2 = domainDataNames_[i].find("_");

        subString[i] = domainDataNames_[i].substr(pos1 + 1, pos2 - pos1 - 1);
        this->domainDataNames_[i] = subString[i];
        this->domainData_[i] = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get(subString[i], 1.e13);

        TEUCHOS_TEST_FOR_EXCEPTION(this->domainData_[i] > 1.e12, std::logic_error, " Parameter not set correctly. Parameter " << this->domainDataNames_[i] << " received default value!!");

        // Pre-compute which parameters need acceleration/deceleration, by exact name: a
        // substring match would also scale kEtaPlus and mEtaPlus (through "Eta").
        static const std::set<std::string> acceleratedNames = {
            "LambdaBarCDotMax", "LambdaBarCDotMin", "Eta", "K3", "K4", "K7", "Beta1", "Gamma6",
            "KDotMin", "KDotMax", "LambdaBarDotPMin", "LambdaBarDotPMax"};
        static const std::set<std::string> deceleratedNames = {"Gamma5", "Gamma2"};
        if (acceleratedNames.count(domainDataNames_[i]) > 0)
            acceleratedParamIndices_.push_back(i);
        else if (deceleratedNames.count(domainDataNames_[i]) > 0)
            deceleratedParamIndices_.push_back(i);
    }

    for (int i = 0; i < this->postDataLength_; i++) {
        this->postDataNames_[i] = std::string(postDataNames[i]);
    }

    // Create map from field name to position in postDataNames_
    for (int i = 0; i < postDataNames_.size(); i++) {
        fieldNameToPosition_[postDataNames_[i]] = i;
    }

    this->residuumRint_.resize(30, 0.0);
    this->residuumRc_.resize(10, 0.0);
    this->residuumRdyn_.resize(30, 0.0);

    this->stiffnessMatrixKuu_.resize(30, vec_dbl_Type(30, 0.0));
    this->stiffnessMatrixKuc_.resize(30, vec_dbl_Type(10, 0.0));
    this->stiffnessMatrixKcu_.resize(10, vec_dbl_Type(30, 0.0));
    this->stiffnessMatrixKcc_.resize(10, vec_dbl_Type(10, 0.0));
    this->massMatrixMc_.resize(10, vec_dbl_Type(10, 0.0));

#endif

    this->subiterationTolerance_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Subiteration Tolerance", 1.e-7);
    this->FEType_ = std::get<1>(this->diskTuple_->at(0));     // FEType of Disk
    this->dofsSolid_ = std::get<2>(this->diskTuple_->at(0));  // Degrees of freedom per node
    this->dofsChem_ = std::get<2>(this->diskTuple_->at(1));   // Degrees of freedom per node

    this->numNodesSolid_ = std::get<3>(this->diskTuple_->at(0));  // Number of nodes of element
    this->numNodesChem_ = std::get<3>(this->diskTuple_->at(1));   // Number of nodes of element

    this->dofsElement_ = this->dofsSolid_ * this->numNodesSolid_ + this->dofsChem_ * this->numNodesChem_;  // "Dimension of return matrix"

    this->positions_ = std::vector<double>(30, 0.0);
    this->displacements_ = std::vector<double>(30, 0.0);
    this->accelerations_ = std::vector<double>(30, 0.0);
    this->concentrations_ = std::vector<double>(10, 0.0);
    this->rates_ = std::vector<double>(10, 0.0);

	// Einlesen durch Parameterdatei irgendwann cool

    // historyGP: Vector of history variables [Order: LambdaBarC1, LambdaBarC2, nA1, nA2, nB1, nB2, nC1, nC2, nD1, nD2, LambdaA1, LambdaA2, k251, k252, LambdaBarP1, LambdaBarP2, Theta1, Theta2, Theta3, Ag11, Ag12, Ag13, Ag21, Ag22, Ag23, Ag31, Ag32, Ag33, a11, a12, a13, a21, a22, a23] (The length must be equal to number of history variables per gauss point(34) * number of gauss points)
    std::vector<double> historyGP = {1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.82758, 1.82758, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    this->history_.clear();
    this->historyUpdated_.clear();

    this->history_.reserve(this->historyLength_);
    for (int i = 0; i < this->numberOfIntegrationPoints_; i++)
        this->history_.insert(this->history_.end(), historyGP.begin(), historyGP.end());

    // Error out if history length is inconsistent
    TEUCHOS_TEST_FOR_EXCEPTION(this->history_.size() != this->historyLength_, std::logic_error, "History input length does not match history size of model! \n History input length: " << this->history_.size() << "\n History size of model: " << this->historyLength_ << "\n");

    this->historyUpdated_.reserve(this->historyLength_);
    for (int i = 0; i < this->numberOfIntegrationPoints_; i++)
        this->historyUpdated_.insert(this->historyUpdated_.end(), historyGP.begin(), historyGP.end());

    this->solutionC_n_.resize(10, 0.);
    this->solutionC_n1_.resize(10, 0.);

    this->postProcessingData_ = Teuchos::rcp(new vec2D_dbl_Type(this->numNodesSolid_, vec_dbl_Type(this->postDataLength_)));
    this->solution_.reset(new vec_dbl_Type(this->dofsElement_, 0.));

#ifdef FEDD_HAVE_ACEGENINTERFACE	   
	
    // Nodal Positions in Reference Coordinates
    int count = 0;
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 3; j++) {
            this->positions_[count] = this->getNodesRefConfig()[i][j];
            count++;
        }

    // Active and Growth Time intervals
    int numSegmentsActive = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").get("Number of Segments", 0);
    int numSegmentsGrowth = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").get("Number of Segments", 0);
    int numSegmentsReorientation = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").get("Number of Segments", 0);

    for (int i = 1; i <= numSegmentsActive; i++) {
        double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").sublist(std::to_string(i)).get("Start Time", 0.);
        double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").sublist(std::to_string(i)).get("End Time", 0.);

        if (i == 1)
            TEUCHOS_TEST_FOR_EXCEPTION(startTime != this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("ActiveStartTime", 0.), std::logic_error, "!!!WARNING:: The ActiveStartTime and the start of Time stepping intervalls for active response do not match!!!");

        vec_dbl_Type segment = {startTime, endTime};
        segmentsActive_.push_back(segment);
    }

    for (int i = 1; i <= numSegmentsGrowth; i++) {
        double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").sublist(std::to_string(i)).get("Start Time", 0.);
        double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").sublist(std::to_string(i)).get("End Time", 0.);

        vec_dbl_Type segment = {startTime, endTime};
        segmentsGrowth_.push_back(segment);
    }

    for (int i = 1; i <= numSegmentsReorientation; i++) {
        double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").sublist(std::to_string(i)).get("Start Time", 0.);
        double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").sublist(std::to_string(i)).get("End Time", 0.);

        vec_dbl_Type segment = {startTime, endTime};
        segmentsReorientation_.push_back(segment);
    }
#endif
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleJacobian() {
    SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(this->dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

    assemble_SCI_SMC_Active_Growth_Reorientation(true);

    for (int i = 0; i < 30; i++)
        for (int j = 0; j < 30; j++)
            (*elementMatrix)[i][j] = this->stiffnessMatrixKuu_[i][j];

    for (int i = 0; i < 30; i++)
        for (int j = 0; j < 10; j++)
            (*elementMatrix)[i][j + 30] = this->stiffnessMatrixKuc_[i][j];

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 30; j++)
            (*elementMatrix)[i + 30][j] = this->stiffnessMatrixKcu_[i][j];

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            (*elementMatrix)[i + 30][j + 30] = this->stiffnessMatrixKcc_[i][j] + (1. / this->getTimeIncrement()) * this->massMatrixMc_[i][j];

#endif

    this->jacobian_ = elementMatrix;
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::checkingReorientationActiveGrowth() {
    bool restart = this->params_->sublist("Timestepping Parameter").get("Restart", false);
    double timeStepRestart = this->params_->sublist("Timestepping Parameter").get("Time step", 0.0);

    // Checking for active response
    for (int i = 0; i < segmentsActive_.size(); i++) {
        // if (time >= startTime && time < endTime)
        if ((this->timeStep_ > segmentsActive_[i][0] || approxEqual(this->timeStep_, segmentsActive_[i][0])) && (this->timeStep_ < segmentsActive_[i][1] && !approxEqual(this->timeStep_, segmentsActive_[i][1]))) {
            this->activeBool_ = 1;
            if (!this->activeInitialized_) {
                // if (restart and timeStepRestart > firstActiveStartTime)
                if (restart && timeStepRestart > segmentsActive_[0][0] && !approxEqual(timeStepRestart, segmentsActive_[0][0]))
                    this->activeInitialized_ = true;
                else
                    this->initializeActiveResponse();
            }
            break;
        } else
            this->activeBool_ = 0;
    }

    // Checking for growth
    for (int i = 0; i < segmentsGrowth_.size(); i++) {
        // if (time >= startTime && time < endTime)
        if ((this->timeStep_ > segmentsGrowth_[i][0] || approxEqual(this->timeStep_, segmentsGrowth_[i][0])) && (this->timeStep_ < segmentsGrowth_[i][1] && !approxEqual(this->timeStep_, segmentsGrowth_[i][1]))) {
            this->growthBool_ = 1;
            if (!this->growthInitialized_) {
                // if(restart and timeStepRestart > firstGrowthStartTime)
                if (restart && timeStepRestart > segmentsGrowth_[0][0] && !approxEqual(timeStepRestart, segmentsGrowth_[0][0]))
                    this->growthInitialized_ = true;
                else
                    this->initializeGrowth();
            }
            break;
        } else
            this->growthBool_ = 0;
    }

    for (int i = 0; i < segmentsReorientation_.size(); i++) {
        // if (time >= startTime && time < endTime)
        if ((this->timeStep_ > segmentsReorientation_[i][0] || approxEqual(this->timeStep_, segmentsReorientation_[i][0])) && (this->timeStep_ < segmentsReorientation_[i][1] && !approxEqual(this->timeStep_, segmentsReorientation_[i][1]))) {
            this->reorientationBool_ = 1;
            break;
        } else
            this->reorientationBool_ = 0;
    }

    if (activeBool_ == 1 && growthBool_ == 1)
        std::cout << " WARNING: Growth and Reorientation occuring simultaneously in element: " << this->globalElementID_ << ". Ignore message if intended.\n";

    if (activeBool_ == 1 && reorientationBool_ == 1)
        std::cout << " WARNING: Active Response and Reorientation occuring simultaneously in element: " << this->globalElementID_ << ". Ignore message if intended.\n";

    if (growthBool_ == 1 && reorientationBool_ == 1)
        std::cout << " WARNING: Growth and Reorientation occuring simultaneously in element: " << this->globalElementID_ << ". Ignore message if intended.\n";

    this->updateDomainData("growthBool", growthBool_);
    this->updateDomainData("ActiveBool", activeBool_);
    this->updateDomainData("reorientationBool", reorientationBool_);
}

// This is called at the beginning of each time step and sets time = t_n+1 along with the correct time increment dt
template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::advanceInTime(Teuchos::RCP<TimeSteppingTools> timeSteppingTool) {
    this->timeStep_ = timeSteppingTool->currentTime();  // Sets time to t_n+1

    this->timeIncrement_ = timeSteppingTool->get_dt();

    // Checking for Active Response Reorientation and Growth
    checkingReorientationActiveGrowth();

    for (int i = 0; i < this->historyLength_; i++)
        this->history_[i] = this->historyUpdated_[i];

    for (int i = 0; i < 10; i++)
        this->solutionC_n_[i] = (*this->solution_)[i + 30];
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleRHS() {
    this->rhsVec_.reset(new vec_dbl_Type(this->dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

    assemble_SCI_SMC_Active_Growth_Reorientation(false);

    for (int i = 0; i < 30; i++)
        (*this->rhsVec_)[i] = this->residuumRint_[i];  //+residuumRDyn[i];

    for (int i = 0; i < 10; i++)
        (*this->rhsVec_)[i + 30] = this->residuumRc_[i];
#endif
		}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assemble_SCI_SMC_Active_Growth_Reorientation(bool computeTangent) {
#ifdef FEDD_HAVE_ACEGENINTERFACE

    // We check this in case of restart. Then the history parameters are already set or we are already in a state of reorientation etc. Then the history does not need to be explicitly updated again.
    bool restart = this->params_->sublist("Timestepping Parameter").get("Restart", false);
    if (restart) {
        double timeStepRestart = this->params_->sublist("Timestepping Parameter").get("Time step", 0.0);
        if (approxEqual(this->timeStep_, timeStepRestart + this->getTimeIncrement()) && this->historyImported_)
            checkingReorientationActiveGrowth();
	}

    double deltaT = this->getTimeIncrement();

    double time = this->getTimeStep();

    for (int i = 0; i < 30; i++)
        this->displacements_[i] = (*this->solution_)[i];

    for (int i = 0; i < 10; i++) {
        this->concentrations_[i] = (*this->solution_)[i + 30];
        this->solutionC_n1_[i] = (*this->solution_)[i + 30];  // in each newtonstep solution for n+1 is updated.
    }

    for (int i = 0; i < 10; i++)
        this->rates_[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / deltaT;

    std::vector<double> domainDataModified(this->domainDataLength_);

    std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

    if (time < this->activeAcceleratedEndTime_) {
        for (int idx : acceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
        }
        for (int idx : deceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
        }
    }

    AceGenElement_Type& elem = aceGenElement(this->displacements_.data(), this->concentrations_.data(), this->accelerations_.data(), this->rates_.data(), domainDataModified.data(), deltaT, time);

    int errorCode = elem.compute(computeTangent);
    // With adaptive time stepping the failure is recorded and the time step repeated with a smaller size
    if (errorCode != 0 && this->recordFailure_)
        this->recordFailure("AssembleFE_SCI_SMC: AceGen element " + std::to_string(this->getGlobalElementID()) + " failed to compute its state, error code " + std::to_string(errorCode));

    double* residuumRint = elem.getResiduumVectorRint();
    for (int i = 0; i < 30; i++)
        this->residuumRint_[i] = residuumRint[i];

    double* residuumRdyn = elem.getResiduumVectorRdyn();
    for (int i = 0; i < 30; i++)
        this->residuumRdyn_[i] = residuumRdyn[i];

    double* residuumRc = elem.getResiduumVectorRc();
    for (int i = 0; i < 10; i++)
        this->residuumRc_[i] = residuumRc[i];

    if (computeTangent) {
        double** stiffnessMatrixKuu = elem.getStiffnessMatrixKuu();
        for (int i = 0; i < 30; i++)
            for (int j = 0; j < 30; j++)
                this->stiffnessMatrixKuu_[i][j] = stiffnessMatrixKuu[i][j];

        double** stiffnessMatrixKuc = elem.getStiffnessMatrixKuc();
        for (int i = 0; i < 30; i++)
            for (int j = 0; j < 10; j++)
                this->stiffnessMatrixKuc_[i][j] = stiffnessMatrixKuc[i][j];

        double** stiffnessMatrixKcu = elem.getStiffnessMatrixKcu();
        for (int i = 0; i < 10; i++)
            for (int j = 0; j < 30; j++)
                this->stiffnessMatrixKcu_[i][j] = stiffnessMatrixKcu[i][j];

        double** massMatrixMc = elem.getMassMatrixMc();
        for (int i = 0; i < 10; i++)
            for (int j = 0; j < 10; j++)
                this->massMatrixMc_[i][j] = massMatrixMc[i][j];

        double** stiffnessMatrixKcc = elem.getStiffnessMatrixKcc();
        for (int i = 0; i < 10; i++)
            for (int j = 0; j < 10; j++)
                this->stiffnessMatrixKcc_[i][j] = stiffnessMatrixKcc[i][j];

        double* historyUpdated = elem.getHistoryUpdated();
        for (int i = 0; i < this->historyLength_; i++)
            this->historyUpdated_[i] = historyUpdated[i];
    }

#endif
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::postProcessing() {
#ifdef FEDD_HAVE_ACEGENINTERFACE

	double displacements[30];
    for (int i = 0; i < 30; i++)
        displacements[i] = (*this->solution_)[i];
   
    double concentrations[10];
    for (int i = 0; i < 10; i++) {
        concentrations[i] = (*this->solution_)[i + 30];
        solutionC_n1_[i] = (*this->solution_)[i + 30];  // in each newtonstep solution for n+1 is updated.
	}	
	
    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   	
    double rates[10];
    for (int i = 0; i < 10; i++)
        rates[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / this->getTimeIncrement();

    double deltaT = this->getTimeIncrement();

    double time = this->getTimeStep();

    std::vector<double> domainDataModified(this->domainDataLength_);

    std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

    if (time < this->activeAcceleratedEndTime_) {
        for (int idx : acceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
        }
        for (int idx : deceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
        }
	}

    AceGenElement_Type& elem = aceGenElement(&displacements[0], &concentrations[0], &accelerations[0], &rates[0], domainDataModified.data(), deltaT, time);
    
    double** postProcessingResults = elem.postProcess(&displacements[0], &concentrations[0], this->history_.data(), &rates[0], &accelerations[0]);

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < this->postDataLength_; j++)
            (*this->postProcessingData_)[i][j] = postProcessingResults[i][j];
#endif
	}

template <class SC, class LO, class GO, class NO>
int AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::findPosition(const std::string& subString, const std::vector<std::string>& stringArray) {
    for (int i = 0; i < stringArray.size(); i++)
        if (stringArray[i] == subString)
            return i;
	return -1;
	}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::initializeGrowth() {
#ifdef FEDD_HAVE_ACEGENINTERFACE

	double displacements[30];
    for (int i = 0; i < 30; i++)
        displacements[i] = (*this->solution_)[i];

    double concentrations[10];
    for (int i = 0; i < 10; i++) {
        concentrations[i] = (*this->solution_)[i + 30];
        solutionC_n1_[i] = (*this->solution_)[i + 30];  // in each newtonstep solution for n+1 is updated.
	}

    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   
    double rates[10];
    for (int i = 0; i < 10; i++)
        rates[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / this->getTimeIncrement();
   
    double deltaT = this->getTimeIncrement();

    double time = this->getTimeStep();

    std::vector<double> domainDataModified(this->domainDataLength_);

    std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

    if (time < this->activeAcceleratedEndTime_) {
        for (int idx : acceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
	}
        for (int idx : deceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
        }
	}

    AceGenElement_Type& elem = aceGenElement(&displacements[0], &concentrations[0], &accelerations[0], &rates[0], domainDataModified.data(), deltaT, time);
			
    std::vector<double> historyNew = elem.initializeGrowthOrientationVectors();

    for (int i = 0; i < this->historyLength_; i++) {
        this->history_[i] = historyNew[i];
        this->historyUpdated_[i] = historyNew[i];
    }
#endif
    growthInitialized_ = true;
		}
		
template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::initializeActiveResponse() {
    double deltaT = this->getTimeIncrement();
    double time = this->getTimeStep();
#ifdef FEDD_HAVE_ACEGENINTERFACE
    std::vector<double> domainDataModified(this->domainDataLength_);

    std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

    if (time < this->activeAcceleratedEndTime_) {
        for (int idx : acceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
        }
        for (int idx : deceleratedParamIndices_) {
            domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
        }
	}

    AceGenElement_Type& elem = aceGenElement(this->displacements_.data(), this->concentrations_.data(), this->accelerations_.data(), this->rates_.data(), domainDataModified.data(), deltaT, time);
			
    std::vector<double> stretches = elem.getGaussPointStretches();

    int historyPerGP = (int)this->historyLength_ / this->numberOfIntegrationPoints_;
    for (int i = 0; i < this->numberOfIntegrationPoints_; i++) {
        this->history_[i * historyPerGP + 10] = stretches[i * 2];
        this->history_[i * historyPerGP + 11] = stretches[i * 2 + 1];
        this->historyUpdated_[i * historyPerGP + 10] = stretches[i * 2];
        this->historyUpdated_[i * historyPerGP + 11] = stretches[i * 2 + 1];
		}
#endif
    activeInitialized_ = true;
	}
			
#ifdef FEDD_HAVE_ACEGENINTERFACE
template <class SC, class LO, class GO, class NO>
typename AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::AceGenElement_Type&
AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::aceGenElement(double* displacements, double* concentrations, double* accelerations, double* rates, double* domainData, double deltaT, double time) {
#ifdef ACEGENINTERFACE_REUSABLE_ELEMENTS
    // Interface2 resets the element's results in compute(), so one element can be
    // reused with the data of each call.
    if (!aceGenElement_.is_null()) {
        aceGenElement_->setPositions(this->positions_.data());
        aceGenElement_->setDisplacements(displacements);
        aceGenElement_->setConcentrations(concentrations);
        aceGenElement_->setAccelerations(accelerations);
        aceGenElement_->setRates(rates);
        aceGenElement_->setDomainData(domainData);
        aceGenElement_->setHistoryVector(this->history_.data());
        aceGenElement_->setSubIterationTolerance(this->subiterationTolerance_);
        aceGenElement_->setTimeIncrement(deltaT);
        aceGenElement_->setTime(time);
        aceGenElement_->setIntegrationCode(this->iCode_);
        aceGenElement_->setElementID(this->getGlobalElementID());
        return *aceGenElement_;
		}
#endif
    // Older Interface2 versions add the results of compute() to those of the previous
    // call: a new element for every call.
    aceGenElement_.reset(new AceGenElement_Type(this->positions_.data(), displacements, concentrations, accelerations, rates, domainData, this->history_.data(), this->subiterationTolerance_, deltaT, time, this->iCode_, this->getGlobalElementID()));
    return *aceGenElement_;
	}
#endif

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::updateDomainData(const std::string& dataName, double dataValue) {
    int position = findPosition(dataName, this->domainDataNames_);
	TEUCHOS_TEST_FOR_EXCEPTION(position == -1, std::logic_error, " Parameter " << dataName << " not found in domain data names!!");
    this->domainData_[position] = dataValue;
}
	
template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::saveState() {
    AssembleFE<SC, LO, GO, NO>::saveState();
    savedDisplacements_ = this->displacements_;
    savedConcentrations_ = this->concentrations_;
    savedRates_ = this->rates_;
    savedAccelerations_ = this->accelerations_;
    savedSolutionC_n_ = this->solutionC_n_;
    savedSolutionC_n1_ = this->solutionC_n1_;
    savedDomainData_ = this->domainData_;
    savedActiveBool_ = this->activeBool_;
    savedGrowthBool_ = this->growthBool_;
    savedReorientationBool_ = this->reorientationBool_;
    savedActiveInitialized_ = this->activeInitialized_;
    savedGrowthInitialized_ = this->growthInitialized_;
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::restoreState() {
    AssembleFE<SC, LO, GO, NO>::restoreState();
    this->displacements_ = savedDisplacements_;
    this->concentrations_ = savedConcentrations_;
    this->rates_ = savedRates_;
    this->accelerations_ = savedAccelerations_;
    this->solutionC_n_ = savedSolutionC_n_;
    this->solutionC_n1_ = savedSolutionC_n1_;
    this->domainData_ = savedDomainData_;
    this->activeBool_ = savedActiveBool_;
    this->growthBool_ = savedGrowthBool_;
    this->reorientationBool_ = savedReorientationBool_;
    this->activeInitialized_ = savedActiveInitialized_;
    this->growthInitialized_ = savedGrowthInitialized_;
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::synchronizeTime(Teuchos::RCP<TimeSteppingTools> timeSteppingTool) {
    this->timeStep_ = timeSteppingTool->t_;
    this->timeIncrement_ = timeSteppingTool->get_dt();
}

} // namespace FEDD
#endif // AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
