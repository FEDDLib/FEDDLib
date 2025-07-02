#ifndef PROBLEM_DEF_hpp
#define PROBLEM_DEF_hpp
#include "Problem_decl.hpp"

/*!
 Definition of Problem
 @brief  Problem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using Teuchos::outArg;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;

namespace FEDD
{
    template <class SC, class LO, class GO, class NO>
    Problem<SC, LO, GO, NO>::Problem(CommConstPtr_Type comm) : dim_(-1),
                                                               comm_(comm),
                                                               system_(),
                                                               rhs_(),
                                                               solution_(),
                                                               preconditioner_(),
                                                               linearSolverBuilder_(),
                                                               verbose_(comm->getRank() == 0),
                                                               parameterList_(),
                                                               domainPtr_vec_(),
                                                               domain_FEType_vec_(),
                                                               variableName_vec_(),
                                                               bcFactory_(),
                                                               feFactory_(),
                                                               dofsPerNode_vec_(),
                                                               sourceTerm_(),
                                                               rhsFuncVec_(),
                                                               parasSourceFunc_(0)
#ifdef FEDD_TIMER
                                                               ,
                                                               solveProblemTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Solve")), bcMatrixTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Set Boundaries Matrix")), bcRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Set Boundaries RHS"))
#endif
    {
        linearSolverBuilder_.reset(new Stratimikos::DefaultLinearSolverBuilder());
        preconditioner_ = Teuchos::rcp(new Preconditioner_Type(this));
        feFactory_.reset(new FEFac_Type());
    }

    template <class SC, class LO, class GO, class NO>
    Problem<SC, LO, GO, NO>::Problem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm) : dim_(-1),
                                                                                                     comm_(comm),
                                                                                                     system_(),
                                                                                                     rhs_(),
                                                                                                     solution_(),
                                                                                                     preconditioner_(),
                                                                                                     linearSolverBuilder_(),
                                                                                                     verbose_(comm->getRank() == 0),
                                                                                                     parameterList_(parameterList),
                                                                                                     domainPtr_vec_(),
                                                                                                     domain_FEType_vec_(),
                                                                                                     variableName_vec_(),
                                                                                                     bcFactory_(),
                                                                                                     feFactory_(),
                                                                                                     dofsPerNode_vec_(),
                                                                                                     sourceTerm_(),
                                                                                                     rhsFuncVec_(),
                                                                                                     parasSourceFunc_(0)
#ifdef FEDD_TIMER
                                                                                                     ,
                                                                                                     solveProblemTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Solve")), bcMatrixTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Set Boundaries Matrix")), bcRHSTimer_(Teuchos::TimeMonitor::getNewCounter("FEDD - Problem - Set Boundaries RHS"))
#endif
    {
        linearSolverBuilder_.reset(new Stratimikos::DefaultLinearSolverBuilder());
        preconditioner_ = Teuchos::rcp(new Preconditioner_Type(this));
        feFactory_.reset(new FEFac_Type());
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::infoProblem()
    {
        bool verbose(comm_->getRank() == 0);
        if (verbose)
        {
            std::cout << "\t ### Problem Information ###" << std::endl;
            std::cout << "\t ### Dimension: " << dim_ << std::endl;
            std::cout << "\t ### Linear problem: "
                      << "??" << std::endl;
            std::cout << "\t ### Number of blocks/equations/variables: " << domainPtr_vec_.size() << std::endl;
            for (int i = 0; i < domainPtr_vec_.size(); i++)
            {
                std::cout << "\t \t # Block " << i + 1 << "\t name: " << variableName_vec_.at(i) << "\t d.o.f.s: " << dofsPerNode_vec_.at(i) << "\t FE type: " << domain_FEType_vec_.at(i) << std::endl;
            }

            // ch 15.04.19: Hier ggf. unterscheiden zwischen Monolithic und Teko bzw. anderen Block-Precs.
            ParameterListPtr_Type pListThyraPrec = sublist(parameterList_, "ThyraPreconditioner");
            std::cout << "\t ### ### ###" << std::endl;
            std::cout << "\t ### Preconditioner Information ###" << std::endl;
            std::cout << "\t ### Type: " << parameterList_->sublist("General").get("Preconditioner Method", "Monolithic") << std::endl;
            std::cout << "\t ### Prec.: " << pListThyraPrec->get("Preconditioner Type", "FROSch") << std::endl;

            if (!pListThyraPrec->get("Preconditioner Type", "FROSch").compare("FROSch") && parameterList_->sublist("General").get("Preconditioner Method", "Monolithic") == "Monolithic")
            {
                std::cout << "\t ### Variant: " << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("FROSch Preconditioner Type", "TwoLevelBlockPreconditioner") << std::endl;
                std::cout << "\t ### Two Level: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("TwoLevel", false)
                          << "\t Overlap: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("Overlap", 0)
                          << "\t Level Combination: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination", "Additive") << std::endl;

                std::cout << "\t OverlappingOperator Type: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("OverlappingOperator Type", "AlgebraicOverlappingOperator") << std::endl;

                std::cout << "\t CoarseOperator Type: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type", "GDSWCoarseOperator") << std::endl;

                for (int i = 0; i < this->parameterList_->get("Number of blocks", 2); i++)
                {
                    std::cout << "\t \t # Block " << i + 1 << "\t d.o.f.s: "
                              << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("DofsPerNode" + std::to_string(i + 1), 1)
                              << "\t d.o.f. ordering: "
                              << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("DofOrdering" + std::to_string(i + 1), "NodeWise") << std::endl;
                }
            }
            else
            {
                std::cout << "\t ### Full preconditioner information only available for Monolithic preconditioner type ###" << std::endl;
            }
        }
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::initializeProblem(int nmbVectors)
    {

        this->system_.reset(new BlockMatrix_Type(1));
        this->initializeVectors(nmbVectors);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::addRhsFunction(RhsFunc_Type func)
    {
        this->rhsFuncVec_.push_back(func);
    }
    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::addRhsFunction(RhsFunc_Type func, int i)
    {
        if (this->rhsFuncVec_.size() <= i + 1)
            this->rhsFuncVec_[i] = func;
        else if (this->rhsFuncVec_.size() == i)
            this->rhsFuncVec_.push_back(func);
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Insertion Index smaller then size of vector");
    }

    /*template<class SC,class LO,class GO,class NO>
    void Problem<SC,LO,GO,NO>::addRhsFunctionAndFlag(RhsFunc_Type func, int i, int flag){

        TEUCHOS_TEST_FOR_EXCEPTION(!parameterList_->sublist("Parameter").get("Source Type","volume").compare("volume") , std::runtime_error, "Source Type Error: RHS Flags only make sence with surface source type.");

        if(this->rhsFuncVec_.size() <= i+1){
            this->rhsFuncVec_[i] = func;
            this->rhsFuncFlagVec_[i] = flag;
        }
        else if(this->rhsFuncVec_.size() == i){
            this->rhsFuncVec_.push_back(func);
            this->rhsFuncFlagVec_.push_back(flag);
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Insertion Index smaller then size of vector");


    }*/

    template <class SC, class LO, class GO, class NO>
    RhsFunc_Type &Problem<SC, LO, GO, NO>::getRhsFunction(int i)
    {
        return rhsFuncVec_[i];
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::addVariable(const DomainConstPtr_Type &domain, std::string FEType, std::string name, int dofsPerNode)
    {

        domainPtr_vec_.push_back(domain);
        domain_FEType_vec_.push_back(FEType);
        variableName_vec_.push_back(name);
        feFactory_->addFE(domain);
        dofsPerNode_vec_.push_back(dofsPerNode);
    }

    // template<class SC,class LO,class GO,class NO>
    // void Problem<SC,LO,GO,NO>::reAssemble(){
    //     if (verbose_) {
    //         cout << "Nothing to reassemble for linear problem." << endl;
    //     }
    // }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::assembleSourceTerm(double time) const
    {

        TEUCHOS_TEST_FOR_EXCEPTION(sourceTerm_.is_null(), std::runtime_error, "Initialize source term before you assemble it - sourceTerm pointer is null");
        this->sourceTerm_->putScalar(0.);
        std::string sourceType = parameterList_->sublist("Parameter").get("Source Type", "volume");

        if (sourceType == "volume")
            assembleVolumeTerm(time);
        else if (sourceType == "surface")
            assembleSurfaceTerm(time);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::assembleVolumeTerm(double time) const
    {
        for (UN i = 0; i < sourceTerm_->size(); i++)
        {
            if (this->rhsFuncVec_.size() > i)
            {
                if (!this->rhsFuncVec_[i].empty())
                {
                    MultiVectorPtr_Type FERhs;
                    // funcParameter[0] is always the time
                    vec_dbl_Type funcParameter(1, 0.);
                    funcParameter[0] = time;

                    // how can we use different parameters for different blocks here?
                    for (int j = 0; j < parasSourceFunc_.size(); j++)
                        funcParameter.push_back(parasSourceFunc_[j]);

                    std::string type;
                    if (this->getDofsPerNode(i) > 1)
                    {
                        FERhs = Teuchos::rcp(new MultiVector_Type(this->domainPtr_vec_.at(i)->getMapVecFieldRepeated()));
                        type = "Vector";
                    }
                    else
                    {
                        FERhs = Teuchos::rcp(new MultiVector_Type(this->domainPtr_vec_.at(i)->getMapRepeated()));
                        type = "Scalar";
                    }

                    this->feFactory_->assemblyRHS(this->dim_,
                                                  this->domain_FEType_vec_.at(i),
                                                  FERhs,
                                                  type,
                                                  this->rhsFuncVec_[i],
                                                  funcParameter);

                    this->sourceTerm_->getBlockNonConst(i)->exportFromVector(FERhs, false, "Add");
                }
            }
        }
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::assembleSurfaceTerm(double time) const
    {
        for (UN i = 0; i < sourceTerm_->size(); i++)
        {
            if (this->rhsFuncVec_.size() > i)
            {

                if (!this->rhsFuncVec_[i].empty())
                {
                    MultiVectorPtr_Type FERhs;
                    // funcParameter[0] is always the time
                    vec_dbl_Type funcParameter(1, 0.);
                    funcParameter[0] = time;

                    // how can we use different parameters for different blocks here?
                    for (int j = 0; j < parasSourceFunc_.size(); j++)
                        funcParameter.push_back(parasSourceFunc_[j]);

                    // we add an additional parameter to place the surface flag of the element there during the assembly and shift the degree of the function to the last place now
                    funcParameter.push_back(funcParameter[funcParameter.size() - 1]);
                    std::string type;
                    if (this->getDofsPerNode(i) > 1)
                    {
                        FERhs = Teuchos::rcp(new MultiVector_Type(this->domainPtr_vec_.at(i)->getMapVecFieldRepeated()));
                        type = "Vector";
                    }
                    else
                    {
                        FERhs = Teuchos::rcp(new MultiVector_Type(this->domainPtr_vec_.at(i)->getMapRepeated()));
                        type = "Scalar";
                    }

                    this->feFactory_->assemblySurfaceIntegral(this->getDomain(i)->getDimension(),
                                                              this->getDomain(i)->getFEType(),
                                                              FERhs,
                                                              "Vector",
                                                              this->rhsFuncVec_[i],
                                                              funcParameter);

                    this->sourceTerm_->getBlockNonConst(i)->exportFromVector(FERhs, false, "Add");
                }
            }
        }
        //    this->sourceTerm_->scale(-1.); // this scaling is needed for TPM problem
    }

    template <class SC, class LO, class GO, class NO>
    int Problem<SC, LO, GO, NO>::solve(BlockMultiVectorPtr_Type rhs)
    {
        int its;
        if (verbose_)
            std::cout << "-- Solve System ..." << std::endl;
        {

#ifdef FEDD_TIMER
            TimeMonitor_Type solveTM(*solveProblemTimer_);
#endif
            LinearSolver<SC, LO, GO, NO> linSolver;
            std::string type = parameterList_->sublist("General").get("Preconditioner Method", "Monolithic");
            its = linSolver.solve(this, rhs, type); // if rhs is null. Then the rhs_ of this is used in the linear solver
        }

        if (verbose_)
            std::cout << " done -- " << std::endl;

        return its;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setupPreconditioner(std::string type) const
    {

        preconditioner_->buildPreconditioner(type);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::initializePreconditioner(std::string type) const
    {

        preconditioner_->initializePreconditioner(type);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::addBoundaries(const BCConstPtr_Type &bcFactory)
    {

        bcFactory_ = bcFactory;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setBoundaries(double time) const
    {
#ifdef FEDD_TIMER
        TimeMonitor_Type bcMatTM(*bcMatrixTimer_);
        TimeMonitor_Type bcRHSTM(*bcRHSTimer_);
#endif
        bcFactory_->set(system_, rhs_, time);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setBoundariesRHS(double time) const
    {
#ifdef FEDD_TIMER
        TimeMonitor_Type bcRHSTM(*bcRHSTimer_);
#endif
        bcFactory_->setRHS(rhs_, time);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setAllDirichletZero(BlockMultiVectorPtr_Type u) const
    {
#ifdef FEDD_TIMER
        TimeMonitor_Type bcRHSTM(*bcRHSTimer_);
#endif
        bcFactory_->setAllDirichletZero(u);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setBoundariesSystem() const
    {
#ifdef FEDD_TIMER
        TimeMonitor_Type bcMatTM(*bcMatrixTimer_);
#endif
        bcFactory_->setSystem(system_);
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::initializeVectors(int nmbVectors)
    {

        UN size = domainPtr_vec_.size();
        solution_.reset(new BlockMultiVector_Type(size));
        rhs_.reset(new BlockMultiVector_Type(size));
        sourceTerm_.reset(new BlockMultiVector_Type(size));
        rhsFuncVec_.resize(size);

        for (UN i = 0; i < size; i++)
        {
            if (dofsPerNode_vec_[i] > 1)
            {
                MapConstPtr_Type map = domainPtr_vec_[i]->getMapVecFieldUnique();
                MultiVectorPtr_Type solutionPart = Teuchos::rcp(new MultiVector_Type(map));
                solution_->addBlock(solutionPart, i);
                MultiVectorPtr_Type rhsPart = Teuchos::rcp(new MultiVector_Type(map));
                rhs_->addBlock(rhsPart, i);
                MultiVectorPtr_Type sourceTermPart = Teuchos::rcp(new MultiVector_Type(map));
                sourceTerm_->addBlock(sourceTermPart, i);
            }
            else
            {
                MapConstPtr_Type map = domainPtr_vec_[i]->getMapUnique();
                MultiVectorPtr_Type solutionPart = Teuchos::rcp(new MultiVector_Type(map));
                solution_->addBlock(solutionPart, i);
                MultiVectorPtr_Type rhsPart = Teuchos::rcp(new MultiVector_Type(map));
                rhs_->addBlock(rhsPart, i);
                MultiVectorPtr_Type sourceTermPart = Teuchos::rcp(new MultiVector_Type(map));
                sourceTerm_->addBlock(sourceTermPart, i);
            }
        }
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BlockMultiVectorPtr_Type Problem<SC, LO, GO, NO>::getRhs()
    {

        return rhs_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BlockMultiVectorPtr_Type Problem<SC, LO, GO, NO>::getRhs() const
    {

        return rhs_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BlockMultiVectorPtr_Type Problem<SC, LO, GO, NO>::getSolution()
    {

        return solution_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BlockMatrixPtr_Type Problem<SC, LO, GO, NO>::getSystem() const
    {

        return system_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::PreconditionerPtr_Type Problem<SC, LO, GO, NO>::getPreconditioner()
    {
        return preconditioner_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::PreconditionerConstPtr_Type Problem<SC, LO, GO, NO>::getPreconditionerConst() const
    {
        return preconditioner_;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::setPreconditionerThyraFromLinOp(ThyraLinOpPtr_Type precLinOp)
    {
        preconditioner_->setPreconditionerThyraFromLinOp(precLinOp);
    }

    template <class SC, class LO, class GO, class NO>
    bool Problem<SC, LO, GO, NO>::getVerbose() const
    {

        return verbose_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::FEFacConstPtr_Type Problem<SC, LO, GO, NO>::getFEFactory()
    {

        return feFactory_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BCConstPtr_Type Problem<SC, LO, GO, NO>::getBCFactory()
    {

        return bcFactory_;
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::DomainConstPtr_Type Problem<SC, LO, GO, NO>::getDomain(int i) const
    {

        return domainPtr_vec_.at(i);
    }

    template <class SC, class LO, class GO, class NO>
    std::string Problem<SC, LO, GO, NO>::getFEType(int i) const
    {

        return domain_FEType_vec_.at(i);
    }

    template <class SC, class LO, class GO, class NO>
    std::string Problem<SC, LO, GO, NO>::getVariableName(int i) const
    {

        return variableName_vec_.at(i);
    }

    template <class SC, class LO, class GO, class NO>
    int Problem<SC, LO, GO, NO>::getDofsPerNode(int i) const
    {

        return dofsPerNode_vec_.at(i);
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::ParameterListPtr_Type Problem<SC, LO, GO, NO>::getParameterList() const
    {

        return parameterList_;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::addToRhs(BlockMultiVectorPtr_Type x) const
    {

        rhs_->update(1., *x, 1.);
    }

    template <class SC, class LO, class GO, class NO>
    typename Problem<SC, LO, GO, NO>::BlockMultiVectorPtr_Type Problem<SC, LO, GO, NO>::getSourceTerm()
    {
        TEUCHOS_TEST_FOR_EXCEPTION(sourceTerm_.is_null(), std::runtime_error, "Problem has no source term to return - source term pointer is null");
        return sourceTerm_;
    }

    template <class SC, class LO, class GO, class NO>
    bool Problem<SC, LO, GO, NO>::hasSourceTerm() const
    {
        for (int i = 0; i < rhsFuncVec_.size(); i++)
        {
            if (!rhsFuncVec_[i].empty())
                return true;
        }
        return false;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::initSolutionWithVector(MultiVector_Type &mv)
    {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "initSolutionWithVector not implemented. DO we need this?");

        //    *solution_ = mv;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::initializeSolverBuilder() const
    {

        ParameterListPtr_Type pListThyraPrec = sublist(parameterList_, "ThyraPreconditioner");

        linearSolverBuilder_->setParameterList(pListThyraPrec);
    }

    // Functions that return the H1 and L2 Norm of a given Vector. The Norms being defined as:
    // || mv ||_L2 = mv^T * M * mv
    // || mv ||_H1 = mv^T * K * mv
    // with M beeing the mass matrix and K beeing the stiffness matrix

    template <class SC, class LO, class GO, class NO>
    double Problem<SC, LO, GO, NO>::calculateL2Norm(MultiVectorConstPtr_Type mv, int domainInd)
    {

        MatrixPtr_Type M;
        M = Teuchos::rcp(new Matrix_Type(this->domainPtr_vec_.at(domainInd)->getMapUnique(), this->getDomain(domainInd)->getApproxEntriesPerRow()));
        this->feFactory_->assemblyMass(this->dim_, this->domain_FEType_vec_.at(domainInd), "Scalar", M);

        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> mvOutput = Teuchos::rcp(new MultiVector_Type(this->domainPtr_vec_.at(domainInd)->getMapUnique()));

        M->apply(*mv, *mvOutput);

        Teuchos::ArrayRCP<SC> vector = mv->getDataNonConst(0);
        Teuchos::ArrayRCP<SC> outputVector = mvOutput->getDataNonConst(0);

        double result = 0;
        for (int i = 0; i < vector.size(); i++)
        {
            result += vector[i] * outputVector[i];
        }
        reduceAll<int, double>(*comm_, REDUCE_SUM, result, outArg(result));

        return result;
    }
    template <class SC, class LO, class GO, class NO>
    double Problem<SC, LO, GO, NO>::calculateH1Norm(MultiVectorConstPtr_Type mv, int blockId1, int blockId2, int domainInd)
    {

        MatrixPtr_Type K = this->getSystem()->getBlock(blockId1, blockId2);

        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> mvOutput = Teuchos::rcp(new MultiVector_Type(K->getMap()));

        K->apply(*mv, *mvOutput); // this represents mvOutput = K * mv ;

        Teuchos::ArrayRCP<SC> vector = mv->getDataNonConst(0);
        Teuchos::ArrayRCP<SC> outputVector = mvOutput->getDataNonConst(0);

        double result = 0;
        for (int i = 0; i < vector.size(); i++)
        {
            result += vector[i] * outputVector[i];
        }
        reduceAll<int, double>(*comm_, REDUCE_SUM, result, outArg(result));

        return result;
    }

    template <class SC, class LO, class GO, class NO>
    void Problem<SC, LO, GO, NO>::infoParameter()
    {
        bool verbose(comm_->getRank() == 0);
        if (verbose)
        {
            std::cout << "#####################################" << std::endl;
            std::cout << " ### Problem Information ###" << std::endl;
            std::cout << " ### Dimension: " << dim_ << std::endl;
            std::cout << " ### Number of blocks/equations/variables: " << domainPtr_vec_.size() << std::endl;
            for (int i = 0; i < domainPtr_vec_.size(); i++)
            {
                std::cout << "  ## Block " << i + 1 << " name: " << variableName_vec_.at(i) << " d.o.f.s: " << dofsPerNode_vec_.at(i) << " FE type: " << domain_FEType_vec_.at(i) << std::endl;
            }
            std::cout << " ####################################" << std::endl;
            ParameterListPtr_Type parameterlist = sublist(parameterList_, "Parameter");
            std::cout << " ### Parameter Information ###" << std::endl;
            std::cout << " ### Mesh Type: " << parameterlist->get("Mesh Type", "???") << std::endl;
            if(parameterlist->get("Mesh Type", "???") == "structured" || parameterlist->get("Mesh Type", "???") == "structured_bfs" )
                std::cout << " ### H/h= " << parameterlist->get("H/h", 0) << std::endl;
            else
                std::cout << " ### Mesh File Name 1: " << parameterList_->sublist("Mesh Partitioner").get("Mesh 1 Name", "???") << std::endl;

            std::cout << " ### Viscosity: " << parameterlist->get("Viscosity", 0.) 
                      << " ### Density: " << parameterlist->get("Density", 0.) << std::endl;

            if(abs(parameterlist->get("MaxVelocity", 0.)) > 0 )
                std::cout << " ### Maximum Velocity: " << parameterlist->get("MaxVelocity", 0.) << std::endl;   
            else if(abs(parameterlist->get("Max Velocity", 0.)) > 0 )
                std::cout << " ### or Maximum Velocity: " << parameterlist->get("Max Velocity", 0.) << std::endl;   

            std::cout << " ### Rel. Tol.: " << parameterlist->get("relNonLinTol", 0.) 
                << " ### Abs. Tol.: " << parameterlist->get("absNonLinTol", 0.) << std::endl;    

            std::cout << " ####################################" << std::endl;

            
            // ch 15.04.19: Hier ggf. unterscheiden zwischen Monolithic und Teko bzw. anderen Block-Precs.
            ParameterListPtr_Type pListThyraPrec = sublist(parameterList_, "ThyraPreconditioner");
            std::cout << " ### Preconditioner Information ###" << std::endl;
            std::cout << "  ## Type: " << parameterList_->sublist("General").get("Preconditioner Method", "Monolithic") << std::endl;

            if (!pListThyraPrec->get("Preconditioner Type", "FROSch").compare("FROSch") && parameterList_->sublist("General").get("Preconditioner Method", "Monolithic") == "Monolithic")
            {
                std::cout << "  ## Variant: " << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("FROSch Preconditioner Type", "TwoLevelBlockPreconditioner") << std::endl;
                std::cout << "  ## Overlap: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("Overlap", 0) << std::endl;
                std::cout << "  ## OverlappingOperator Type: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("OverlappingOperator Type", "AlgebraicOverlappingOperator") << std::endl;

                std::cout << "  ## CoarseOperator Type: "
                          << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type", "GDSWCoarseOperator") << std::endl;

                if(pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type", "GDSWCoarseOperator") == "IPOUHarmonicCoarseOperator"){
                    for (int i = 0; i < this->parameterList_->get("Number of blocks", 2); i++)
                    {
                        std::cout << "   # IPOU Block "<< std::to_string(i + 1) <<":  " << pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").sublist("IPOUHarmonicCoarseOperator").sublist("Blocks").sublist(std::to_string(i + 1)).sublist("InterfacePartitionOfUnity").get("Type","NOTFOUND") << std::endl;
                    }
                }
                    
            }
            else if (!parameterList_->sublist("General").get("Preconditioner Method", "Monolithic").compare("Teko")){
                std::cout << "  ### Block Preconditioner Type: " << parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE") << std::endl;
                std::cout << "  ### Velocity Preconditioner:   " << parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("FROSch-Velocity").get("CoarseOperator Type","GDSW#") << std::endl;
                std::cout << "  ### Pressure Preconditioner:   " << parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("FROSch-Pressure").get("CoarseOperator Type","GDSW#") << std::endl;

            }
            else if (!parameterList_->sublist("General").get("Preconditioner Method", "Monolithic").compare("Diagonal") ||
                !parameterList_->sublist("General").get("Preconditioner Method", "Monolithic").compare("Triangular") ||
                !parameterList_->sublist("General").get("Preconditioner Method", "Monolithic").compare("PCD") ||
                !parameterList_->sublist("General").get("Preconditioner Method", "Monolithic").compare("LSC"))
            {
                // ParameterListPtr_Type pListThyraPrecBlock = sublist(parameterList_, "Block Preconditioner");

                std::cout << "  ### Block Preconditioner Type: " << parameterList_->sublist("General").get("Preconditioner Method", "Monolithic") << std::endl;
                std::cout << "   ## Velocity Preconditioner:   " << parameterList_->sublist("Velocity preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type", "Nada") << std::endl;
                std::cout << "    # Variant: " << parameterList_->sublist("Velocity preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("FROSch Preconditioner Type", "TwoLevelBlockPreconditioner") << std::endl;
                std::cout << "    # Overlap: "
                          << parameterList_->sublist("Velocity preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Overlap", 0) << std::endl;
               
                std::cout << "   ## Schur Complement Preconditioner: " << parameterList_->sublist("Schur complement preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type", "Nada") << std::endl;
                std::cout << "    # Variant: " << parameterList_->sublist("Schur complement preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("FROSch Preconditioner Type", "TwoLevelBlockPreconditioner") << std::endl;
                std::cout << "    # Overlap: "
                          << parameterList_->sublist("Schur complement preconditioner").sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Overlap", 0) << std::endl;
            }
            else
            {
                std::cout << " ### Full preconditioner information only available for Monolithic/Teko preconditioner type ###" << std::endl;
            }
            std::cout << " ####################################" << std::endl;

            std::cout << " ### Git commit hash: " << GIT_COMMIT_HASH << std::endl;
            std::cout << " ### Git branch name: " << GIT_BRANCH_NAME << std::endl;

            std::cout << "#####################################" << std::endl;

        }
    }


}
#endif
