#ifndef NONLINEARPROBLEM_DEF_hpp
#define NONLINEARPROBLEM_DEF_hpp
#include "NonLinearProblem_decl.hpp"

/*!
 Definition of NonLinearProblem

 @brief  NonLinearProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD
{

    template <class SC, class LO, class GO, class NO>
    NonLinearProblem<SC, LO, GO, NO>::NonLinearProblem(CommConstPtr_Type comm) : Problem<SC, LO, GO, NO>(comm),
                                                                                 previousSolution_(),
                                                                                 residualVec_(),
                                                                                 nonLinearTolerance_(1.e-6),
                                                                                 coeff_(0)
    {
    }

    template <class SC, class LO, class GO, class NO>
    NonLinearProblem<SC, LO, GO, NO>::NonLinearProblem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm) : Problem<SC, LO, GO, NO>(parameterList, comm),
                                                                                                                       previousSolution_(),
                                                                                                                       residualVec_(),
                                                                                                                       nonLinearTolerance_(1.e-6),
                                                                                                                       coeff_(0)
    {
    }
    template <class SC, class LO, class GO, class NO>
    NonLinearProblem<SC, LO, GO, NO>::~NonLinearProblem()
    {
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initializeProblem(int nmbVectors)
    {
        this->system_.reset(new BlockMatrix_Type(1));

        this->initializeVectors(nmbVectors);

        this->initializeVectorsNonLinear(nmbVectors);

        // Init ThyraVectorSpcaes for NOX.
        this->initVectorSpaces();
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initializeVectorsNonLinear(int nmbVectors)
    {

        UN size = this->domainPtr_vec_.size();
        this->previousSolution_.reset(new BlockMultiVector_Type(size));
        this->residualVec_.reset(new BlockMultiVector_Type(size));

        for (UN i = 0; i < size; i++)
        {
            if (this->dofsPerNode_vec_[i] > 1)
            {
                MapConstPtr_Type map = this->domainPtr_vec_[i]->getMapVecFieldUnique();
                MultiVectorPtr_Type prevSolutionPart = Teuchos::rcp(new MultiVector_Type(map));
                this->previousSolution_->addBlock(prevSolutionPart, i);
                MultiVectorPtr_Type residualPart = Teuchos::rcp(new MultiVector_Type(map));
                this->residualVec_->addBlock(residualPart, i);
            }
            else
            {
                MapConstPtr_Type map = this->domainPtr_vec_[i]->getMapUnique();
                MultiVectorPtr_Type prevSolutionPart = Teuchos::rcp(new MultiVector_Type(map));
                this->previousSolution_->addBlock(prevSolutionPart, i);
                MultiVectorPtr_Type residualPart = Teuchos::rcp(new MultiVector_Type(map));
                this->residualVec_->addBlock(residualPart, i);
            }
        }

        this->residualVec_->putScalar(0.);
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::infoNonlinProblem()
    {

        bool verbose(this->comm_->getRank() == 0);
        if (verbose)
        {
            std::cout << "\t ### ### ###" << std::endl;
            std::cout << "\t ### Nonlinear Problem Information ###" << std::endl;
            std::cout << "\t ### Linearization: " << this->parameterList_->sublist("General").get("Linearization", "default") << std::endl;
            std::cout << "\t ### Relative tol: " << this->parameterList_->sublist("Parameter").get("relNonLinTol", 1.e-6) << "\t absolute tol: " << this->parameterList_->sublist("Parameter").get("absNonLinTol", 1.e-4) << "(not used for Newton or Fixed-Point)" << std::endl;
        }
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::calculateNonLinResidualVec(SmallMatrix<double> &coeff, std::string type, double time)
    {

        coeff_ = coeff;
        this->calculateNonLinResidualVec(type, time);
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::reAssembleAndFill(BlockMatrixPtr_Type bMat, std::string type)
    {

        this->assemble(type);
        TEUCHOS_TEST_FOR_EXCEPTION(bMat->size() != this->system_->size(), std::logic_error, "Sizes of BlockMatrices are differen. reAssembleAndFill(...)");

        for (int i = 0; i < this->system_->size(); i++)
        {
            for (int j = 0; j < this->system_->size(); j++)
            {
                if (this->system_->blockExists(i, j))
                {
                    //                MatrixPtr_Type mat = Teuchos::rcp(new Matrix_Type ( this->system_->getBlock( i, j ) ) );
                    bMat->addBlock(this->system_->getBlock(i, j), i, j);
                }
            }
        }
    }

    template <class SC, class LO, class GO, class NO>
    double NonLinearProblem<SC, LO, GO, NO>::calculateResidualNorm() const
    {

        Teuchos::Array<SC> residual(1);
        residualVec_->norm2(residual());
        TEUCHOS_TEST_FOR_EXCEPTION(residual.size() != 1, std::logic_error, "We need to change the code for numVectors>1.");
        return residual[0];
    }

    template <class SC, class LO, class GO, class NO>
    int NonLinearProblem<SC, LO, GO, NO>::solveUpdate()
    {

        // solution COPY!
        *previousSolution_ = *this->solution_;
       
        int its = this->solve(residualVec_);

        return its;
    }

    template <class SC, class LO, class GO, class NO>
    int NonLinearProblem<SC, LO, GO, NO>::solveAndUpdate(const std::string &criterion, double &criterionValue)
    {
        //    BlockMatrixPtr_Type system
        int its = solveUpdate();

        if (criterion == "Update")
        {
            Teuchos::Array<SC> updateNorm(1);
            this->solution_->norm2(updateNorm());
            criterionValue = updateNorm[0];
        }
        this->solution_->update(1., *previousSolution_, 1.);

        return its;
    }

    template <class SC, class LO, class GO, class NO>
    typename NonLinearProblem<SC, LO, GO, NO>::BlockMultiVectorPtr_Type NonLinearProblem<SC, LO, GO, NO>::getResidualVector() const
    {

        return residualVec_;
    }

    template <class SC, class LO, class GO, class NO>
    Thyra::ModelEvaluatorBase::InArgs<SC>
    NonLinearProblem<SC, LO, GO, NO>::getNominalValues() const
    {
        return nominalValues_;
    }

    template <class SC, class LO, class GO, class NO>
    Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> NonLinearProblem<SC, LO, GO, NO>::get_x_space() const
    {
        return xSpace_;
    }

    template <class SC, class LO, class GO, class NO>
    Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> NonLinearProblem<SC, LO, GO, NO>::get_f_space() const
    {
        return fSpace_;
    }

    template <class SC, class LO, class GO, class NO>
    Thyra::ModelEvaluatorBase::InArgs<SC>
    NonLinearProblem<SC, LO, GO, NO>::createInArgs() const
    {
        return prototypeInArgs_;
    }

    // Private functions overridden from ModelEvaulatorDefaultBase

    template <class SC, class LO, class GO, class NO>
    Thyra::ModelEvaluatorBase::OutArgs<SC>
    NonLinearProblem<SC, LO, GO, NO>::createOutArgsImpl() const
    {
        return prototypeOutArgs_;
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initNOXParameters()
    {

        using Teuchos::RCP;
        using Teuchos::rcp;
        using ::Thyra::VectorBase;
        typedef ::Thyra::ModelEvaluatorBase MEB;

        MEB::InArgsSetup<SC> inArgs;
        inArgs.setModelEvalDescription(this->description());
        inArgs.setSupports(MEB::IN_ARG_x);
        this->prototypeInArgs_ = inArgs;

        MEB::OutArgsSetup<SC> outArgs;
        outArgs.setModelEvalDescription(this->description());
        outArgs.setSupports(MEB::OUT_ARG_f);
        outArgs.setSupports(MEB::OUT_ARG_W_op);
        outArgs.setSupports(MEB::OUT_ARG_W_prec);
        this->prototypeOutArgs_ = outArgs;

        this->nominalValues_ = inArgs;
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initVectorSpaces()
    {
        this->initNOXParameters();

        std::string type = this->parameterList_->sublist("General").get("Preconditioner Method", "Monolithic");
        if (!type.compare("Monolithic"))
            initVectorSpacesMonolithic();
        else if (!type.compare("Teko") || type == "FaCSI" || type == "FaCSI-Teko" || type == "Diagonal")
            initVectorSpacesBlock();
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unkown preconditioner/solver type.");
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initVectorSpacesMonolithic()
    {

        BlockMapPtr_Type map = Teuchos::rcp_const_cast<BlockMap_Type>(this->solution_->getMap());

       // Teuchos::RCP<const XTpetra_Type> xTpetraMap = Teuchos::rcp_dynamic_cast<const XTpetra_Type>(map->getMergedMap()->getXpetraMap()->getMap());

        TpetraMapConstPtr_Type tpetraMap = map->getMergedMap()->getTpetraMap();

        this->xSpace_ = Thyra::createVectorSpace<SC, LO, GO, NO>(tpetraMap);
        this->fSpace_ = Thyra::createVectorSpace<SC, LO, GO, NO>(tpetraMap);

        typedef Teuchos::ScalarTraits<SC> ST;
        x0_ = ::Thyra::createMember(this->xSpace_);
        V_S(x0_.ptr(), ST::zero());

        this->nominalValues_.set_x(x0_);
    }

    template <class SC, class LO, class GO, class NO>
    void NonLinearProblem<SC, LO, GO, NO>::initVectorSpacesBlock()
    {

        BlockMapPtr_Type map = Teuchos::rcp_const_cast<BlockMap_Type>(this->solution_->getMap());
    
        TpetraMap_Type tpetra_map;

        Teuchos::Array<ThyraVecSpaceConstPtr_Type> vecSpaceArray(map->size());
        for (int i = 0; i < map->size(); i++)
        {
            //Teuchos::RCP<const XTpetra_Type> xTpetraMap =
            //    Teuchos::rcp_dynamic_cast<const XTpetra_Type>(map->getBlock(i)->getXpetraMap()->getMap());
            TpetraMapConstPtr_Type tpetraMap =map->getBlock(i)->getTpetraMap();
            ThyraVecSpaceConstPtr_Type vecSpace = Thyra::createVectorSpace<SC, LO, GO, NO>(tpetraMap);
            vecSpaceArray[i] = vecSpace;
        }
        this->xSpace_ = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<SC>(vecSpaceArray()));
        this->fSpace_ = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<SC>(vecSpaceArray()));
    
        typedef Teuchos::ScalarTraits<SC> ST;
        x0_ = ::Thyra::createMember(this->xSpace_);
        V_S(x0_.ptr(), ST::zero());

        this->nominalValues_.set_x(x0_);
    }

}

#endif
