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

        this->newtonStep_=0;
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
        else if (!type.compare("Teko") || type == "FaCSI" || type == "FaCSI-Teko" || type == "Diagonal" || type == "Triangular" || type == "PCD"|| type == "LSC")
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


    template<class SC,class LO,class GO,class NO>
    void NonLinearProblem<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                                ) const
    {
        std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        if ( !type.compare("Monolithic"))
            evalModelImplMonolithic( inArgs, outArgs );
        else if ( !type.compare("Teko")|| !type.compare("Diagonal") || type == "Triangular"|| !type.compare("PCD") || !type.compare("LSC")){
    #ifdef FEDD_HAVE_TEKO
            evalModelImplBlock( inArgs, outArgs );
    #else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
    #endif
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
    }

    /*!
        \brief Monolithic Approach for Nonlinear Solver NOX. Input. Includes calculation of the residual vector and update (reAssembly) of non constant matrices with new solution.
            ResidualVec and SystemMatrix of this class are then converted into the corresponding Thyra/Tpetra objects for Solver.
    */
    template<class SC,class LO,class GO,class NO>
    void NonLinearProblem<SC,LO,GO,NO>::evalModelImplMonolithic(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                            const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
    {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_dynamic_cast;
        using Teuchos::rcp_const_cast;
        using Teuchos::ArrayView;
        using Teuchos::Array;
        RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");

        RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
        RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

        RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);

        this->solution_->fromThyraMultiVector(vecThyraNonConst);

        const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
        const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
        const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();

        typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
        typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
        typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
        typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;
    
        const bool fill_f = nonnull(f_out);
        const bool fill_W = nonnull(W_out);
        const bool fill_W_prec = nonnull(W_prec_out);


        if ( fill_f || fill_W || fill_W_prec ) {

            // ****************
            // Get the underlying xpetra objects
            // ****************
            if (fill_f) {

                this->calculateNonLinResidualVec("standard"); // Calculating residual Vector

                // Changing the residualVector into a ThyraMultivector

                Teuchos::RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
                f_out->assign(*f_thyra);
            }

            TpetraMatrixPtr_Type W;
            if (fill_W) {
                this->assemble("Newton"); // ReAssembling matrices with updated u  in this class

                this->setBoundariesSystem(); // setting boundaries to the system

                Teuchos::RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
                Teuchos::RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
                
                TpetraMatrixConstPtr_Type W_systemTpetra = this->getSystem()->getMergedMatrix()->getTpetraMatrix();           
                TpetraMatrixPtr_Type W_systemTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_systemTpetra);
                
                //Tpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
                //Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                
                Teuchos::RCP<TpetraMatrix_Type> tpetraMatTpetra = W_systemTpetraNonConst; //xTpetraMat.getTpetra_CrsMatrixNonConst();
                
                W_tpetraMat->resumeFill();

                for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
                    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                    tpetraMatTpetra->getLocalRowView( i, indices, values);
                    W_tpetraMat->replaceLocalValues( i, indices, values);
                }
                W_tpetraMat->fillComplete();

            }

            if (fill_W_prec) {
                
                if (precInitOnly_){
                    int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);
                    if(this->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                    {
                        this->setupPreconditioner( "Monolithic" );
                    }
                    else{
                        if (this->verbose_)
                            std::cout << " NonLinearProblem<SC,LO,GO,NO>::evalModelImplMonolithic:: Skipping preconditioner reconstruction" << std::endl;
                    }
                }
                else
                    precInitOnly_ = true;

                // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
                // If this is the case, we need to pre apply the coarse level to the residual(f_out).

                std::string levelCombination = this->parameterList_->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive");
                if (!levelCombination.compare("Multiplicative")) {
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX.");
                }

                this->newtonStep_ ++; 

            }
        }
    }
    /*!
        \brief Block Approach for Nonlinear Solver NOX. Input. Includes calculation of the residual vector and update (reAssembly) of non constant matrices with new solution.
            ResidualVec and SystemMatrix of this class are then converted into the corresponding Thyra/Tpetra objects for Solver.
    */
    #ifdef FEDD_HAVE_TEKO
    template<class SC,class LO,class GO,class NO>
    void NonLinearProblem<SC,LO,GO,NO>::evalModelImplBlock(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                    const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
    {
        using Teuchos::RCP;
        using Teuchos::rcp;
        using Teuchos::rcp_dynamic_cast;
        using Teuchos::rcp_const_cast;
        using Teuchos::ArrayView;
        using Teuchos::Array;

        RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");

        RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
        RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

        RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);

        RCP< Thyra::ProductVectorBase< SC > > vecThyraBlock = rcp_dynamic_cast<Thyra::ProductVectorBase< SC > > (vecThyraNonConst);

        this->solution_->getBlockNonConst(0)->fromThyraMultiVector( vecThyraBlock->getNonconstVectorBlock(0) );
        this->solution_->getBlockNonConst(1)->fromThyraMultiVector( vecThyraBlock->getNonconstVectorBlock(1) );

        const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
        const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
        const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();

        typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
        typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
        typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
        typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;

        const bool fill_f = nonnull(f_out);
        const bool fill_W = nonnull(W_out);
        const bool fill_W_prec = nonnull(W_prec_out);

        if ( fill_f || fill_W || fill_W_prec ) {

            // ****************
            // Get the underlying xpetra objects
            // ****************
            if (fill_f) {

                this->calculateNonLinResidualVec("standard");

                Teko::MultiVector f0;
                Teko::MultiVector f1;
                f0 = this->getResidualVector()->getBlockNonConst(0)->getThyraMultiVector();
                f1 = this->getResidualVector()->getBlockNonConst(1)->getThyraMultiVector();

                std::vector<Teko::MultiVector> f_vec; f_vec.push_back(f0); f_vec.push_back(f1);

                Teko::MultiVector f = Teko::buildBlockedMultiVector(f_vec);

                f_out->assign(*f);
            }

            TpetraMatrixPtr_Type W;
            if (fill_W) {

                typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;

                this->assemble("Newton");

                this->setBoundariesSystem();

                RCP<ThyraBlockOp_Type> W_blocks = rcp_dynamic_cast<ThyraBlockOp_Type>(W_out);
                RCP<const ThyraOp_Type> W_block00 = W_blocks->getBlock(0,0);
                RCP<ThyraOp_Type> W_block00NonConst = rcp_const_cast<ThyraOp_Type>( W_block00 );
                RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator( W_block00NonConst );

                RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

                TpetraMatrixConstPtr_Type W_matrixTpetra = this->getSystem()->getBlock(0,0)->getTpetraMatrix();
                TpetraMatrixPtr_Type W_matrixTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_matrixTpetra);
                RCP<TpetraMatrix_Type> tpetraMatTpetra = W_matrixTpetraNonConst;

                W_tpetraMat->resumeFill();

                for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
                    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                    tpetraMatTpetra->getLocalRowView( i, indices, values);
                    W_tpetraMat->replaceLocalValues( i, indices, values);
                }
                W_tpetraMat->fillComplete();

            }

            if (fill_W_prec) {
                std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
                if (precInitOnly_){
                    int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);
                    if(this->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                    {
                        this->setupPreconditioner( type );
                    }
                    else{
                        if (this->verbose_)
                            std::cout << " \n NonLinearProblem<SC,LO,GO,NO>::evalModelImplBlock:: Skipping preconditioner reconstruction \n " << std::endl;
                    }
                }
                else
                    precInitOnly_ = true;
                // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
                // If this is the case, we need to pre apply the coarse level to the residual(f_out).

                ParameterListPtr_Type tmpSubList = sublist( sublist( sublist( sublist( this->parameterList_, "Teko Parameters" ) , "Preconditioner Types" ) , "Teko" ) , "Inverse Factory Library" );

                std::string levelCombination1 = tmpSubList->sublist( "FROSch-Velocity" ).get("Level Combination","Additive");
                std::string levelCombination2 = tmpSubList->sublist( "FROSch-Pressure" ).get("Level Combination","Additive");

                if ( !levelCombination1.compare("Multiplicative") || !levelCombination2.compare("Multiplicative") ) {

                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX.");
                    ParameterListPtr_Type solverPList = this->getLinearSolverBuilder()->getNonconstParameterList();
                }
                this->newtonStep_ ++; 

            }
        }
    }
    #endif

        
    template<class SC,class LO,class GO,class NO>
    Teuchos::RCP<Thyra::LinearOpBase<SC> > NonLinearProblem<SC,LO,GO,NO>::create_W_op() const
    {
        std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        if ( !type.compare("Monolithic"))
            return create_W_op_Monolithic( );
        else if ( !type.compare("Teko") || !type.compare("Diagonal") || !type.compare("Triangular") || !type.compare("PCD") || !type.compare("LSC") ){
    #ifdef FEDD_HAVE_TEKO
            return create_W_op_Block( );
    #else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
    #endif
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
    }

    template<class SC,class LO,class GO,class NO>
    Teuchos::RCP<Thyra::LinearOpBase<SC> > NonLinearProblem<SC,LO,GO,NO>::create_W_op_Monolithic() const
    {
        Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
        Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
        return W_op;
    }

    #ifdef FEDD_HAVE_TEKO
    template<class SC,class LO,class GO,class NO>
    Teuchos::RCP<Thyra::LinearOpBase<SC> > NonLinearProblem<SC,LO,GO,NO>::create_W_op_Block() const
    {

        Teko::LinearOp thyraF = this->system_->getBlock(0,0)->getThyraLinOp();
        Teko::LinearOp thyraBT = this->system_->getBlock(0,1)->getThyraLinOp();
        Teko::LinearOp thyraB = this->system_->getBlock(1,0)->getThyraLinOp();

        if (!this->system_->blockExists(1,1)){
            MatrixPtr_Type dummy = Teuchos::rcp( new Matrix_Type( this->system_->getBlock(1,0)->getMap(), 1 ) );
            dummy->fillComplete();
            this->system_->addBlock( dummy, 1, 1 );
        }

        Teko::LinearOp thyraC = this->system_->getBlock(1,1)->getThyraLinOp();

        Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = Thyra::block2x2(thyraF,thyraBT,thyraB,thyraC);
        Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
        return W_op;
    }
    #endif

    template<class SC,class LO,class GO,class NO>
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > NonLinearProblem<SC,LO,GO,NO>::create_W_prec() const
    {
        this->initializeSolverBuilder();

        std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        this->setBoundariesSystem();

        if (!type.compare("Teko") || !type.compare("Diagonal") || !type.compare("Triangular") || !type.compare("PCD") || !type.compare("LSC")) { //
            this->setupPreconditioner( type );
            precInitOnly_ = false;
        }
        else{
            this->setupPreconditioner( type ); // initializePreconditioner( type );
            precInitOnly_ = false;
        }
        

        Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
        Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst = Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);

        return thyraPrecNonConst;

    }

}

#endif
