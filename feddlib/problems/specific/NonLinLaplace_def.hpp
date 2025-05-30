#ifndef NonLinLaplace_def_hpp
#define NonLinLaplace_def_hpp
#include "NonLinLaplace_decl.hpp"
/*!
 Definition of NonLinLaplace

 @brief NonLinLaplace
 @author Kyrill Ho
 @version 1.0
 @copyright KH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
NonLinLaplace<SC, LO, GO, NO>::NonLinLaplace(
    const DomainConstPtr_Type &domain, std::string FEType,
    ParameterListPtr_Type parameterList)
    : NonLinearProblem<SC, LO, GO, NO>(parameterList, domain->getComm()),
      u_rep_() {
    this->nonLinearTolerance_ =
        this->parameterList_->sublist("Parameter").get("relNonLinTol", 1.0e-6);
    this->initNOXParameters();
    this->addVariable(domain, FEType, "u", 1);
    this->dim_ = this->getDomain(0)->getDimension();
    u_rep_ = Teuchos::rcp(
        new MultiVector_Type(this->getDomain(0)->getMapRepeated()));
}

template <class SC, class LO, class GO, class NO>
NonLinLaplace<SC, LO, GO, NO>::~NonLinLaplace() {}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::info() {
    this->infoProblem();
    this->infoNonlinProblem();
}

/*
 * General assemble function that is called whenever assembly of any kind is
 * required type specifies whether the residual or the Jacobian is needed
 */
template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::assemble(std::string type) const {
    if (type == "") {
        if (this->verbose_) {
            std::cout << "-- Assembly nonlinear laplace ... " << std::flush;
        }
        this->initAssemble();
    } else {
        this->reAssemble(type);
    }
    if (this->verbose_) {
        std::cout << "done -- " << std::endl;
    }
}

/*
 * Initial assembly of the system matrix (Jacobian)
 */
template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::initAssemble() const {

    if (this->verbose_) {
        std::cout << "-- Initial assembly " << std::flush;
    }

    if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(1));

    if (this->residualVec_.is_null())
        this->residualVec_.reset(new BlockMultiVector_Type(1));

    this->feFactory_->assemblyNonlinearLaplace(
        this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
        this->system_, this->residualVec_, this->parameterList_, "Jacobian");

    // Initialise solution to 1 everywhere
    this->solution_->putScalar(1.);

    if (this->verbose_) {
        std::cout << "done -- " << std::endl;
    }
}

/*
 * Reassembly of the residual or Jacobian
 * Required by the nonlinear solver
 */
template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::reAssemble(std::string type) const {

    if (this->verbose_)
        std::cout << "-- Reassembly nonlinear laplace"
                  << " (" << type << ") ... " << std::flush;
    // Update the locally stored solution to the problem
    MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
    this->u_rep_->importFromVector(u, true);

    if (type == "Rhs") {

        this->feFactory_->assemblyNonlinearLaplace(
            this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
            this->system_, this->residualVec_, this->parameterList_, "Rhs");

    } else if (type == "Newton") {

        this->feFactory_->assemblyNonlinearLaplace(
            this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
            this->system_, this->residualVec_, this->parameterList_,
            "Jacobian");
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::reAssembleExtrapolation(
    BlockMultiVectorPtrArray_Type previousSolutions) {

    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Only Newton/NOX implemented for nonlinear Laplace!");
}


/*
 * Updates the residual vector by reassembling and applying boundary conditions
 */
template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::calculateNonLinResidualVec(
    std::string type, double time) const {

    this->reAssemble("Rhs");
    // rhs_ contains the dirichlet boundary conditions
    // Apparently so does bcFactory_. Not sure why both do.
    if (!type.compare("standard")) {
        // this = 1*this - 1*rhs
        this->residualVec_->update(-1., *this->rhs_, 1.);
        this->bcFactory_->setVectorMinusBC(this->residualVec_, this->solution_,
                                           time);
    } else if (!type.compare("reverse")) {
        // this = -1*this + 1*rhs

        this->residualVec_->update(1., *this->rhs_, -1.);

        // Sets the residualVec_ at the boundary = boundary condition - solution
        // Necessary since reAssemble("Rhs") only assembles the residual on
        // internal nodes
        this->bcFactory_->setBCMinusVector(this->residualVec_, this->solution_,
                                           time);
    } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                                   "Unknown type for residual computation.");
    }
}
} // namespace FEDD
#endif
