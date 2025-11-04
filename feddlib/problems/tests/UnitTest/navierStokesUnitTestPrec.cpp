#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"

#include "feddlib/problems/specific/NavierStokes.hpp"

/*!
 main of Navier Stokes unit test

 @brief Navier Stokes main

 This stokes unit test compares the current solution of the stokes problem in
 2D and 3D to older already stored solutions. The folling test constructs a
 structured mesh with H/h = 4. The stored solutions make sense for the following
 configurations:
 - 2D on 4 procs with P1 or P2 elements
 - 3D on 8 procs with P1 or P2 elements

 */

void zeroBC(double *x, double *res, double t, const double *parameters) { res[0] = 0.; }
void oneFunc(double *x, double *res, double *parameters) { res[0] = 1.; }

void zeroDirichlet(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

// For Lid Driven Cavity Test
void ldcFunc2D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.*parameters[0];
    res[1] = 0.;
    
    return;
}

// For Lid Driven Cavity Test
void ldcFunc3D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.*parameters[0];
    res[1] = 0.;
    res[2] = 0.;

    return;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    int boolExportSolution = 0;
    myCLP.setOption("exportSolution", &boolExportSolution, "Export Solution");

    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    std::string FETypeV = "P2";
    myCLP.setOption("FEType", &FETypeV, "Discretization");

    std::string precMethod = "Monolithic";
    myCLP.setOption("precMethod", &precMethod, "Preconditioning Method");
    
    std::string FEType = "P1";

    std::string xmlSolverFile = "parametersSolver_PrecTest.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);

    
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }
    
    std::string xmlPrecFile; 
    if(precMethod=="Monolithic")
        xmlPrecFile = "parametersPrec_NavierStokes_Mono.xml";
    else if(precMethod=="Diagonal" || precMethod=="Triangular" || precMethod=="PCD" || precMethod=="LSC")
        xmlPrecFile = "parametersPrec_NavierStokes_Block.xml";
    else if(precMethod=="Teko")
        xmlPrecFile = "parametersPrec_NavierStokes_Teko.xml";

    {
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListPrec));
        parameterListAll->setParameters(*parameterListSolver);
        parameterListAll->setParameters(*parameterListPrec);
        // Mesh
        int m = 3;
        std::string meshType = "structured";
        std::string meshDelimiter = " ";

        int n;
        int numProcsCoarseSolve =0;
        int size = comm->getSize();

        int minNumberSubdomains = 1;

        //We exclude any other tests, than the one prescribed
        TEUCHOS_TEST_FOR_EXCEPTION(!(size == 9 && m == 3), std::logic_error, "The 2D test solutions are only sensible for 4 processors.");
       
        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainPressure;
        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainVelocity;

        if (dim == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(std::floor(std::sqrt(size)) != std::sqrt(size), std::logic_error, "Wrong number of processors for structured squared mesh.");
            n = (int)(std::pow(size / minNumberSubdomains, 1 / 2.) + 100 * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(2);
            x[0] = 0.0;
            x[1] = 0.0;
            domainPressure.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
            domainVelocity.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
        } 
        else 
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only 2D Test.");
            
        domainPressure->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        domainVelocity->buildMesh(1, "Square", dim, FETypeV, n, m, numProcsCoarseSolve);

        // ####################
        std::vector<double> parameter_vec(1, 1.0); // Velocity
        parameter_vec.push_back(1.);               // height of inflow region

        // Lid Driven Cavity Test 
        Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());
        bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim); // wall
        bcFactory->addBC(ldcFunc2D, 2, 0, domainVelocity, "Dirichlet", dim,parameter_vec); // lid
        bcFactory->addBC(zeroDirichlet, 3, 1, domainPressure, "Dirichlet", 1); // pressure node


        // Setting parameters that are needed for stokes 
        // !! This avoids that errors occure when default values are changed along the way.
        parameterListAll->sublist("Parameter").set("Density",1.0); // Value 1.0
        parameterListAll->sublist("Parameter").set("Viscosity",1.0e-2); // Value 0.01

        parameterListAll->sublist("Parameter").set("relNonLinTol",1.0e-8); // Value 1.0e-8
        parameterListAll->sublist("General").set("Preconditioner Method",precMethod); // Value 1.0e-8

        vec_dbl_ptr_Type its = Teuchos::rcp(new vec_dbl_Type ( 2, 0. ) ); //0:linear iterations, 1: nonlinear iterations

        NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, FETypeV, domainPressure, FEType, parameterListAll );

        {
            navierStokes.addBoundaries(bcFactory);
            navierStokes.initializeProblem();
            navierStokes.assemble();
            navierStokes.setBoundaries();

            std::string nlSolverType = "Newton";
            NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
            nlSolver.solve( navierStokes,its );
        }

        // std::cout << " (*its)[0] " << (*its)[0] << std::endl;
       
        // HDF5Export<SC, LO, GO, NO> exporterV(navierStokes.getSolution()->getBlock(0)->getMap(),
        //     "ReferenceSolutions/solution_navier_stokes_velocity_" + precMethod +"_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores"); //  Map and file name
        // exporterV.writeVariablesHDF5("velocity",
        //     navierStokes.getSolution()->getBlock(0)); // VariableName and Variable

        // HDF5Export<SC, LO, GO, NO> exporterP(navierStokes.getSolution()->getBlock(1)->getMap(),
        //     "ReferenceSolutions/solution_navier_stokes_pressure_" + precMethod +"_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores"); //  Map and file name
        // exporterP.writeVariablesHDF5("pressure",
        //     navierStokes.getSolution()->getBlock(1)); // VariableName and Variable


        HDF5Import<SC, LO, GO, NO> importerV(navierStokes.getSolution()->getBlock(0)->getMap(),
            "ReferenceSolutions/solution_navier_stokes_velocity_" + precMethod +"_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores"); //  Map and file name
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImportedVeloctiy = importerV.readVariablesHDF5("velocity");

        HDF5Import<SC, LO, GO, NO> importerP(navierStokes.getSolution()->getBlock(1)->getMap(),
            "ReferenceSolutions/solution_navier_stokes_pressure_" + precMethod +"_" + std::to_string(dim) + "d_" + FETypeV + "_" + std::to_string(size) + "cores"); //  Map and file name
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImportedPressure = importerP.readVariablesHDF5("pressure");

        // We compare the imported solution to the current one
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionStokesVelocity = navierStokes.getSolution()->getBlock(0);
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionStokesPressure = navierStokes.getSolution()->getBlock(1);

        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValuesVelocity = Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(navierStokes.getSolution()->getBlock(0)->getMap()));
        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValuesPressure = Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(navierStokes.getSolution()->getBlock(1)->getMap()));

        // this = alpha*A + beta*B + gamma*this
        errorValuesVelocity->update(1., solutionStokesVelocity, -1., solutionImportedVeloctiy, 0.);
        errorValuesPressure->update(1., solutionStokesPressure, -1., solutionImportedPressure, 0.);

        // Computing norm
        Teuchos::Array<SC> normV(1);
        Teuchos::Array<SC> normP(1);

        errorValuesVelocity->normInf(normV);
        errorValuesPressure->normInf(normP);

        double normErrorV = normV[0];
        double normErrorP = normP[0];

        solutionStokesVelocity->normInf(normV);
        solutionStokesPressure->normInf(normP);

        double normVe = normV[0];
        double normPr = normP[0];
        // Output of error
        if (comm->getRank() == 0) {
            std::cout << " --------------------------------------------------" << std::endl;
            std::cout << "  Error Report " << std::endl;
            std::cout << "   || velocity_current - velocity_stored||_inf = " << normErrorV << std::endl;
            std::cout << "   || pressure_current - pressure_stored||_inf = " << normErrorP << std::endl;
            std::cout << "   || velocity_current - velocity_stored||_inf/|| velocity_current||_inf = " << normErrorV/normVe << std::endl;
            std::cout << "   || pressure_current - pressure_stored||_inf/|| pressure_current ||_inf = " << normErrorP/normPr << std::endl;
            std::cout << "   || velocity_current||_inf = " << normVe << std::endl;
            std::cout << "   || pressure_current ||_inf = " << normPr << std::endl;
            std::cout << " --------------------------------------------------" << std::endl;
        }

        // Throwing exception, if error is too great.
        TEUCHOS_TEST_FOR_EXCEPTION(normErrorV/normVe > 1.e-9 || normErrorP/normPr > 1.e-9 , std::logic_error,
                                    "Difference between current solution and "
                                    "stored solution greater than 1e-9.");

        if(precMethod=="Monolithic"){
            TEUCHOS_TEST_FOR_EXCEPTION((*its)[0] - 18.25 > 1.e-9 , std::logic_error,
                                    "Itertion count for Monolithic changed compared to initial test.");
        }                                
        else if(precMethod=="Diagonal"){
            TEUCHOS_TEST_FOR_EXCEPTION((*its)[0] - 82.0 > 1.e-9 , std::logic_error,
                                    "Itertion count for Diagonal changed compared to initial test.");
        }                                    
        else if(precMethod=="PCD"){
            TEUCHOS_TEST_FOR_EXCEPTION((*its)[0] - 73.25 > 1.e-9 , std::logic_error,
                                    "Itertion count for PCD changed compared to initial test.");
        }
        else if(precMethod=="Teko"){
            TEUCHOS_TEST_FOR_EXCEPTION((*its)[0] - 21.25 > 1.e-9 , std::logic_error,
                                    "Itertion count for Teko changed compared to initial test.");
        }

        if (boolExportSolution) {
            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());

            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolution = navierStokes.getSolution()->getBlock(1);

            exPara->setup("solutionStokes", domainPressure->getMesh(), FEType);

            exPara->addVariable(exportSolution, "p_current", "Scalar", 1, domainPressure->getMapUnique());
            exPara->addVariable(solutionImportedPressure, "p_imported", "Scalar", 1, domainPressure->getMapUnique());

            exPara->save(0.0);
        }
    }
    return (EXIT_SUCCESS);
}
