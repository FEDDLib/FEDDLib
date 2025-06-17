#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A, S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"

#include "feddlib/problems/specific/NavierStokesAssFE.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"

/*!
 main of Generalized Newtonian Power Law unit test

 @brief Generalized Newtonian Power Law main

 This generalized newtonian power law unit test compares the current solution of the generalized newtonian power law problem in
 2D to already stored solutions. The stored solutions make sense for the following
 configurations:
 - 2D on 4 procs with P2/ P1 elements (P2 Velocity, P1 Pressure)
 */

void zeroBC(double *x, double *res, double t, const double *parameters) { res[0] = 0.; }
void oneFunc(double *x, double *res, double *parameters) { res[0] = 1.; }

/*For a 2D-Poiseulle-Flow of a Power-Law fluid with eta=K*gamma^(n-1) an analytical solution for the velocity field
  can be derived depending on K, n, H and the pressure gradient dp/dx
  So a simple test case for the generalized-Newtonian fluid solver is a 2D-Poiseuille Flow, prescribing the analytical velocity profile u(y) and
  checking if correct pressure gradient is recovered
*/
void inflowPowerLaw2D(double *x, double *res, double t, const double *parameters)
{

    double K = parameters[0];  // Parameter in Power-Law Model
    double n = parameters[1];  // For n=1.0 we have parabolic inflow profile (Newtonian case)
    double H = parameters[2];  // Height of Channel
    double dp = parameters[3]; // dp/dx constant pressure gradient along channel

    // This corresponds to the analytical solution of a Poiseuille like Plug-flow of a Power-Law fluid
    res[0] = (n / (n + 1.0)) * std::pow(dp / (K), 1.0 / n) * (std::pow(H / (2.0), (n + 1.0) / n) - std::pow(std::abs((H / 2.0) - x[1]), (n + 1.0) / n));
    res[1] = 0.;

    return;
}

void zeroDirichlet2D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;

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
    typedef Teuchos::RCP<Domain<SC, LO, GO, NO>> DomainPtr_Type;

    typedef Matrix<SC, LO, GO, NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    std::string xmlPrecFile;
    std::string xmlProblemFile;

    std::string xmlSolverFile = "parametersSolver.xml";
    xmlPrecFile = "parametersPrec_GNF_2D.xml";
    xmlProblemFile = "parametersProblem_GNF_2D.xml";
    
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");
    myCLP.setOption("problemfile", &xmlProblemFile, ".xml file with Inputparameters.");


    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    {

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        int dim = parameterListProblem->sublist("Parameter").get("Dimension", 2);
        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity", "P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure", "P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type", "unstructured");
        std::string meshName = parameterListProblem->sublist("Mesh Partitioner").get("Mesh 1 Name", "Meshes/rectangle_200.mesh");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter", " ");
        std::string linearization = parameterListProblem->sublist("General").get("Linearization", "Newton");
        std::string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method", "Monolithic");
        std::string FEType = "P2P1";


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        int size = comm->getSize();
        // Inside we construct the mesh so define connectivities etc.
        {
            DomainPtr_Type domainPressure;
            DomainPtr_Type domainVelocity;

            domainPressure.reset(new Domain<SC, LO, GO, NO>(comm, dim));
            domainVelocity.reset(new Domain<SC, LO, GO, NO>(comm, dim));

            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainPressure;

            ParameterListPtr_Type pListPartitioner = sublist(parameterListProblem, "Mesh Partitioner");
            MeshPartitioner<SC, LO, GO, NO> partitionerP1(domainP1Array, pListPartitioner, "P1", dim);

            partitionerP1.readAndPartition();

            if (discVelocity == "P2")
                domainVelocity->buildP2ofP1Domain(domainPressure); // jumps
            else
                domainVelocity = domainPressure;

            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Material").get("PowerLawParameter K", 0.035));
            parameter_vec.push_back(parameterListProblem->sublist("Material").get("PowerLaw index n", 0.7));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Height Inflow", 0.1));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Constant Pressure Gradient", 10.0));

            // Boundary conditions

            Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());

            bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);                // wall
            bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // wall
            bcFactory->addBC(inflowPowerLaw2D, 4, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // Inlet

            NavierStokesAssFE<SC, LO, GO, NO> navierStokesAssFE(domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll);

            {
                MAIN_TIMER_START(NavierStokesAssFE, " AssFE:   Assemble System and solve");
                navierStokesAssFE.addBoundaries(bcFactory);
                navierStokesAssFE.initializeProblem();
                navierStokesAssFE.assemble();

                navierStokesAssFE.setBoundariesRHS();
                NonLinearSolver<SC, LO, GO, NO> nlSolverAssFE(linearization);
                nlSolverAssFE.solve(navierStokesAssFE); 

                MAIN_TIMER_STOP(NavierStokesAssFE);
                comm->barrier();
            }

    
        //Exporting the solution to a .h5 file 
        /*HDF5Export<SC, LO, GO, NO> exporter(navierStokesAssFE.getSolution()->getBlock(0)->getMap(),
                                             "ReferenceSolutions/solution_GNF_velocity_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores"); //  Map and file name
        exporter.writeVariablesHDF5("solution",
                                     navierStokesAssFE.getSolution()->getBlock(0)); // VariableName and Variable

        HDF5Export<SC, LO, GO, NO> exporterP(navierStokesAssFE.getSolution()->getBlock(1)->getMap(),
                                             "ReferenceSolutions/solution_GNF_pressure_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores"); //  Map and file name
         exporterP.writeVariablesHDF5("solution",
                                 navierStokesAssFE.getSolution()->getBlock(1));
         */   

         // We exclude any other tests, than the one prescribed
        TEUCHOS_TEST_FOR_EXCEPTION(!(size == 4), std::logic_error, "The 2D reference solution were generated using 4 processors.");

        HDF5Import<SC, LO, GO, NO> importer(navierStokesAssFE.getSolution()->getBlock(0)->getMap(),
                                            "ReferenceSolutions/solution_GNF_velocity_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores");
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported = importer.readVariablesHDF5("solution");
        HDF5Import<SC, LO, GO, NO> importer_pressure(navierStokesAssFE.getSolution()->getBlock(1)->getMap(),
                                            "ReferenceSolutions/solution_GNF_pressure_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores");
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported_pressure = importer_pressure.readVariablesHDF5("solution");

        // We compare the imported solution to the current one
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solution_GNF_velocity = navierStokesAssFE.getSolution()->getBlock(0);
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solution_GNF_pressure = navierStokesAssFE.getSolution()->getBlock(1);

        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues_velocity = Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(navierStokesAssFE.getSolution()->getBlock(0)->getMap()));
        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues_pressure= Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(navierStokesAssFE.getSolution()->getBlock(1)->getMap()));

        // Calculating the error per node
        // this = alpha*A + beta*B + gamma*this
        errorValues_velocity->update(1., solution_GNF_velocity, -1., solutionImported, 0.);
        // Computing norm
        Teuchos::Array<SC> norm(1);
        errorValues_velocity->normInf(norm);
        double normError_velocity = norm[0];


        errorValues_pressure->update(1., solution_GNF_pressure, -1., solutionImported_pressure, 0.);
        Teuchos::Array<SC> norm_pressure(1);
        errorValues_pressure->normInf(norm_pressure);
        double normError_pressure = norm_pressure[0];



        // Output of error
        if (comm->getRank() == 0) {
            std::cout << " --------------------------------------------------" << std::endl;
            std::cout << "  Error Report " << std::endl;
            std::cout << "   || velocity_current - velocity_stored||_inf = " << normError_velocity << std::endl;
            std::cout << "   || pressure_current - pressure_stored||_inf = " << normError_pressure << std::endl;
            std::cout << " --------------------------------------------------" << std::endl;
        }

        // Throwing exception, if error is too great.
        TEUCHOS_TEST_FOR_EXCEPTION(normError_velocity > 1.e-9 || normError_pressure > 1.e-9 , std::logic_error,
                                    "Difference between current solution and "
                                  "stored solution greater than 1e-9.");

        int boolExportSolution = 0;
        if (boolExportSolution) {
            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());

            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolution = navierStokesAssFE.getSolution()->getBlock(0);

            exPara->setup("solutionGNF", domainVelocity->getMesh(), "P2");

            exPara->addVariable(exportSolution, "u_current", "Vector", dim, domainVelocity->getMapUnique());
            exPara->addVariable(solutionImported, "u_imported", "Vector", dim, domainVelocity->getMapUnique());

            exPara->save(0.0);
        }

        
        

        }
    }
    return (EXIT_SUCCESS);
}
