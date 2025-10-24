#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A, S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokesAssFE.hpp"


#include <boost/function.hpp>

/*!
 Main of steady-state Generalized Newtonian fluid flow problem with generalized Newtonian shear stress tensor assumption
 where we use, for example, the Power-law viscosity model eta( gamma_dot) = K* gamma_dot(n-1)

***********************************************************************************************


 Rectangular grid present for Poiseuille-Flow using, e.g., a Power-Law GNF model: Set the analytical velocity profile for a Power-Law fluid and see if correct pressure gradient is recovered


 @brief steady-state generalized-Newtonian Flow main
 @author Natalie Kubicki
 @version 1.0
 @copyright NK
 */

using namespace std;
using namespace Teuchos;
using namespace FEDD;

// These are the boundary conditions already defined
void zeroDirichlet(double *x, double *res, double t, const double *parameters)
{

    res[0] = 0.0;
    return;
}
void zeroDirichlet2D(double *x, double *res, double t, const double *parameters)
{

    res[0] = 0.;
    res[1] = 0.;

    return;
}

void one(double *x, double *res, double t, const double *parameters)
{

    res[0] = 1.;
    res[1] = 1.;

    return;
}

void onex(double *x, double *res, double t, const double *parameters)
{

    res[0] = 1.;
    res[1] = 0.;

    return;
}

void constx3D(double *x, double *res, double t, const double *parameters)
{

    res[0] = 0.1;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void oney(double *x, double *res, double t, const double *parameters)
{

    res[1] = 1.;
    res[0] = 0.;

    return;
}

void two(double *x, double *res, double t, const double *parameters)
{

    res[0] = 2.;
    res[1] = 2.;

    return;
}
void three(double *x, double *res, double t, const double *parameters)
{

    res[0] = 3.;
    res[1] = 3.;
    return;
}
void four(double *x, double *res, double t, const double *parameters)
{
    res[0] = 4.;
    res[1] = 4.;
    return;
}

void zeroDirichlet3D(double *x, double *res, double t, const double *parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

/*For a 2D-Poiseulle-Flow of a Power-Law fluid with eta=K*gamma^(n-1) an analytical solution for the velocity field
  can be derived depending on K, n, H and the pressure gradient dp/dx
  So a simple test case for the generalized-Newtonian fluid solver is a 2D-Poiseuille Flow, prescribing the analytical velocity profile u(y) and
  checking if correct pressure gradient is recovered
*/
void inflowPowerLaw2D(double *x, double *res, double t, const double *parameters)
{

    double K = parameters[0];  //
    double n = parameters[1];  // For n=1.0 we have parabolic inflow profile (Newtonian case)
    double H = parameters[2];  // Height of Channel
    double dp = parameters[3]; // dp/dx constant pressure gradient along channel

    // This corresponds to the analytical solution of a Poiseuille like Plug-flow of a Power-Law fluid
    res[0] = (n / (n + 1.0)) * std::pow(dp / (K), 1.0 / n) * (std::pow(H / (2.0), (n + 1.0) / n) - std::pow(std::abs((H / 2.0) - x[1]), (n + 1.0) / n));
    res[1] = 0.;

    return;
}

void inflowParabolic2D(double *x, double *res, double t, const double *parameters)
{

    double H = parameters[2];
    double v_max = parameters[4];
    res[0] = 4 * v_max * x[1] * (H - x[1]) / (H * H);
    res[1] = 0.;

    return;
}

void inflowParabolic3D(double *x, double *res, double t, const double *parameters)
{

    double H = parameters[2];
    double mu = 0.035;
    double dp = parameters[3];

    res[0] = (1 / (2 * mu)) * (-1.0) * dp * (x[1] * x[1] - x[1] * H);
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void dummyFunc(double *x, double *res, double t, const double *parameters)
{

    return;
}



// Add here predefined switching strategies for FixedPoint-Newtons Method

// Simple switch to the other linearization method after a defined number of nonlinear iterations
bool switchingAfterMaxIter(std::string& linearization, int nlIts, double criterionValue, ParameterListPtr_Type parameterList )
{
    if (nlIts == parameterList->get("MaximumIteration_Switch", 1) )
    {
        if( (parameterList->get("Initial_Linearization", "FixedPoint") == "FixedPoint") && (linearization == "FixedPoint") )   
        {
            linearization = "Newton";
            return true;
        }
        else if(parameterList->get("Initial_Linearization", "FixedPoint")  == "Newton" && (linearization == "Newton") )
        {
            linearization = "FixedPoint";
            return true;
        }
        else
        {
            // Linearization not changed
            return false;
        }
    }
    else // Stay with the current linearization if number of iterations to switch not yet reached
    {
        return false;
    }
}

// Switch to the other linearization method if the nonlinear criterion value is below a certain tolerance
bool switchingAfterNonlinearTolerance(std::string& linearization, int nlIts, double criterionValue, ParameterListPtr_Type parameterList)
{
    // Switch to Newtons Method if relative residual is below a certain value
    if (criterionValue < parameterList->get("NonlinearTolerance_Switch", 0.1))
    {
        if( (parameterList->get("Initial_Linearization", "FixedPoint") == "FixedPoint") && (linearization == "FixedPoint") )   
        {
            linearization = "Newton";
            return true;
        }
        else if(parameterList->get("Initial_Linearization", "FixedPoint")  == "Newton" && (linearization == "Newton") )
        {
            linearization = "FixedPoint";
            return true;
        }
        else
        {
            // Linearization not changed
            return false;
        }
    }
    else // Stay with the current linearization if tolerance not yet reached
    {
        return false;
    }
}


// Switch linearization method after every iteration
bool switchingAlternating(std::string& linearization, int nlIts, double criterionValue, ParameterListPtr_Type parameterList)
{
    if (nlIts > 0) // Switch only after first iteration, because one solution step should be done with the initial linearization strategy
    {
        if (linearization == "FixedPoint")
        {
            linearization = "Newton";
            return true;
        }
        else if (linearization == "Newton")
        {
            linearization = "FixedPoint";
            return true;
        }
        else
        {
            std::cerr << "Error: Current Linearization not FixedPoint or Newton - Staying with current one!" << std::endl;
            return false;
        }
    }
    else
    {
        return false; // Do not change in first iteration, because we want to start one computation with the initial linearization strategy
    }
}



typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

// Now the main starts
int main(int argc, char *argv[])
{

    typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC, LO, GO, NO>> DomainPtr_Type;

    typedef Matrix<SC, LO, GO, NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    // boost function type for switching strategy: 
    typedef boost::function<bool( std::string& , int, double, ParameterListPtr_Type)> SwitchingStrategyFunc;

    
    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    bool verbose(comm->getRank() == 0);

    if (verbose)
    {
        cout << "###############################################################" << endl;
        cout << "##################### Steady Flow of Generalized Newtonian Fluid ####################" << endl;
        cout << "###############################################################" << endl;
    }

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    string xmlProblemFile = "parametersProblem_PicardNewton.xml";
    myCLP.setOption("problemfile", &xmlProblemFile, ".xml file with Inputparameters.");

    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");

    string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile", &xmlTekoPrecFile, ".xml file with Inputparameters.");

    double length = 4.;
    myCLP.setOption("length", &length, "length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        return EXIT_SUCCESS;
    }
    // Einlesen von Parameterwerten
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);

        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListPrecTeko = Teuchos::getParametersFromXmlFile(xmlTekoPrecFile);

        int dim = parameterListProblem->sublist("Parameter").get("Dimension", 3);

        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity", "P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure", "P1");

        string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type", "structured");
        string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name", "circle2D_1800.mesh");
        string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter", " ");
        int m = parameterListProblem->sublist("Parameter").get("H/h", 5);
        string linearization = parameterListProblem->sublist("General").get("Linearization", "FixedPoint");

        string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method", "Monolithic");
        int mixedFPIts = parameterListProblem->sublist("General").get("MixedFPIts", 1);
        int n;

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        else
            parameterListAll->setParameters(*parameterListPrecTeko);
        parameterListAll->setParameters(*parameterListSolver);

        int minNumberSubdomains;
        if (!meshType.compare("structured"))
        {
            minNumberSubdomains = 1;
        }
        else if (!meshType.compare("structured_bfs"))
        {
            minNumberSubdomains = (int)2 * length + 1;
        }

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse", 0);
        int size = comm->getSize() - numProcsCoarseSolve;

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

            //          **********************  BOUNDARY CONDITIONS ***********************************
            std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type", "parabolic");
            // We can here differentiate between different problem cases and boundary conditions

            /***** For boundary condition read parameter values */
            // So parameter[0] is K, [1] is n, [2] is H, [3] is dp/dx
            /*******************************************************************************/
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Material").get("PowerLawParameter K", 0.));
            parameter_vec.push_back(parameterListProblem->sublist("Material").get("PowerLaw index n", 1.));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Height Inflow", 1.));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("Constant Pressure Gradient", 1.));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter").get("MaxVelocity", 1.));


            /*** Switching Strategy Parameters ******************************************************************************************************************************/
            // Define switching strategy function based on user input -  <Parameter name="SwitchingStrategy" type="string" value="AfterMaxIter"/>  Predefined switching strategies are: AfterMaxIter, AfterResidual -->
            std::string switchingStrategyType = parameterListProblem->sublist("General").get("SwitchingStrategy", "AfterMaxIter"); // Default is AfterMaxIter
            SwitchingStrategyFunc switchingStrategy = []( std::string& currentLinearization , int nlIts, double criterionValue, ParameterListPtr_Type parameterList ) {return false;};
            if (switchingStrategyType == "AfterMaxIter")
            {
                if (verbose) std::cout << "Selected Switching Strategy for FixedPointNewton: " << switchingStrategyType << std::endl;
                switchingStrategy = switchingAfterMaxIter;
            }
            else if (switchingStrategyType == "Alternating")
            {
                if (verbose) std::cout << "Selected Switching Strategy for FixedPointNewton: " << switchingStrategyType << std::endl;
                switchingStrategy = switchingAlternating;
            }
            else if (switchingStrategyType == "AfterResidual")
            {
                if (verbose) std::cout << "Selected Switching Strategy for FixedPointNewton: " << switchingStrategyType << std::endl;
                switchingStrategy = switchingAfterNonlinearTolerance;
            }
            else
            {
                if (verbose) std::cerr << "Error: Selected Switching Strategy for FixedPointNewton not recognized - No switching will be performed!" << std::endl;
            }



            /*******************************************************************************/

            Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());

            if (dim == 2)
            {
                // Flags for rectangular channel inside this folder rectangular_200.mesh
                //            1
                //   ------------------
                //   |                |
                // 4 |                | 3
                //   |                |
                //   ------------------
                //            2
                //************* POISEUILLE FLOW  *************
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);                // wall
                bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // wall

                bcFactory->addBC(inflowPowerLaw2D, 4, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // Analytical solution Power-Law 2D Poiseuille flow based on parameters and pressure gradient
            }
            else if (dim == 3)
            {


            }

            //          **********************  CALL SOLVER ***********************************
            NavierStokesAssFE<SC, LO, GO, NO> navierStokesAssFE(domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll);

            {
                MAIN_TIMER_START(NavierStokesAssFE, " AssFE:   Assemble System and solve");
                navierStokesAssFE.addBoundaries(bcFactory);
                navierStokesAssFE.initializeProblem();
                navierStokesAssFE.assemble();

                navierStokesAssFE.setBoundariesRHS();

                std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization", "FixedPoint");

                NonLinearSolver<SC, LO, GO, NO> nlSolverAssFE(nlSolverType);

                nlSolverAssFE.addSwitchingStrategyFunction( switchingStrategy); // Add switching strategy function of NonLinearSolver to User defined one

                nlSolverAssFE.solve(navierStokesAssFE); // jumps into NonLinearSolver_def.hpp

                MAIN_TIMER_STOP(NavierStokesAssFE);
                comm->barrier();
            }

            //          **********************  POST-PROCESSING ***********************************
            Teuchos::TimeMonitor::report(cout, "Main");

            //****************************************************************************************
            //          **********************  POST-PROCESSING - VISCOSITY COMPUTATION FOR GENERALIZEDN-NEWTONIAN FLUID ***********************************
            // We only write out viscosity field if we consider generalized-Newtonian fluid because otherwise it is constant
            if ((parameterListProblem->sublist("Material").get("Newtonian", true) == false) && (parameterListProblem->sublist("Material").get("WriteOutViscosity", false)) == true)
            {

                Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exParaViscsoity(new ExporterParaView<SC, LO, GO, NO>());
                DomainPtr_Type domV = domainVelocity;

                navierStokesAssFE.computeSteadyPostprocessingViscosity_Solution();

                //**************** Write out viscosity ******************
                Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolutionViscosityAssFE = navierStokesAssFE.viscosity_element_;
                exParaViscsoity->setup("viscosity", domV->getMesh(), "P0"); // Viscosity averaged therefore P0 value
                exParaViscsoity->addVariable(exportSolutionViscosityAssFE, "viscosityAssFE", "Scalar", 1, domV->getElementMap());
                exParaViscsoity->save(0.0);
		        exParaViscsoity->closeExporter();
            }

            //****************************************************************************************
            //          **********************  POST-PROCESSING - WRITE OUT VELOCITY AND PRESSURE ***********************************
            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exParaVelocity(new ExporterParaView<SC, LO, GO, NO>());
            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exParaPressure(new ExporterParaView<SC, LO, GO, NO>());

            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolutionVAssFE = navierStokesAssFE.getSolution()->getBlock(0);
            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolutionPAssFE = navierStokesAssFE.getSolution()->getBlock(1);

            DomainPtr_Type dom = domainVelocity;
            exParaVelocity->setup("velocity", dom->getMesh(), dom->getFEType());
            UN dofsPerNode = dim;
            exParaVelocity->addVariable(exportSolutionVAssFE, "uAssFE", "Vector", dofsPerNode, dom->getMapUnique());

            dom = domainPressure;
            exParaPressure->setup("pressure", dom->getMesh(), dom->getFEType());
            exParaPressure->addVariable(exportSolutionPAssFE, "pAssFE", "Scalar", 1, dom->getMapUnique());

            exParaVelocity->save(0.0);
            exParaPressure->save(0.0);

	        exParaVelocity->closeExporter();
	        exParaPressure->closeExporter();


            //**************************** Plot Subdomains without overlap ****************************************
            if (parameterListAll->sublist("General").get("ParaView export subdomains", false))
            {

                if (verbose)
                    std::cout << "\t### Exporting fluid subdomains ###\n";

                typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
                typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
                typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
                typedef BlockMultiVector<SC, LO, GO, NO> BlockMultiVector_Type;
                typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

                {
                    MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type(domainVelocity->getElementMap()));
                    MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                    vecDecomposition->putScalar(comm->getRank() + 1.);

                    Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());

                    exPara->setup("subdomains_fluid", domainVelocity->getMesh(), "P0");

                    exPara->addVariable(vecDecompositionConst, "subdomains", "Scalar", 1, domainVelocity->getElementMap());
                    exPara->save(0.0);
                    exPara->closeExporter();
                }
            }

        }
    }
    Teuchos::TimeMonitor::report(cout);
    }
    
    return EXIT_SUCCESS;
}
