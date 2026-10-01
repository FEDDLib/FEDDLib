#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include "feddlib/core/General/AceGenInterfaceCheck.hpp"
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>
#include "feddlib/core/General/HDF5Import.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

void reactionTerm(double *, double *, double *);
void loadFunction(double *, double *, double *);
void zeroDirichlet3D(double *, double *, double, const double *);
void inflowChem(double *, double *, double, const double *);
using namespace FEDD;


int main(int argc, char *argv[])
{
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
    if (!FEDD::aceGenInterfaceAvailable())
        return EXIT_FAILURE;

    Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    Teuchos::CommandLineProcessor commandLineProcessor;
    std::string underlyingLibrary = "Tpetra";
    std::string simulationParametersXML = "simulationParameters.xml";
    std::string materialParametersXML = "materialParameters_dan.xml";
    std::string solverParametersXML = "solverParameters.xml";
    std::string structurePreconditionerParametersXML = "preconditionerParameters_Structure.xml";
    std::string chemistryPreconditionerParametersXML = "preconditionerParameters_Chemistry.xml";
    commandLineProcessor.setOption("underlyingLibrary", &underlyingLibrary, "Underlying Library");
    commandLineProcessor.setOption("simulationsParameters", &simulationParametersXML, "xml file with simulation parameters");
    commandLineProcessor.setOption("materialParameters", &materialParametersXML, "xml file with material parameters");
    commandLineProcessor.setOption("solverParameters", &solverParametersXML, "xml file with solver parameters");
    commandLineProcessor.setOption("preconditionerParametersStructure", &structurePreconditionerParametersXML, "xml file with structure preconditoner parameters");
    commandLineProcessor.setOption("preconditionerParametersChemistry", &chemistryPreconditionerParametersXML, "xml file with chemistry preconditioner parameters");

    commandLineProcessor.recogniseAllOptions(true);
    commandLineProcessor.throwExceptions(true);

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = commandLineProcessor.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    Teuchos::RCP<Teuchos::StackedTimer> stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Structure-chemical interaction", true));
    bool verbose(comm->getRank() == 0);
    bool testPassed = true;

    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
    {
        Teuchos::RCP<Teuchos::ParameterList> simulationParameters = Teuchos::getParametersFromXmlFile(simulationParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> materialParameters = Teuchos::getParametersFromXmlFile(materialParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> solverParameters = Teuchos::getParametersFromXmlFile(solverParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> structurePreconditionerParameters = Teuchos::getParametersFromXmlFile(structurePreconditionerParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> chemistryPreconditionerParamerters = Teuchos::getParametersFromXmlFile(chemistryPreconditionerParametersXML); // Parameters for solving Chemistry explicit

        int dimension = simulationParameters->sublist("Simulation Parameters").get("Dimension", 3);
        std::string discretizationType = simulationParameters->sublist("Simulation Parameters").get("Discretization", "P2");
        // string preconditionerType = simulationParameters->sublist("Simulation Parameters").get("Preconditioner Type");

        Teuchos::RCP<Teuchos::ParameterList> allParameters = Teuchos::rcp(new Teuchos::ParameterList(*simulationParameters));

        allParameters->sublist("Parameter").set("Chemistry Explicit", false );  // We set chemistry explicit to false again here, since this main only considers chem explicit 

        Teuchos::RCP<Teuchos::ParameterList> preconditionerParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));

        allParameters->setParameters(*materialParameters); // Adding Material Parameters
        allParameters->setParameters(*structurePreconditionerParameters); // Adding Preconditioning Parameters
        allParameters->setParameters(*solverParameters); // Adding Solver Parameters
        //allParameters->setParameters(*chemistryPreconditionerParamerters); // Adding Chemistry Preconditioner Parameters

        Teuchos::RCP<Teuchos::ParameterList> allDiffusionParameters = Teuchos::rcp(new Teuchos::ParameterList(*chemistryPreconditionerParamerters));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Chem"));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Simulation Parameters"));
        allDiffusionParameters->setParameters(*solverParameters);
        allDiffusionParameters->setParameters(*chemistryPreconditionerParamerters);

        Teuchos::RCP<Teuchos::ParameterList> allStructureParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));
        Teuchos::sublist(allStructureParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Solid"));
        allStructureParameters->setParameters(*materialParameters); // Adding Material Parameters
        allStructureParameters->setParameters(*solverParameters); // Adding Material Parameters
      
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Diffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Structure;
        // Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP2Diffusion;
        // Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP2Structure;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainDiffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainStructure;

        domainP1Diffusion.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainP1Structure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainDiffusion.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));
        domainStructure.reset(new FEDD::Domain<SC, LO, GO, NO>(comm, dimension));

        std::vector<Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>>> domainP1Array(1);

        domainP1Array[0] = domainP1Structure;

        Teuchos::RCP<Teuchos::ParameterList> partitionerParameters = Teuchos::sublist(allParameters, "Mesh Partitioner");
        partitionerParameters->set("Build Edge List", true);
        partitionerParameters->set("Build Surface List", true);

        // Read the given mesh and partition it
        FEDD::MeshPartitioner<SC, LO, GO, NO> p1Partitioner(domainP1Array, partitionerParameters, "P1", dimension);
        p1Partitioner.readAndPartition(15);

        // Paraview export of node and element flags
        domainP1Structure->exportElementFlags();
        domainP1Structure->exportNodeFlags();

        // Build a P2 mesh using the P1 mesh
        domainDiffusion->buildP2ofP1Domain(domainP1Structure);
        domainStructure->buildP2ofP1Domain(domainP1Structure);

        // Set the degrees of freedem per node
        domainStructure->setDofs(dimension);
        domainDiffusion->setDofs(1);

        // Setting the reference configuration for both domains
        domainStructure->setReferenceConfiguration();
        domainDiffusion->setReferenceConfiguration();

        // Define the block system
        Teuchos::RCP<FEDD::SmallMatrix<int>> defTS;
        defTS.reset(new FEDD::SmallMatrix<int>(2));
        (*defTS)[0][0] = 1; // Structure
        (*defTS)[1][1] = 1; // Diffusion

        vec2D_dbl_Type diffusionTensor(dimension, vec_dbl_Type(3));
        double D0 = allParameters->sublist("Parameter Diffusion").get("D0", 1.);
        for (int i = 0; i < dimension; i++)
        {
            diffusionTensor[0][0] = D0;
            diffusionTensor[1][1] = D0;
            diffusionTensor[2][2] = D0;

            if (i > 0)
            {
                diffusionTensor[i][i - 1] = 0;
                diffusionTensor[i - 1][i] = 0;
            }
            else
                diffusionTensor[i][i + 1] = 0;
        }

        // Creating an object of the SCI system
        FEDD::SCI<SC, LO, GO, NO> sci(domainStructure, discretizationType,
                                      domainDiffusion, discretizationType,
                                      diffusionTensor, reactionTerm,
                                      allStructureParameters, allDiffusionParameters,
                                      allParameters, defTS);
        // Printing information
        sci.info();

        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactory(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryDiffusion(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryStructure(new FEDD::BCBuilder<SC, LO, GO, NO>());

        // Getting the surface load parameters
        double pressure = -allParameters->sublist("Parameter").get("Pressure", 0.016);
        double rampTimeStep = allParameters->sublist("Parameter").get("Load Step Size", 0.05);
        double timeRampEnd = allParameters->sublist("Parameter").get("Ramp End Time", 1.0);

        // Setting the surface load parameters
        sci.problemStructureNonLin_->addParemeterRhs(pressure);
        sci.problemStructureNonLin_->addParemeterRhs(rampTimeStep);
        sci.problemStructureNonLin_->addParemeterRhs(timeRampEnd);
        sci.problemStructureNonLin_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter

        // Set load function (the second parameter is the block index)
        sci.problemStructureNonLin_->addRhsFunction(loadFunction, 0);

        // Structure dirichtlet boundary conditions
        bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0
        bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0.5
        bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0, inner ring
        bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0.5 inner ring
        bcFactoryStructure->addBC(zeroDirichlet3D, 6, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0, outer ring
        bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dimension);    // z=0.5, outer ring
        bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dimension); // additional point(s) on outer ring held in x-z direction
        bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dimension); // additional point(s) on outer ring held in y-z direction

        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);

        bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 6, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dimension);

        // Diffusion boundary conditions
        std::vector<double> parameter_vec(1, allParameters->sublist("Parameter").get("Inflow Start Time", 0.));
        bcFactoryDiffusion->addBC(inflowChem, 5, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // Inflow through inner wall
        bcFactoryDiffusion->addBC(inflowChem, 7, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // z=0, inner ring on innter wall
        bcFactoryDiffusion->addBC(inflowChem, 8, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // z=0.5 inner ring

        bcFactory->addBC(inflowChem, 5, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);
        bcFactory->addBC(inflowChem, 7, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);
        bcFactory->addBC(inflowChem, 8, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);

        // Inflow through outer wall
        bcFactoryDiffusion->addBC(inflowChem, 4, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // Inflow through outer wall
        bcFactoryDiffusion->addBC(inflowChem, 6, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // z=0, outer ring on outer wall
        bcFactoryDiffusion->addBC(inflowChem, 9, 0, domainDiffusion, "Dirichlet", 1, parameter_vec); // z=0.5 inner ring

        bcFactory->addBC(inflowChem, 4, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);
        bcFactory->addBC(inflowChem, 6, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);
        bcFactory->addBC(inflowChem, 9, 1, domainDiffusion, "Dirichlet", 1, parameter_vec);

        sci.problemChem_->addBoundaries(bcFactoryDiffusion);

        sci.addBoundaries(bcFactory);

        sci.initializeProblem();

        sci.initializeCE();

        sci.assemble();

        FEDD::DAESolverInTime<SC, LO, GO, NO> daeTimeSolver(allParameters, comm);

        daeTimeSolver.defineTimeStepping(*defTS);

        daeTimeSolver.setProblem(sci);

        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();

        // Restart test: a restarted run (phase 2) must reproduce the uninterrupted run
        // (phase 1): its solution at the final time is compared with the checkpoint
        // phase 1 wrote at that time.
        if (allParameters->sublist("Timestepping Parameter").get("Restart", false))
        {
            double finalTime = allParameters->sublist("Timestepping Parameter").get("Final time compare", 0.0);
            double tolerance = allParameters->sublist("Timestepping Parameter").get("Restart tolerance", 1.e-8);

            bool anyNonzeroReference = false; // comparing identically zero solutions checks nothing
            for (int block = 0; block < 2; block++)
            {
                std::string variable = sci.getVariableName(block);
                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solution = sci.getSolution()->getBlock(block);

                HDF5Import<SC,LO,GO,NO> importer(solution->getMap(), restartFile(allParameters, "Solution" + variable));
                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > uninterrupted = importer.readVariablesHDF5(std::to_string(finalTime));

                Teuchos::RCP<MultiVector<SC,LO,GO,NO> > difference = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>(solution->getMap()));
                difference->update(1., solution, -1., uninterrupted, 0.);

                Teuchos::Array<SC> norm(1), reference(1);
                difference->norm2(norm);
                uninterrupted->norm2(reference);
                double relativeError = reference[0] > 0. ? norm[0] / reference[0] : norm[0];
                anyNonzeroReference = anyNonzeroReference || reference[0] > 0.;

                if (comm->getRank() == 0)
                    std::cout << " Restart test: " << variable << " at t = " << finalTime << ", relative 2-norm of the difference to the uninterrupted run: " << relativeError << " (tolerance " << tolerance << ")" << std::endl;

                if (!(relativeError <= tolerance))
                    testPassed = false;
            }
            if (!anyNonzeroReference)
            {
                if (comm->getRank() == 0)
                    std::cout << " Restart test: the uninterrupted solution is zero at t = " << finalTime << ", the comparison checks nothing" << std::endl;
                testPassed = false;
            }
        }
    }


    FEDD::TimeMonitor_Type::report(std::cout);
    stackedTimer->stop("Structure-chemical interaction");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report((std::cout), comm, options);

    return testPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}

// @brief Reaction Term in the reaction diffusion equation
void reactionTerm(double *x, double *res, double *parameters)
{
    double m = 0.0;
    res[0] = m * x[0];
}

/* The values of parameters are set in feddlib/core/FE/FE_def.hpp assemblySurfaceIntegralExternal()
 * x* for some reason gives the middle point of the surface element
 * parameters[0] is always the time, but all others depend on the order in which they *are added using the addParametersRhs() (which is in Problem_decl.hpp) in the main function
 * Note that the time starts from 0, i.e. the first time step is 0
 */
void loadFunction(double *x, double *res, double *parameters)
{
    res[0] = 0.0;

    double currentTime = parameters[0];
    double pressure = parameters[1];
    double rampTimeStep = parameters[2];
    double timeRampEnd = parameters[3];
    double lambda = 0.0;

    if (currentTime < timeRampEnd)
        lambda = 0.875 * (currentTime + rampTimeStep);
    else
        lambda = 0.875;

    if (parameters[5] == 5) // If the surface flag is 5
        res[0] = pressure * lambda;
}

// @brief Fix all degrees of freedom
void zeroDirichlet3D(double *x, double *res, double t, const double *parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
}

void inflowChem(double *x, double *res, double t, const double *parameters)
{
    if (t >= parameters[0])
        res[0] = 1.;
    else
        res[0] = 0.;
}