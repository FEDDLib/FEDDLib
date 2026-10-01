// Hollow cylinder (Ri=1.0, Ro=1.3, L=2.0) with the SCI SMC element
// (SCI_SMC_Active_Growth_Reorientation, see materialParameters.xml), BC scheme
// (bottom/top faces axially fixed, 3 circumferential pins for static
// determinacy, ramped inner-wall pressure 0->20 over 10 steps), all
// active/growth/reorientation bools off (purely passive, load-driven
// response only). Adapted from feddlib/problems/examples/arteries/
// simpleTest/main.cpp -- only the mesh, BC flags, and pressure ramp are
// changed; the overall SCI/BCBuilder/DAESolverInTime driving pattern is
// identical.
//
// Mesh face/vertex flags (see meshes/ccb_hollow_cylinder/
// hollow_cylinder_p1.mesh):
//   15 - volume (single material)
//    2 - bottom face (z=0)                -> Dirichlet_Z
//    3 - top face (z=2)                   -> Dirichlet_Z
//    4 - outer face (r=1.3, no pins)       -> traction-free
//    5 - inner face (r=1.0)                -> ramped pressure
//   13 - pin vertex at theta=0   (outer)   -> Dirichlet_Y (tangential dir)
//   14 - pin vertex at theta=90  (outer)   -> Dirichlet_X (tangential dir)
//   16 - pin vertex at theta=180 (outer)   -> Dirichlet_Y (tangential dir)

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

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

void reactionTerm(double *, double *, double *);
void loadFunction(double *, double *, double *);
void zeroDirichlet3D(double *, double *, double, const double *);

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
    std::string materialParametersXML = "materialParameters.xml";
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

    Teuchos::RCP<Teuchos::StackedTimer> stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("CCB hollow cylinder", true));
    bool verbose(comm->getRank() == 0);

    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
    {
        Teuchos::RCP<Teuchos::ParameterList> simulationParameters = Teuchos::getParametersFromXmlFile(simulationParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> materialParameters = Teuchos::getParametersFromXmlFile(materialParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> solverParameters = Teuchos::getParametersFromXmlFile(solverParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> structurePreconditionerParameters = Teuchos::getParametersFromXmlFile(structurePreconditionerParametersXML);
        Teuchos::RCP<Teuchos::ParameterList> chemistryPreconditionerParamerters = Teuchos::getParametersFromXmlFile(chemistryPreconditionerParametersXML);

        int dimension = simulationParameters->sublist("Simulation Parameters").get("Dimension", 3);
        std::string discretizationType = simulationParameters->sublist("Simulation Parameters").get("Discretization", "P2");

        Teuchos::RCP<Teuchos::ParameterList> allParameters = Teuchos::rcp(new Teuchos::ParameterList(*simulationParameters));
        allParameters->sublist("Parameter").set("Chemistry Explicit", false);

        Teuchos::RCP<Teuchos::ParameterList> preconditionerParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));

        allParameters->setParameters(*materialParameters);
        allParameters->setParameters(*structurePreconditionerParameters);
        allParameters->setParameters(*solverParameters);

        Teuchos::RCP<Teuchos::ParameterList> allDiffusionParameters = Teuchos::rcp(new Teuchos::ParameterList(*chemistryPreconditionerParamerters));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Chem"));
        Teuchos::sublist(allDiffusionParameters, "Parameter")->setParameters(simulationParameters->sublist("Simulation Parameters"));
        allDiffusionParameters->setParameters(*solverParameters);
        allDiffusionParameters->setParameters(*chemistryPreconditionerParamerters);

        Teuchos::RCP<Teuchos::ParameterList> allStructureParameters = Teuchos::rcp(new Teuchos::ParameterList(*structurePreconditionerParameters));
        Teuchos::sublist(allStructureParameters, "Parameter")->setParameters(simulationParameters->sublist("Parameter Solid"));
        allStructureParameters->setParameters(*materialParameters);
        allStructureParameters->setParameters(*solverParameters);

        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Diffusion;
        Teuchos::RCP<FEDD::Domain<SC, LO, GO, NO>> domainP1Structure;
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

        // Read the hollow-cylinder P1 mesh, no unit conversion.
        FEDD::MeshPartitioner<SC, LO, GO, NO> p1Partitioner(domainP1Array, partitionerParameters, "P1", dimension);
        p1Partitioner.readAndPartition(15);

        domainP1Structure->exportElementFlags();
        domainP1Structure->exportNodeFlags();

        domainDiffusion->buildP2ofP1Domain(domainP1Structure);
        domainStructure->buildP2ofP1Domain(domainP1Structure);

        domainStructure->setDofs(dimension);
        domainDiffusion->setDofs(1);

        domainStructure->setReferenceConfiguration();
        domainDiffusion->setReferenceConfiguration();

        Teuchos::RCP<FEDD::SmallMatrix<int>> defTS;
        defTS.reset(new FEDD::SmallMatrix<int>(2));
        (*defTS)[0][0] = 1;
        (*defTS)[1][1] = 1;

        FEDD::vec2D_dbl_Type diffusionTensor(dimension, FEDD::vec_dbl_Type(3));
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

        FEDD::SCI<SC, LO, GO, NO> sci(domainStructure, discretizationType,
                                      domainDiffusion, discretizationType,
                                      diffusionTensor, reactionTerm,
                                      allStructureParameters, allDiffusionParameters,
                                      allParameters, defTS);
        sci.info();

        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactory(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryDiffusion(new FEDD::BCBuilder<SC, LO, GO, NO>());
        Teuchos::RCP<FEDD::BCBuilder<SC, LO, GO, NO>> bcFactoryStructure(new FEDD::BCBuilder<SC, LO, GO, NO>());

        double rampTimeStep = allParameters->sublist("Parameter").get("Load Step Size", 0.1);
        double timeRampEnd = allParameters->sublist("Parameter").get("Ramp End Time", 1.0);

        // Plain pressure units, no mmHg conversion -- "Max Pressure mmHg"
        // is used directly as the target pressure (ramped 0 -> 20).
        // The negative sign is
        // kept from simpleTest's convention (outward/inflating internal
        // pressure).
        double maxPressureMmHg = allParameters->sublist("Parameter").get("Max Pressure mmHg", 20.0);
        double targetPressureMmHg = allParameters->sublist("Parameter").get("Target Pressure mmHg", 20.0);
        double pressureReductionStartTime = allParameters->sublist("Parameter").get("Pressure Reduction Start Time", 1000000.0);
        double pressureReductionEndTime = allParameters->sublist("Parameter").get("Pressure Reduction End Time", 1000000.0);
        double pressureReductionAmountMmHg = allParameters->sublist("Parameter").get("Pressure Reduction Amount mmHg", 0.0);

        double pressure = -maxPressureMmHg;
        double initialLambda = targetPressureMmHg / maxPressureMmHg;
        double lambdaReduction = pressureReductionAmountMmHg / maxPressureMmHg;

        sci.problemStructureNonLin_->addParemeterRhs(pressure);
        sci.problemStructureNonLin_->addParemeterRhs(rampTimeStep);
        sci.problemStructureNonLin_->addParemeterRhs(timeRampEnd);
        sci.problemStructureNonLin_->addParemeterRhs(initialLambda);
        sci.problemStructureNonLin_->addParemeterRhs(pressureReductionStartTime);
        sci.problemStructureNonLin_->addParemeterRhs(pressureReductionEndTime);
        sci.problemStructureNonLin_->addParemeterRhs(lambdaReduction);
        sci.problemStructureNonLin_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter

        sci.problemStructureNonLin_->addRhsFunction(loadFunction, 0);

        // Structure Dirichlet BCs: bottom/top faces axially fixed, 3
        // circumferential pins (single-component, aligned with global
        // X/Y since the pins sit at theta=0/90/180) for static
        // determinacy.
        bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Y", dimension);
        bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_X", dimension);
        bcFactoryStructure->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet_Y", dimension);

        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);

        bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dimension);
        bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Y", dimension);
        bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_X", dimension);
        bcFactory->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet_Y", dimension);

        // No diffusion Dirichlet BCs
        // (concentration is left with natural/zero-flux boundaries
        // everywhere -- the model is passive, so this dof only matters
        // through Kuc/Kcu/Kcc's own internal dynamics).
        sci.problemChem_->addBoundaries(bcFactoryDiffusion);

        sci.addBoundaries(bcFactory);

        sci.initializeProblem();
        sci.initializeCE();

        FEDD::DAESolverInTime<SC, LO, GO, NO> daeTimeSolver(allParameters, comm);
        daeTimeSolver.defineTimeStepping(*defTS);
        daeTimeSolver.setProblem(sci);
        sci.assemble();
        daeTimeSolver.setupTimeStepping();
        daeTimeSolver.advanceInTime();
    }
    FEDD::TimeMonitor_Type::report(std::cout);
    stackedTimer->stop("CCB hollow cylinder");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report((std::cout), comm, options);

    return (EXIT_SUCCESS);
}

// No separate reaction term outside the AceGen kernel's own internal Rc
// computation (m=0 is a no-op).
void reactionTerm(double *x, double *res, double *parameters)
{
    double m = 0.0;
    res[0] = m * x[0];
}

/* Parameter order (see Problem_decl.hpp addParemeterRhs / FE_def.hpp
 * assemblySurfaceIntegralExternal):
 * parameters[0]: current time
 * parameters[1]: pressure
 * parameters[2]: rampTimeStep
 * parameters[3]: timeRampEnd
 * parameters[4]: initialLambda
 * parameters[5]: pressureReductionStartTime
 * parameters[6]: pressureReductionEndTime
 * parameters[7]: lambdaReduction
 * parameters[8]: degree of the load function in space (0)
 * parameters[9]: surface flag
 */
void loadFunction(double *x, double *res, double *parameters)
{
    res[0] = 0.0;

    double currentTime = parameters[0];
    double pressure = parameters[1];
    double timeRampEnd = parameters[3];
    double initialLambda = parameters[4];
    double pressureReductionStartTime = parameters[5];
    double pressureReductionEndTime = parameters[6];
    double lambdaReduction = parameters[7];
    double surfaceFlag = parameters[9];

    double lambda = 0.0;
    double currentLambdaReduction = 0.0;

    // currentTime is already t_{n+1} (DAESolverInTime advances time before the solve), so
    // adding rampTimeStep here would apply every load one step early.
    if (currentTime < timeRampEnd)
        lambda = initialLambda * currentTime / timeRampEnd;
    else
        lambda = initialLambda;

    if (currentTime >= pressureReductionStartTime && currentTime <= pressureReductionEndTime) {
        double reductionProgress = (currentTime - pressureReductionStartTime) /
                                 (pressureReductionEndTime - pressureReductionStartTime);
        currentLambdaReduction = lambdaReduction * reductionProgress;
    } else if (currentTime > pressureReductionEndTime) {
        currentLambdaReduction = lambdaReduction;
    }

    lambda = lambda - currentLambdaReduction;
    if (lambda < 0.0) lambda = 0.0;

    if (surfaceFlag == 5) // inner wall
        res[0] = pressure * lambda;
}

// Fix all degrees of freedom (BCBuilder picks out only the component(s)
// its BC type string names, e.g. "Dirichlet_Z" only uses res[2]).
void zeroDirichlet3D(double *x, double *res, double t, const double *parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
}
