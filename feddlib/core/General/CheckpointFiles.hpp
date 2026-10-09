#ifndef CHECKPOINTFILES_hpp
#define CHECKPOINTFILES_hpp

#include "feddlib/core/FEDDCore.hpp"

#include <algorithm>
#include <string>
#include <vector>

/*!
 Paths of the HDF5 checkpoint files of time-dependent problems.

 With "Checkpointing" (sublist "Timestepping Parameter") a run writes, at every
 checkpoint time t of "Checkpoint Times" (a list, e.g.
 <Parameter name="Checkpoint Times" type="Array(double)" value="{0.05, 0.1}"/>), the files Solution<variable>, SolutionNewmark<variable>,
 ds_Velocity, ds_Acceleration, Rhs<variable> and Solutiond_f, whose variables
 are named by t, and the element history History<t>. With "Restart" a run reads
 them back at the time "Time step".

 The files are written to "Checkpoint directory" and read from
 "Restart directory" (both in "Timestepping Parameter"; by default the run
 directory), so a restart can read the checkpoints of another run in place.
 */

namespace FEDD {

inline std::string joinPath(const std::string& directory, const std::string& file)
{
    if (directory.empty() || directory == ".")
        return file;
    return directory.back() == '/' ? directory + file : directory + "/" + file;
}

/// The checkpoint times ("Checkpoint Times" in "Timestepping Parameter"), in increasing order.
inline std::vector<double> checkpointTimes(const ParameterListPtr_Type& parameters)
{
    Teuchos::ParameterList& timestepping = parameters->sublist("Timestepping Parameter");
    TEUCHOS_TEST_FOR_EXCEPTION(timestepping.isParameter("Number Checkpoints") || timestepping.isSublist("Checkpoints"), std::runtime_error,
        "The checkpoints are a list now: replace \"Number Checkpoints\" and the sublist \"Checkpoints\" with "
        "<Parameter name=\"Checkpoint Times\" type=\"Array(double)\" value=\"{t1, t2, ...}\"/> in \"Timestepping Parameter\".");
    Teuchos::Array<double> times = timestepping.get("Checkpoint Times", Teuchos::Array<double>());
    std::vector<double> sorted(times.begin(), times.end());
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

/// Path of the checkpoint file 'file' a run writes.
inline std::string checkpointFile(const ParameterListPtr_Type& parameters, const std::string& file)
{
    return joinPath(parameters->sublist("Timestepping Parameter").get("Checkpoint directory", std::string("")), file);
}

/// Path of the checkpoint file 'file' a restart reads.
inline std::string restartFile(const ParameterListPtr_Type& parameters, const std::string& file)
{
    return joinPath(parameters->sublist("Timestepping Parameter").get("Restart directory", std::string("")), file);
}

}
#endif
