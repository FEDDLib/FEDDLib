# PT: primary tested (high priority TPL)
# ST: secondary tested (medium priority TPL)
# EX: experimental TPL
tribits_repository_define_tpls(
    MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"     PT
    Trilinos        "cmake/TPLs/FindTPLTrilinos.cmake"                  PT
    AceGENInterface "cmake/TPLs/FindTPLAceGENInterface.cmake"           ST
    HDF5            "cmake/TPLs/FindTPLHDF5.cmake"                      PT
    Z               "cmake/TPLs/FindTPLZ.cmake"                         PT
    DL              "cmake/TPLs/FindTPLDL.cmake"                        PT
)
