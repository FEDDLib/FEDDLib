tribits_repository_define_tpls(
    MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"     PT
    Trilinos        "cmake/TPLs/FindTPLTrilinos.cmake"                  PT
    AceGENInterface "cmake/TPLs/FindTPLAceGENInterface.cmake"           PT
)
