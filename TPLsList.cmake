TRIBITS_REPOSITORY_DEFINE_TPLS(
    MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/"     PT
    Trilinos        "cmake/TPLs/"                                       PT
    AceGENInterface "cmake/TPLs/FindTPLAceGENInterface.cmake"           PT
)

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
#
