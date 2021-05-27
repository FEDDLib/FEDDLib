MACRO(TRIBITS_REPOSITORY_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_REPOSITORY_DEFINE_PACKAGING() called for XLib!")
 
  # We need to make sure that these excludes only apply to XLib, not the global
  # project.
  SET(XLib_SOURCE_EXCLUDE_DIR ${XLib_SOURCE_DIR})
  #PRINT_VAR(Trilinos_SOURCE_EXCLUDE_DIR)

    SET(CPACK_SOURCE_IGNORE_FILES
    ${CPACK_SOURCE_IGNORE_FILES}
    /.git/
    ".gitignore"
    )

  APPEND_SET(TRIBITS_CPACK_PACKAGES_TO_NOT_IGNORE TriBITS)
  
ENDMACRO()
