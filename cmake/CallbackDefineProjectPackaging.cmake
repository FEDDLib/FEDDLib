MACRO(TRIBITS_PROJECT_DEFINE_PACKAGING)

  #MESSAGE("TRIBITS_PROJECT_DEFINE_PACKAGING() called for Trilinos!")
  
  # The CPACK_RESOURCE_FILE_[LICENSE|README] files must end in one of
  # .txt .rtf .html. Copy the pertinant file to the binary directory with
  # a .txt extension. This is only the case with the PackageMaker 
  # generator, but it doesn't hurt to do it for other generators as
  # well.
  TRIBITS_COPY_INSTALLER_RESOURCE(XLib_README
    "${XLib_SOURCE_DIR}/README"
    "${XLib_BINARY_DIR}/README.txt")
  TRIBITS_COPY_INSTALLER_RESOURCE(XLib_LICENSE
    "${XLib_SOURCE_DIR}/LICENSE"
    "${XLib_BINARY_DIR}/LICENSE.txt")
  
    SET(CPACK_PACKAGE_DESCRIPTION "XLib!")
    SET(CPACK_PACKAGE_FILE_NAME "xlib-setup-${XLib_VERSION}")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "XLib ${XLib_VERSION}")
    SET(CPACK_PACKAGE_REGISTRY_KEY "XLib ${XLib_VERSION}")
    SET(CPACK_PACKAGE_NAME "xlib")
    SET(CPACK_PACKAGE_VENDOR "Christian Hochmuth")
    SET(CPACK_PACKAGE_VERSION "${XLib_VERSION}")
    SET(CPACK_RESOURCE_FILE_README "${XLib_README}")
    SET(CPACK_RESOURCE_FILE_LICENSE "${XLib_LICENSE}")
    SET(${PROJECT_NAME}_CPACK_SOURCE_GENERATOR_DEFAULT "TGZ;TBZ2")
    SET(CPACK_SOURCE_FILE_NAME "XLib-source-${XLib_VERSION}")
    SET(CPACK_COMPONENTS_ALL ${XLib_PACKAGES} Unspecified)
  
ENDMACRO()
