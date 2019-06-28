
find_package(fccedm)

find_package(ROOT COMPONENTS REQUIRED ROOTDataFrame)

file(GLOB sources *.cc)
file(GLOB headers *.h)

add_library(fcc_analyzers SHARED ${sources} ${headers} )
target_link_libraries(fcc_analyzers podio datamodel ROOT::ROOTDataFrame)

set_target_properties(fcc_analyzers PROPERTIES
  PUBLIC_HEADER "${headers}")

#install(TARGETS utilities
#  # IMPORTANT: Add the podio library to the "export-set"
#  EXPORT fccedmTargets
#  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
#  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/utilities"
#    COMPONENT dev)
