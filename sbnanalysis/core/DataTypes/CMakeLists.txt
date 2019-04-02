include_directories(${CMAKE_CURRENT_SOURCE_DIR})

cet_make(
  LIBRARY_NAME sbncode_sbnanalysis_core_DataTypes
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_CORE}
    ${ROOT_EG}
)

install_headers()
install_fhicl()
install_source()

