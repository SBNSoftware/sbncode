INSTALL(FILES services_sbnd_3drift.fcl NumuRecoSBND.fcl NumuRecoICARUS.fcl selectors.fcl
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/fcl
)

add_subdirectory(reprocess)
