# cet_set_compiler_flags(DIAGS CAUTIOUS
#   WERROR
#   NO_UNDEFINED
#   EXTRA_FLAGS  -Wconversion -Wno-sign-conversion
#   )

art_make(
  MODULE_LIBRARIES

                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          art_Persistency_Common canvas
                          art_Persistency_Provenance canvas
                          art_Utilities canvas

                          canvas
                          cetlib
                          cetlib_except

                          larcorealg_Geometry
                          larcore_Geometry_Geometry_service
                          larevt_Filters
                          lardataobj_AnalysisBase
                          lardataobj_RawData
                          lardataobj_RecoBase
                          lardata_Utilities
                          lardata_RecoObjects
                          larpandora_LArPandoraInterface
                          larreco_RecoAlg
                          larsim_MCCheater_BackTrackerService_service
                          larsim_MCCheater_ParticleInventoryService_service
                          larsim_Simulation lardataobj_Simulation

                          nusimdata_SimulationBase
                          ${MF_MESSAGELOGGER}

                          ${FHICLCPP}
                          ${ROOT_GEOM}
                          ${ROOT_XMLIO}
                          ${ROOT_GDML}
                          ${ROOT_BASIC_LIB_LIST}

                          sbndcode_RecoUtils
                          sbndcode_OpDetSim
                          # sbndcode_FlashMatch
        )
install_headers()
install_fhicl()
install_source()
# install_scripts()
cet_enable_asserts()

add_subdirectory(template_generators)
