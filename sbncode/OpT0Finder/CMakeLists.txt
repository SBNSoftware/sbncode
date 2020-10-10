# include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Turning off WERROR only for this package, because it needs to
# call OpFlastScintillation which has an unused variable.
# Once the semi-analytic photon propogation will be moved
# to its own service, and we won't need to call OpFlastScintillation
# in larsim, we will be able to remobe this. -- Marco Del Tutto, Oct 2020
#cet_set_compiler_flags(DIAGS CAUTIOUS
#  # WERROR
#  NO_UNDEFINED
#  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
#)
include_directories($ENV{GEANT4_FQ_DIR}/include)

#art_make(
#    LIB_LIBRARIES
#        sbndcode_OpT0Finder_flashmatch_Base
#        sbndcode_OpT0Finder_flashmatch_Algorithms
#        sbndcode_Geometry
#        sbndcode_OpDetSim
#        larcorealg_Geometry
#        larcore_Geometry_Geometry_service
#        lardata_Utilities
#        lardataobj_RecoBase
#        larsim_Simulation
#        # larsim_LegacyLArG4
#        lardataobj_AnalysisBase
#        lardataobj_Simulation
#        nusimdata_SimulationBase
#        ${LARDATA_LIB}
#        canvas
#        cetlib
#        cetlib_except
#        ${ART_FRAMEWORK_ART}
#        ${ART_FRAMEWORK_CORE}
#        ${ART_FRAMEWORK_PRINCIPAL}
#        ${ART_FRAMEWORK_UTILITIES}
#        ${ART_FRAMEWORK_BASIC}
#        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
#        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#        art_Persistency_Common
#        art_Persistency_Provenance
#        art_Utilities
#        ${MF_MESSAGELOGGER}
#        ${MF_UTILITIES}
#        ${FHICLCPP}
#        ${ROOT_GEOM}
#        ${ROOT_XMLIO}
#        ${ROOT_GDML}
#        ${ROOT_BASIC_LIB_LIST}
#    MODULE_LIBRARIES
#        sbndcode_OpT0Finder_flashmatch_Base
#        sbndcode_OpT0Finder_flashmatch_Algorithms
#        sbndcode_OpDetReco_OpFlash_FlashFinder
#        sbndcode_OpDetSim
#        sbndcode_Geometry
#        larcorealg_Geometry
#        larcore_Geometry_Geometry_service
#        lardata_Utilities
#        lardataobj_AnalysisBase
#        lardataobj_RecoBase
#        larsim_Simulation
#        # larsim_LegacyLArG4
#        lardataobj_Simulation
#        larpandora_LArPandoraInterface
#        nusimdata_SimulationBase
#        ${LARDATA_LIB}
#        canvas
#        cetlib
#        cetlib_except
#        ${ART_FRAMEWORK_ART}
#        ${ART_FRAMEWORK_CORE}
#        ${ART_FRAMEWORK_PRINCIPAL}
#        ${ART_FRAMEWORK_BASIC}
#        ${ART_FRAMEWORK_SERVICES_REGISTRY}
#        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
#        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#        art_Persistency_Common
#        art_Persistency_Provenance
#        art_Utilities
#        ${MF_MESSAGELOGGER}
#        ${MF_UTILITIES}
#        ${FHICLCPP}
#        ${ROOT_GEOM}
#        ${ROOT_XMLIO}
#        ${ROOT_GDML}
#        ${ROOT_BASIC_LIB_LIST}
#)

#install_headers()
#install_fhicl()
#install_source()

add_subdirectory(flashmatch)
