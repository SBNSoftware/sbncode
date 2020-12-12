cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS  -Wall -Wextra  -pedantic  #-Wconversion
  )

art_make(
  MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Persistency_Provenance
  art_Utilities canvas

  canvas
  cetlib
  cetlib_except

  larcore_Geometry_Geometry_service
  larcorealg_Geometry

  lardata_Utilities
  lardataobj_AnalysisBase
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_Simulation

  larpandora_LArPandoraInterface

  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${ROOT_GEOM}
  ${ROOT_XMLIO}
  ${ROOT_GDML}
  ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()
cet_enable_asserts()

add_subdirectory(template_generators)
