simple_plugin(TrackSplitter module
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
  ${ART_FRAMEWORK_PRINCIPAL}
  art_Persistency_Common
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_BASIC_LIB_LIST}
  MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
  lardataalg_DetectorInfo

  sbn_LArReco
  sbnobj_Common_Reco
)

simple_plugin(MergedTrackIdentifier module
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
  ${ART_FRAMEWORK_PRINCIPAL}
  art_Persistency_Common
  art_Utilities canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_BASIC_LIB_LIST}
  MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larpandora_LArPandoraInterface

  sbn_LArReco
  sbnobj_Common_Reco
)

install_headers()
install_fhicl()
install_source()
