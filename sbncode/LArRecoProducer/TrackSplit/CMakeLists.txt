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
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service

  sbn_LArReco_PCA
  sbn_LArReco
  sbncode_LArRecoProducer_Products_dict
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
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larpandora_LArPandoraInterface

  sbn_LArReco
  sbncode_LArRecoProducer_Products_dict
)
install_headers()
install_fhicl()
install_source()
