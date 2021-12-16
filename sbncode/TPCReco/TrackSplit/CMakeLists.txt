simple_plugin(TrackSplitter module
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  messagefacility::MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
  lardataalg_DetectorInfo

  sbn_LArReco
  sbnobj_Common_Reco
)

simple_plugin(MergedTrackIdentifier module
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  messagefacility::MF_MessageLogger
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larpandora_LArPandoraInterface

  sbn_LArReco
  sbnobj_Common_Reco
)

install_headers()
install_fhicl()
install_source()
