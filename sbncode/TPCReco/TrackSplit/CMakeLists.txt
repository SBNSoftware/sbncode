simple_plugin(TrackSplitter module
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  cetlib::cetlib cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
  lardataalg::DetectorInfo
  sbn_LArReco
  sbnobj::Common_Reco
)

simple_plugin(MergedTrackIdentifier module
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  cetlib::cetlib cetlib_except::cetlib_except
  messagefacility::MF_MessageLogger
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larpandora::LArPandoraInterface
  sbn_LArReco
  sbnobj::Common_Reco
)

install_headers()
install_fhicl()
install_source()
