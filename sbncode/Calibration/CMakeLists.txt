
cet_build_plugin( TrackCaloSkimmer art::module
                       LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            sbnobj::Common_Reco
                            larcorealg::GeoAlgo
                            lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
                            lardataalg::DetectorInfo
                            ROOT::Minuit
                            sbncode_CAFMaker
                            sbnobj::Common_Calibration_dict
                            larevt::SpaceCharge
)

cet_build_plugin(TrackCaloSkimmerSelectStoppingTrack art::tool
                      LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
                            larcore::Geometry_Geometry_service
                            larcorealg::Geometry
                            larcorealg::GeoAlgo
                            lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
                            lardataalg::DetectorInfo
                            larevt::SpaceCharge
                            sbnobj::Common_Calibration_dict
)

cet_build_plugin(TrackCaloSkimmerSelectAnode2CathodeTrack art::tool
                      LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
                            sbnobj::Common_Calibration_dict
                            larevt::SpaceCharge
)

cet_build_plugin(TrackCaloSkimmerSelectAllTrack art::tool
                      LIBRARIES
			    sbnobj::Common_Calibration_dict
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            larcore::Geometry_Geometry_service
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
)

cet_build_plugin ( DataSelect art::module
           LIBRARIES
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               cetlib::cetlib cetlib_except::cetlib_except
               ROOT::X3d
               sbnobj::Common_Calibration_dict
               larevt::SpaceCharge
               lardataobj::RecoBase
)

add_subdirectory(fcl)

install_headers()
install_fhicl()
install_source()

