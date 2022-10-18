
simple_plugin( TrackCaloSkimmer module
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
                            sbncode::CAFMaker
                            sbnobj::Common_Calibration_dict
                            larevt::SpaceCharge
)

simple_plugin(TrackCaloSkimmerSelectStoppingTrack tool
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

simple_plugin(TrackCaloSkimmerSelectAnode2CathodeTrack tool
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

simple_plugin(TrackCaloSkimmerSelectAllTrack tool
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            larcore::Geometry_Geometry_service
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
)

simple_plugin ( DataSelect module
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

