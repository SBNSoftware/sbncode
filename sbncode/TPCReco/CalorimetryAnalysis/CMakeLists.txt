cet_build_plugin( CalorimetryAnalysis art::module
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::EG
                            ROOT::X3d
                            ROOT::Core
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            caf_RecoUtils
                            larreco::RecoAlg
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            larsim::MCCheater_BackTrackerService_service
                            nusimdata::SimulationBase
                            larsim::MCCheater_ParticleInventoryService_service
                            sbnobj::Common_Reco
                            larcorealg::GeoAlgo
                            lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
                            lardataalg::DetectorInfo
                            larevt::SpaceCharge
)


install_headers()
install_fhicl()
install_source()

