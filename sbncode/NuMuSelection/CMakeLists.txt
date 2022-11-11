cet_build_plugin( MuPVertexStudy art::module
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            caf_RecoUtils
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            lardataalg::DetectorInfo
                            lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                            larsim::MCCheater_BackTrackerService_service
                            nusimdata::SimulationBase
                            larsim::MCCheater_ParticleInventoryService_service
                            larcorealg::GeoAlgo
                            larreco::RecoAlg
)
cet_build_plugin( MuonS2NStudy art::module
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            caf_RecoUtils
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            lardataalg::DetectorInfo
                            lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                            larsim::MCCheater_BackTrackerService_service
                            nusimdata::SimulationBase
                            larsim::MCCheater_ParticleInventoryService_service
                            larcorealg::GeoAlgo
                            larreco::RecoAlg
)
cet_build_plugin( NuMuEfficiencyStudy art::module
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            caf_RecoUtils
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            lardataalg::DetectorInfo
                            lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                            larsim::MCCheater_BackTrackerService_service
                            nusimdata::SimulationBase
                            larsim::MCCheater_ParticleInventoryService_service
                            larcorealg::GeoAlgo
                            larreco::RecoAlg
)

install_headers()
install_fhicl()
install_source()
