cet_build_plugin(G4InfoReducer art::module
                       LIBRARIES
                           larsim::Simulation
                           larcore::Geometry_Geometry_service
                           lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities canvas::canvas
                           BASENAME_ONLY)

install_fhicl()
install_source()
