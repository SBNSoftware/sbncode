art_make(    LIBRARY_NAME sbncode_GeoWrappers
             LIB_LIBRARIES larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::Simulation lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           larevt::Filters
                           lardataobj::RawData
                           lardataobj::RecoBase
                           lardata::RecoObjects  
                           larpandora::LArPandoraInterface
                           nusimdata::SimulationBase
                           cetlib::cetlib cetlib_except::cetlib_except
                           art::Persistency_Common canvas::canvas
                           art::Persistency_Provenance
                           art::Utilities
        )

install_headers()
install_source()
