find_package(larpandora)

art_make_library(
    LIBRARY_NAME sbncode_SinglePhotonAnalysis_SEAview
    LIBRARIES
      larcorealg::Geometry
      larcore::Geometry_Geometry_service
      lardata::Utilities
      larevt::Filters
      lardataobj::RawData
      lardataobj::RecoBase
      lardataobj::MCBase
      larcoreobj::SummaryData
      larsim::MCCheater_BackTrackerService_service
      larsim::MCCheater_ParticleInventoryService_service
      larpandora::LArPandoraInterface
      larpandoracontent::LArPandoraContent 
      nusimdata::SimulationBase
      nurandom::RandomUtils_NuRandomService_service
      art::Persistency_Common canvas::canvas
      art::Persistency_Provenance
      art::Utilities
      cetlib::cetlib cetlib_except::cetlib_except
      ROOT::Core
      ROOT::Geom 
      ROOT::Gpad
)
