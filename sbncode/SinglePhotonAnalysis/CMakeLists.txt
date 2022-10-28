
art_make_library(
        LIB_LIBRARIES
          larcorealg::Geometry
          larcore::Geometry_Geometry_service
          lardata::Utilities
          larevt::Filters
          larpandoraobj
          lardataobj::RawData
          lardataobj::RecoBase
          lardataobj::MCBase
          larcoreobj::SummaryData
          larsim::MCCheater_BackTrackerService_service
          larsim::MCCheater_ParticleInventoryService_service
          larpandora::LArPandoraInterface
          LArPandoraContent
          nusimdata::SimulationBase
          art::Persistency_Common canvas::canvas
          art::Persistency_Provenance
          art::Utilities
          cetlib::cetlib cetlib_except::cetlib_except
          ROOT::Geom
          ROOT::Core
)

cet_build_plugin(NCDeltaRadiative art::module LIBRARIES
          sbncode::SinglePhotonAnalysis_Libraries
          sbncode::SinglePhotonAnalysis_SEAview
          art_root_io::TFileService_service
          ROOT::Tree
          )

cet_build_plugin(NCRadiativeResonant art::module LIBRARIES
          sbncode::SinglePhotonAnalysis_Libraries
          sbncode::SinglePhotonAnalysis_SEAview
          art_root_io::TFileService_service
          ROOT::Tree
          )

cet_build_plugin(SinglePhoton art::module LIBRARIES
          sbncode::SinglePhotonAnalysis_Libraries
          sbncode::SinglePhotonAnalysis_SEAview
          art_root_io::TFileService_service
          ROOT::Tree
          )


add_subdirectory(SEAview)
add_subdirectory(jobs)
add_subdirectory(Libraries)
add_subdirectory(HelperFunctions)
install_headers()
install_fhicl()
install_source()
