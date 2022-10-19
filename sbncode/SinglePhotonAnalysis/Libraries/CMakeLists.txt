FILE( GLOB local_src_libs *.cxx )

art_make_library( LIBRARY_NAME sbncode_SinglePhotonAnalysis_Libraries
                  SOURCE ${local_src_libs}
                  LIBRARIES
                    sbncode::SinglePhotonAnalysis_HelperFunctions
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
                    art::Persistency_Common canvas::canvas
                    art::Persistency_Provenance
                    art::Utilities
                    cetlib::cetlib cetlib_except::cetlib_except
                    )

