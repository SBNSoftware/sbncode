FILE( GLOB local_src_helpers *.cxx )

art_make_library( LIBRARY_NAME sbncode_SinglePhotonAnalysis_HelperFunctions
                  SOURCE  ${local_src_helpers}
                  LIBRARIES
                    art::Persistency_Common
                    art::Utilities canvas::canvas
                    cetlib::cetlib cetlib_except::cetlib_except
                    larsim::MCCheater_BackTrackerService_service
                    larsim::MCCheater_ParticleInventoryService_service
                    larcorealg::Geometry
                    larpandora::LArPandoraInterface
		    sbnobj::Common_CRT
                    )
