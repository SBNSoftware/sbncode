
art_make_library( LIBRARY_NAME caf_RecoUtils
                  SOURCE  RecoUtils.cc
                  LIBRARIES
                  PRIVATE
                    art::Framework_Services_Registry
                    larsim::MCCheater_BackTrackerService_service
                    larsim::MCCheater_ParticleInventoryService_service
                  PUBLIC
                    lardataalg::DetectorInfo
                    lardataobj::RecoBase
                    canvas::canvas
                )
install_headers()
