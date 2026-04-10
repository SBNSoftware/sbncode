
art_make_library( LIBRARY_NAME caf_RecoUtils
                  SOURCE  RecoUtils.cc
                  LIBRARIES
                    art::Framework_Core
                    art::Framework_Services_Registry
                    art::Framework_Principal
                    art::Persistency_Common
                    art::Utilities canvas::canvas
                    cetlib::cetlib cetlib_except::cetlib_except
                    larsim::MCCheater_BackTrackerService_service
                    larsim::MCCheater_ParticleInventoryService_service
                    larcorealg::Geometry
                )
