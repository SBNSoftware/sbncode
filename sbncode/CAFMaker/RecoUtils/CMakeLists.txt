
art_make_library( LIBRARY_NAME caf_RecoUtils
                  SOURCE  RecoUtils.cc
                  LIBRARIES 
		    art::Framework_Core
		    art::Framework_Services_Registry
		    ${ART_FRAMEWORK_SERVICES}
		    art::Framework_Principal
		    art::Persistency_Common
		    art::Utilities canvas::canvas
                    ${ROOT_BASIC_LIB_LIST}
		    messagefacility::MF_MessageLogger
		    
		    fhiclcpp::fhiclcpp
                    cetlib::cetlib cetlib_except
		    larsim_MCCheater_BackTrackerService_service
		    larsim_MCCheater_ParticleInventoryService_service
                    larcorealg_Geometry
                )
