art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Dalitz
          SOURCE HNLDecayDalitz.cpp
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_X3d}
			sbncode_EventGenerator_MeVPrtl_Constants
        )

simple_plugin( Kaon2HNLFlux tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            messagefacility::MF_MessageLogger
                            
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_X3d}
                            ${ROOT_BASIC_LIB_LIST}
			EGPythia6
			CLHEP::CLHEP
			sbnobj_Common_EventGen_MeVPrtl
			nugen::EventGeneratorBase_GENIE
			nurandom::RandomUtils_NuRandomService_service
			GSL
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( HNLMakeDecay tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            messagefacility::MF_MessageLogger
                            
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_X3d}
                            ${ROOT_BASIC_LIB_LIST}
			CLHEP::CLHEP
			nurandom::RandomUtils_NuRandomService_service
			nusimdata::SimulationBase
			GSL
			sbncode_EventGenerator_MeVPrtl_Dalitz
			sbncode_EventGenerator_MeVPrtl_Constants
)


install_headers()
install_fhicl()
install_source()
