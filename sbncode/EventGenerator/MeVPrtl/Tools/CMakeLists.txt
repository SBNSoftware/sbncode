link_directories( $ENV{ROOTSYS}/lib )

art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Constants
          SOURCE Constants.cpp
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_X3d}
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except
        )


simple_plugin( EmptyKaonGen tool
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
			sbncode_EventGenerator_MeVPrtl_Constants
)

cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
##include_directories( $ENV{IFDHC_FQ_DIR}/inc )
##include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )


simple_plugin( NuMiKaonGen tool
                            ifdh_art::IFDH_service
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
			ifdh
			dk2nuTree
			nurandom::RandomUtils_NuRandomService_service
			nusimdata::SimulationBase
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( BNBKaonGen tool
                            ifdh_art::IFDH_service
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
			ifdh
			dk2nuTree
			nurandom::RandomUtils_NuRandomService_service
			nusimdata::SimulationBase
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( WeightedRayTraceBox tool
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
			larcorealg_Geometry
			nurandom::RandomUtils_NuRandomService_service
			dk2nuGenie
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( MixedWeightRayTraceBox tool
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
			larcorealg_Geometry
			nurandom::RandomUtils_NuRandomService_service
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)




simple_plugin( ReThrowRayTraceBox tool
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
			larcorealg_Geometry
			nurandom::RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( RayTraceBox tool
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
			larcorealg_Geometry
			nurandom::RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

add_subdirectory(Higgs)
add_subdirectory(HNL)

install_headers()
install_fhicl()
install_source()
