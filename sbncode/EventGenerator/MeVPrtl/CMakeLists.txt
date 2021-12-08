find_ups_product(ifdhc)
find_ups_product(dk2nudata)
find_ups_product(dk2nugenie)

link_directories( $ENV{ROOTSYS}/lib )
link_directories( $ENV{DK2NUDATA_LIB} )
link_directories( $ENV{DK2NUGENIE_LIB} )

add_subdirectory(Tools)
add_subdirectory(config)

simple_plugin( MeVPrtlGen module
              art::Framework_Core
              art::Framework_Services_Registry
              art_root_io::TFileService_service
              ${ART_FRAMEWORK_SERVICES}
				      
				      art_root_io::tfile_support
              art::Framework_Principal
              art::Persistency_Common
              art::Utilities canvas::canvas
              messagefacility::MF_MessageLogger
              
              fhiclcpp::fhiclcpp
              cetlib::cetlib cetlib_except
              ${ROOT_EVE_LIB_LIST}
              ${ROOT_X3d}
              ROOT::Core
              ${ROOT_BASIC_LIB_LIST}
        			ROOT::XMLIO
        			CLHEP::CLHEP
		nusimdata::SimulationBase
              larcorealg_Geometry
        			nurandom::RandomUtils_NuRandomService_service
        			larcoreobj_SummaryData
        			sbnobj_Common_EventGen_MeVPrtl
        			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( MeVPrtlTestRayTrace module
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
				
				art_root_io::tfile_support
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            messagefacility::MF_MessageLogger
                            
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_X3d}
                            ROOT::Core
                            ${ROOT_BASIC_LIB_LIST}
			ROOT::XMLIO
			CLHEP::CLHEP
			nurandom::RandomUtils_NuRandomService_service
			larcoreobj_SummaryData
        		sbnobj_Common_EventGen_MeVPrtl
			sbncode_EventGenerator_MeVPrtl_Constants
)



install_headers()
install_fhicl()
install_source()
