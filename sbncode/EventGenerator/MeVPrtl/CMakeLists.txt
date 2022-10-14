find_package(ifdhc) 
find_package(dk2nudata) 
find_package(dk2nugenie) 

add_subdirectory(Tools)
add_subdirectory(config)

simple_plugin( MeVPrtlGen module
              art::Framework_Core
              art::Framework_Services_Registry
              art_root_io::TFileService_service
				      art_root_io::TFileService_service
				      art_root_io::tfile_support
              art::Framework_Principal
              art::Persistency_Common
              art::Utilities canvas::canvas
              cetlib::cetlib cetlib_except::cetlib_except
        			CLHEP::CLHEP
              nugen::EventGeneratorBase_GENIE
		nusimdata::SimulationBase
              larcorealg::Geometry
              lardata::Utilities
        			nurandom::RandomUtils_NuRandomService_service
        			larcoreobj::SummaryData
        			sbnobj::Common_EventGen_MeVPrtl
        			sbncode::EventGenerator_MeVPrtl_Constants
)

simple_plugin( MeVPrtlTestRayTrace module
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
				art_root_io::TFileService_service
				art_root_io::tfile_support
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
			CLHEP::CLHEP
                        nugen::EventGeneratorBase_GENIE
			nurandom::RandomUtils_NuRandomService_service
			larcoreobj::SummaryData
        		sbnobj::Common_EventGen_MeVPrtl
			sbncode::EventGenerator_MeVPrtl_Constants
                        lardata::Utilities
)



install_headers()
install_fhicl()
install_source()
