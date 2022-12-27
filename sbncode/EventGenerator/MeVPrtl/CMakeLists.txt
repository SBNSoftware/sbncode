add_subdirectory(Tools)
add_subdirectory(config)

cet_build_plugin( MeVPrtlGen art::module
          LIBRARIES
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
              sbncode_EventGenerator_MeVPrtl_Constants
)

cet_build_plugin( MeVPrtlTestRayTrace art::module
                    LIBRARIES
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
                        nusimdata::SimulationBase
                        larcoreobj::SummaryData
                        sbnobj::Common_EventGen_MeVPrtl
                        sbncode_EventGenerator_MeVPrtl_Constants
                        lardata::Utilities
)



install_headers()
install_fhicl()
install_source()
