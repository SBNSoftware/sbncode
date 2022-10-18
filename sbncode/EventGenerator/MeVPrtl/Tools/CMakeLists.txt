
art_make_library( LIBRARY_NAME sbncode::EventGenerator_MeVPrtl_Constants
          SOURCE Constants.cpp
          LIBRARIES
                            ROOT::Core
                            ROOT::Eve
                            ROOT::X3d
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
        )


simple_plugin( EmptyKaonGen tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            dk2nu::Tree
                            CLHEP::CLHEP
                            nugen::EventGeneratorBase_GENIE
                            nurandom::RandomUtils_NuRandomService_service
                            nusimdata::SimulationBase
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)


simple_plugin( NuMiKaonGen tool
                            ifdh_art::IFDH_service
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            CLHEP::CLHEP
                            ifdh::ifdh
                            dk2nu::Tree
                            nurandom::RandomUtils_NuRandomService_service
                            nusimdata::SimulationBase
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

simple_plugin( WeightedRayTraceBox tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            CLHEP::CLHEP
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            nurandom::RandomUtils_NuRandomService_service
                            dk2nu::Tree
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

simple_plugin( MixedWeightRayTraceBox tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            dk2nu::Tree
                            CLHEP::CLHEP
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            nurandom::RandomUtils_NuRandomService_service
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)




simple_plugin( ReThrowRayTraceBox tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            CLHEP::CLHEP
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            nurandom::RandomUtils_NuRandomService_service
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

simple_plugin( RayTraceBox tool
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::Eve
                            ROOT::X3d
                            ROOT::Core
                            CLHEP::CLHEP
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            nurandom::RandomUtils_NuRandomService_service
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

add_subdirectory(Higgs)
add_subdirectory(HNL)

install_headers()
install_fhicl()
install_source()
