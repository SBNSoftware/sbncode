
art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Constants
          SOURCE Constants.cpp
          LIBRARIES
                            ROOT::Core
                            ROOT::Eve
                            ROOT::X3d
                            fhiclcpp::fhiclcpp
                            cetlib::cetlib cetlib_except::cetlib_except
        )


cet_build_plugin( EmptyKaonGen art::tool
                       LIBRARIES
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
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

# these can be removed once icaruscode moves to cetmodules
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )

cet_build_plugin( NuMiKaonGen art::tool
                       LIBRARIES
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
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( TxtFileGen art::tool
                       LIBRARIES
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
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( BNBKaonGen art::tool
                       LIBRARIES
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
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( WeightedRayTraceBox art::tool
                       LIBRARIES
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
                            nusimdata::SimulationBase
                            dk2nu::Tree
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( MixedWeightRayTraceBox art::tool
                       LIBRARIES
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
                            nusimdata::SimulationBase
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)




cet_build_plugin( ReThrowRayTraceBox art::tool
                       LIBRARIES
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
                            nusimdata::SimulationBase
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( RayTraceBox art::tool
                       LIBRARIES
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
                            nusimdata::SimulationBase
                            sbncode_EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

add_subdirectory(Higgs)
add_subdirectory(HNL)

install_headers()
install_fhicl()
install_source()
