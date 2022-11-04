art_make_library( LIBRARY_NAME sbncode::EventGenerator_MeVPrtl_Dalitz
          SOURCE HNLDecayDalitz.cpp
          LIBRARIES
                        sbncode::EventGenerator_MeVPrtl_Constants
                        sbnobj::Common_EventGen_MeVPrtl
        )


cet_build_plugin( Kaon2HNLFlux art::tool
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::EGPythia6
                            CLHEP::CLHEP
                            sbnobj::Common_EventGen_MeVPrtl
                            nugen::EventGeneratorBase_GENIE
                            nusimdata::SimulationBase
                            nurandom::RandomUtils_NuRandomService_service
                            GSL::gsl
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( HNLMakeDecay art::tool
                        LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            CLHEP::CLHEP
                            nurandom::RandomUtils_NuRandomService_service
                            nusimdata::SimulationBase
                            GSL::gsl
                            larcorealg::Geometry
                            sbncode::EventGenerator_MeVPrtl_Dalitz
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)


install_headers()
install_fhicl()
install_source()
