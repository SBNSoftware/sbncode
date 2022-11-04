cet_build_plugin( MonoEnergyHiggsFlux art::tool
                      LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            CLHEP::CLHEP
                            nugen::EventGeneratorBase_GENIE
                            nurandom::RandomUtils_NuRandomService_service
                            sbncode::EventGenerator_MeVPrtl_Constants
                            nusimdata::SimulationBase
                            sbnobj::Common_EventGen_MeVPrtl
)


cet_build_plugin( Kaon2HiggsFlux art::tool
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
                            nurandom::RandomUtils_NuRandomService_service
                            nusimdata::SimulationBase
                            GSL::gsl
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

cet_build_plugin( HiggsMakeDecay art::tool
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
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
)

install_headers()
install_fhicl()
install_source()
