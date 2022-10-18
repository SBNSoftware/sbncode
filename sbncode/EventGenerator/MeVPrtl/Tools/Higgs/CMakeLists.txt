simple_plugin( MonoEnergyHiggsFlux tool
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
                            sbnobj::Common_EventGen_MeVPrtl
)


simple_plugin( Kaon2HiggsFlux tool
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
                            GSL::gsl
                            sbncode::EventGenerator_MeVPrtl_Constants
                            sbnobj::Common_EventGen_MeVPrtl
)

simple_plugin( HiggsMakeDecay tool
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
