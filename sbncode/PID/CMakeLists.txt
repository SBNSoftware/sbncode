#link_directories(
# $ENV{ROOTSYS}/lib
#)

art_make( MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ${ROOT_EVE_LIB_LIST}
  ${ROOT_X3d}
  ${ROOT_BASIC_LIB_LIST}
  ROOT::Geom
  ROOT::TMVA
  messagefacility::MF_MessageLogger
  sbn_LArReco
  sbnobj_Common_CRT
  sbnobj_Common_Reco
  lardataobj_RecoBase_AssnsDicts_dict
  lardataobj_RecoBase_dict
  lardataobj_RecoBase
  larsim_Utils
  larsim_Simulation
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nusimdata::SimulationBase
  caf_RecoUtils
)

install_headers()
install_fhicl()
install_source()

