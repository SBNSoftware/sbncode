
art_make( MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Persistency_Common
  art::Utilities canvas::canvas
  cetlib::cetlib cetlib_except::cetlib_except
  ROOT::TMVA
  messagefacility::MF_MessageLogger
  sbn_LArReco
  sbnobj::Common_CRT
  sbnobj::Common_Reco
  lardataobj::RecoBase_AssnsDicts_dict
  lardataobj::RecoBase_dict
  lardataobj::RecoBase
  larsim::Utils
  larsim::Simulation
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  nusimdata::SimulationBase
  caf_RecoUtils
)

install_headers()
install_fhicl()
install_source()

