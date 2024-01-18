
cet_build_plugin( Dazzle art::module
  LIBRARIES
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

cet_build_plugin( Razzle art::module
  LIBRARIES
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

cet_build_plugin( Razzled art::module
  LIBRARIES
  ROOT::TMVA
  canvas::canvas
  art_root_io::TFileService_service
  art::Persistency_Provenance
  larsim::Utils
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  sbnobj::Common_Reco
)

install_headers()
install_fhicl()
install_source()

