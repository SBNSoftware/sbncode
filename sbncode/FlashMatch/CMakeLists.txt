#cet_set_compiler_flags(DIAGS CAUTIOUS
#  WERROR
#  NO_UNDEFINED
#  EXTRA_FLAGS  -Wall -Wextra  -pedantic  #-Wconversion
#  )

art_make( 
  MODULE_LIBRARIES
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support ROOT::Core ROOT::Tree
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities canvas::canvas
  canvas::canvas
  cetlib::cetlib
  cetlib_except::cetlib_except
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata::Utilities
  lardataalg::DetectorInfo
  lardataobj::AnalysisBase
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::Simulation
  larpandora::LArPandoraInterface
  larsim::Utils
  sbnobj::Common_Reco
  nusimdata::SimulationBase
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  )

install_headers()
install_fhicl()
install_source()
cet_enable_asserts()

add_subdirectory(template_generators)
