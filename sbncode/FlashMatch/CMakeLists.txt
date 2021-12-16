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
  art_root_io::tfile_support ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities canvas::canvas
  canvas::canvas
  cetlib::cetlib
  cetlib_except::cetlib_except
  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata_Utilities
  lardataobj_AnalysisBase
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_Simulation
  larpandora_LArPandoraInterface
  larsim_Utils
  nusimdata::SimulationBase
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_fhicl()
install_source()
cet_enable_asserts()

add_subdirectory(template_generators)
