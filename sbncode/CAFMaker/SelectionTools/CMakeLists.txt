set( TOOL_LIBRARIES
               sbncode_CAFMaker
               sbnanaobj::StandardRecord
               larcorealg::Geometry
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               cetlib::cetlib cetlib_except::cetlib_except
               ROOT::Core ROOT::Tree
               art_root_io::RootDB
               hep_concurrency::hep_concurrency
               lardataobj::RecoBase
               nurandom::RandomUtils_NuRandomService_service
               lardata::DetectorClocksService
)

cet_build_plugin(CAFSelectionBlinded art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(CAFSelectionBlindEnergy art::tool LIBRARIES ${TOOL_LIBRARIES})
cet_build_plugin(CAFSelectionNtrkNproton art::tool LIBRARIES ${TOOL_LIBRARIES})

install_headers()
install_fhicl()
install_source()
