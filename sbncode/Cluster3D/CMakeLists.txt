
# where should the scripts/..xml file be installed?  Perhaps in bin?

art_make_library(
        LIBRARIES
           larevt::Filters
           lardataalg::DetectorInfo
           lardataobj::RecoBase
           larcorealg::Geometry
           lardata::ArtDataHelper
           lardata::RecoObjects
           art::Framework_Core
           art::Framework_Principal
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art::Utilities
           ROOT::Core
           ROOT::Gpad
           ROOT::Hist
           ROOT::Matrix
           ROOT::Physics
           ROOT::Tree
           art_root_io::tfile_support
           art_root_io::TFileService_service
           art::Persistency_Provenance
           canvas::canvas
           messagefacility::MF_MessageLogger
           fhiclcpp::fhiclcpp
           cetlib::cetlib
           cetlib_except::cetlib_except)
set(    TOOL_LIBRARIES
           larevt::Filters
           larreco::RecoAlg_Cluster3DAlgs
           lardataalg::DetectorInfo
           lardataobj::RecoBase
           lardata::ArtDataHelper
           art::Framework_Services_Registry
           art_root_io::tfile_support
           art_root_io::TFileService_service
           ROOT::Tree)
set(    MODULE_LIBRARIES
          larreco::Calorimetry
          larreco::RecoAlg_Cluster3DAlgs
          larreco::RecoAlg_TCAlg
          larreco::RecoAlg
          larreco::ClusterFinder
          larsim::MCCheater_ParticleInventoryService_service
          lardataobj::AnalysisBase
          lardataobj::RecoBase
          nurandom::RandomUtils_NuRandomService_service
          nusimdata::SimulationBase
          art::Framework_Services_Registry
          art_root_io::tfile_support
          ROOT::Core
          ROOT::Physics
          art_root_io::TFileService_service
          messagefacility::MF_MessageLogger
       )
cet_build_plugin(Cluster3DSBN art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SnippetHit3DBuilderSBN art::tool LIBRARIES ${TOOL_LIBRARIES} SOURCE SnippetHit3DBuilderSBN_tool.cc)
set_property(SOURCE SnippetHit3DBuilderSBN_tool.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)

install_headers()
install_fhicl()
install_source()