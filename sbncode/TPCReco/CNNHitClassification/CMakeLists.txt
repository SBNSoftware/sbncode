cet_build_plugin(CNNID art::EDProducer
  LIBRARIES PRIVATE
    larrecodnn::ImagePatternAlgs_Tensorflow_PointIdAlg
    lardata::ArtDataHelper
    lardata::DetectorClocksService
    lardata::DetectorPropertiesService
    lardata::AssociationUtil
    lardataobj::RecoBase
    sbnobj::Common_Reco
    larcoreobj::SimpleTypesAndConstants
    art::Framework_Services_System_TriggerNamesService_service
    art::Framework_Services_Registry
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::types
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
)

cet_build_plugin(InspectCNNScores art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::ArtDataHelper
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(EmTrackMichelIdLocal art::EDProducer
  LIBRARIES PRIVATE
  larrecodnn::ImagePatternAlgs_Tensorflow_PointIdAlg
  lardata::ArtDataHelper
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  lardata::AssociationUtil
  lardataobj::RecoBase
  larcoreobj::SimpleTypesAndConstants
  art::Framework_Services_System_TriggerNamesService_service
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()
