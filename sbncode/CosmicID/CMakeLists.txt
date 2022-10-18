art_make( MODULE_LIBRARIES
          art_root_io::TFileService_service
          fhiclcpp::fhiclcpp
          art::Persistency_Provenance canvas
          sbnobj::Common_Reco
          sbncode::GeoWrappers
          sbncode::LArRecoProducer
          lardataobj::RecoBase_AssnsDicts_dict
          lardataobj::RecoBase_dict
          lardataobj::RecoBase
          larsim::Utils
          ROOT::TMVA
)

install_fhicl()
install_source()

