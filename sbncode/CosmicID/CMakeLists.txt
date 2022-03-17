#link_directories(
# $ENV{ROOTSYS}/lib
#)

art_make( MODULE_LIBRARIES
	  ${ART_FRAMEWORK_SERVICES}
	  art_root_io::TFileService_service
	  fhiclcpp::fhiclcpp
	  art::Persistency_Provenance canvas
	  ${ROOT_BASIC_LIB_LIST}
	  ${ROOT_TMVA}
	  sbnobj_Common_Reco
	  sbncode_GeoWrappers
	  sbncode_LArRecoProducer
	  lardataobj_RecoBase_AssnsDicts_dict
	  lardataobj_RecoBase_dict	
	  lardataobj_RecoBase
	  larsim_Utils
)

install_fhicl()
install_source()

