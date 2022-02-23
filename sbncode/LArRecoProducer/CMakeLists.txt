#link_directories(
# $ENV{ROOTSYS}/lib
#)

add_subdirectory(LArReco)
art_make( LIB_LIBRARIES
		fhiclcpp::fhiclcpp
		art::Persistency_Provenance canvas
		${ROOT_BASIC_LIB_LIST}
		sbncode_GeoWrappers

	  MODULE_LIBRARIES
		art::Framework_Core
		art::Framework_Services_Registry
		art_root_io::TFileService_service
		${ART_FRAMEWORK_SERVICES}
		art::Utilities canvas::canvas
		art::Framework_Principal
		art::Persistency_Common
		messagefacility::MF_MessageLogger
		cetlib::cetlib cetlib_except
		${ROOT_EVE_LIB_LIST}
		${ROOT_X3d}
		${ROOT_BASIC_LIB_LIST}
		${ROOT_TMVA}
		messagefacility::MF_MessageLogger
		sbn_LArReco	    
		sbncode_LArRecoProducer
		sbnobj_Common_CRT
 		sbnobj_Common_Reco
		lardataobj_RecoBase_AssnsDicts_dict
		lardataobj_RecoBase_dict	
		lardataobj_RecoBase
		larsim_Utils
)

install_headers()
install_fhicl()
install_source()

