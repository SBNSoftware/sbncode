link_directories(
 $ENV{ROOTSYS}/lib
)

art_make_library(
	LIBRARY_NAME sbn_LArReco_PCA
	SOURCE PCA.cc
	LIBRARIES
	sbnobj_Common_Reco
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase_dict	
	sbndcode_OpDetSim
)

simple_plugin(PCAnglePlaneMaker module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
	art::Framework_Principal
	art::Persistency_Common
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger
	
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	messagefacility::MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	sbn_LArReco_PCA
)

simple_plugin(PCAngleKinkFinder module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
	art::Framework_Principal
	art::Persistency_Common
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger
	
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	messagefacility::MF_MessageLogger
	Minuit

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
)

simple_plugin(PCAngleKinkTree module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
	art::Framework_Principal
	art::Persistency_Common
	art::Utilities canvas::canvas
	messagefacility::MF_MessageLogger
	
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	messagefacility::MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        
        art_root_io_tfile_support
        art_root_io
        art_root_io_dict
        caf_RecoUtils
        larsim_MCCheater_BackTrackerService_service
        nusimdata::SimulationBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	sbn_LArReco_PCA
)

install_headers()
install_fhicl()
install_source()

