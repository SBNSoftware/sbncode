art_make_library(
        SOURCE StubBuilder.cxx StubMergeAlgorithms.cxx
        LIBRARIES
	art::Framework_Core
	art::Framework_Principal
	art::Persistency_Common
	art::Persistency_Provenance
	canvas::canvas art::Utilities
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	${ROOT_BASIC_LIB_LIST}
	messagefacility::MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
	larreco_Calorimetry
	lardataobj_RecoBase
	lardataalg_DetectorInfo
)
simple_plugin(VertexChargeVacuum module
	art::Framework_Core
	art::Framework_Principal
	${ART_FRAMEWORK_SERVICES}
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Utilities canvas::canvas
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	messagefacility::MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
	larreco_Calorimetry
	lardataobj_RecoBase
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(VertexStubTracker module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
	art::Framework_Principal
	art::Persistency_Common
	art::Persistency_Provenance
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
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(PlaneStubMerge tool
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service                ${ART_FRAMEWORK_SERVICES}
        art::Framework_Principal
        art::Persistency_Common
	art::Persistency_Provenance
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
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(TwoPlaneStubMerge tool
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
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

add_subdirectory(config)

install_headers()
install_fhicl()
install_source()
