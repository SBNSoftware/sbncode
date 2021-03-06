art_make_library(
        SOURCE StubBuilder.cxx PlaneTransform.cxx StubMergeAlgorithms.cxx
        LIBRARIES
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_PRINCIPAL}
	art_Persistency_Common
	art_Persistency_Provenance
	canvas art_Utilities
	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_BASIC_LIB_LIST}
	MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
	larreco_Calorimetry
	lardataobj_RecoBase
	lardataalg_DetectorInfo
)
simple_plugin(VertexChargeVacuum module
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_PRINCIPAL}
	${ART_FRAMEWORK_SERVICES}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	art_Persistency_Common
	art_Utilities canvas
	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
	larreco_Calorimetry
	lardataobj_RecoBase
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(VertexStubTracker module
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
	${ART_FRAMEWORK_PRINCIPAL}
	art_Persistency_Common
	art_Persistency_Provenance
	art_Utilities canvas
	${MF_MESSAGELOGGER}
	${MF_UTILITIES}
	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(PlaneStubMerge tool
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common
	art_Persistency_Provenance
        art_Utilities canvas
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        cetlib cetlib_except
        ${ROOT_EVE_LIB_LIST}
        ${ROOT_X3d}
        ${ROOT_BASIC_LIB_LIST}
        MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

simple_plugin(TwoPlaneStubMerge tool
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common
        art_Utilities canvas
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
        cetlib cetlib_except
        ${ROOT_EVE_LIB_LIST}
        ${ROOT_X3d}
        ${ROOT_BASIC_LIB_LIST}
        MF_MessageLogger

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	larreco_Calorimetry
	lardataalg_DetectorInfo
	sbncode_TPCReco_VertexStub
)

install_headers()
install_fhicl()
install_source()
