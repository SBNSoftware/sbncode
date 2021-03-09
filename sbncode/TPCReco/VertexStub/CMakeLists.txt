art_make_library(LIBRARY_NAME sbncode_StubBuilder
        SOURCE StubBuilder.cxx
        LIBRARIES
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
	lardataobj_RecoBase
	lardataalg_DetectorInfo
)

art_make_library(LIBRARY_NAME sbncode_StubMergeAlgorithms
        SOURCE StubMergeAlgorithms.cxx
        LIBRARIES
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
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	lardataalg_DetectorInfo
)

simple_plugin(VertexChargeVacuum module
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
	lardataobj_RecoBase
	lardataalg_DetectorInfo
)

simple_plugin(VertexStubTracker module
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
	sbncode_StubBuilder
)

simple_plugin(PlaneStubMerge tool
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
	sbncode_StubMergeAlgorithms
)

install_headers()
install_fhicl()
install_source()
