simple_plugin(TrackAreaHit module
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
        lardataobj_RecoBase_AssnsDicts_dict
        lardataobj_RecoBase     
        lardataobj_RecoBase_dict        

        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
)

simple_plugin(TrackHitFilter module
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
        lardataobj_RecoBase_AssnsDicts_dict
        lardataobj_RecoBase     
        lardataobj_RecoBase_dict        

        larcorealg_Geometry
        larcore_Geometry_Geometry_service
)

simple_plugin(TrackHitDumper module
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
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase	
	lardataobj_RecoBase_dict	

	larcorealg_Geometry
	larcore_Geometry_Geometry_service
)

install_headers()
install_fhicl()
install_source()

