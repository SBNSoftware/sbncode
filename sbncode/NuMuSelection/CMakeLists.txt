simple_plugin( MuPVertexStudy module
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
	art_root_io_tfile_support
	art_root_io
	art_root_io_dict
        caf_RecoUtils
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
        larsim_MCCheater_BackTrackerService_service
        nusimdata::SimulationBase
        larsim_MCCheater_ParticleInventoryService_service
	larcorealg_GeoAlgo
)
simple_plugin( MuonS2NStudy module
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
	art_root_io_tfile_support
	art_root_io
	art_root_io_dict
        caf_RecoUtils
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
        larsim_MCCheater_BackTrackerService_service
        nusimdata::SimulationBase
        larsim_MCCheater_ParticleInventoryService_service
	larcorealg_GeoAlgo
)
simple_plugin( NuMuEfficiencyStudy module
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
	art_root_io_tfile_support
	art_root_io
	art_root_io_dict
        caf_RecoUtils
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
        larsim_MCCheater_BackTrackerService_service
        nusimdata::SimulationBase
        larsim_MCCheater_ParticleInventoryService_service
	larcorealg_GeoAlgo
)

install_headers()
install_fhicl()
install_source()
