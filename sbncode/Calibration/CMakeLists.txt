link_directories(
 $ENV{ROOTSYS}/lib
)

simple_plugin( TrackCaloSkimmer module
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
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
	sbnobj_Common_Reco
	larcorealg_GeoAlgo
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
	Minuit
	sbncode_CAFMaker
)

simple_plugin(TrackCaloSkimmerSelectStoppingTrack tool
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
        larcorealg_Geometry
	larcorealg_GeoAlgo
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
)

simple_plugin(TrackCaloSkimmerSelectAnode2CathodeTrack tool
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
)

simple_plugin ( DataSelect module
                ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               ${ART_FRAMEWORK_SERVICES_BASIC}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               art::Utilities canvas::canvas
               messagefacility::MF_MessageLogger
               fhiclcpp::fhiclcpp
               cetlib::cetlib cetlib_except
        ${ROOT_EVE_LIB_LIST}
        ${ROOT_X3d}
        ${ROOT_BASIC_LIB_LIST}
)


add_subdirectory(fcl)

install_headers()
install_fhicl()
install_source()

