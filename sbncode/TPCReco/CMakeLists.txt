add_subdirectory(TrackHit)
add_subdirectory(TrackSplit)
add_subdirectory(VertexStub)
add_subdirectory(CalorimetryAnalysis)

link_directories(
 $ENV{ROOTSYS}/lib
)

simple_plugin(NuVertexChargeTree module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service
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
	larpandora_LArPandoraInterface
        larsim_MCCheater_ParticleInventoryService_service
	sbn_LArReco	    
	sbncode_TPCReco_VertexStub
)
simple_plugin ( EventSelect module
                ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               messagefacility::MF_MessageLogger
               fhiclcpp::fhiclcpp
               cetlib::cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()

