add_subdirectory(TrackHit)
add_subdirectory(TrackSplit)
add_subdirectory(VertexStub)
add_subdirectory(CalorimetryAnalysis)
add_subdirectory(BlipReco)


cet_build_plugin(NuVertexChargeTree art::module
    LIBRARIES
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::X3d
        messagefacility::MF_MessageLogger
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        sbnobj::Common_Reco
        lardataobj::AnalysisBase
        lardataalg::DetectorInfo
        art_root_io::tfile_support
        art_root_io::art_root_io
        art_root_io::dict
        caf_RecoUtils
        larsim::MCCheater_BackTrackerService_service
        nusimdata::SimulationBase
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        larpandora::LArPandoraInterface
        larsim::MCCheater_ParticleInventoryService_service
        sbn_LArReco
        sbncode_TPCReco_VertexStub
)


cet_build_plugin( EventSelect art::module
     LIBRARIES
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               cetlib::cetlib cetlib_except::cetlib_except
               lardataalg::DetectorInfo
               ROOT::X3d
)


install_headers()
install_fhicl()
install_source()

