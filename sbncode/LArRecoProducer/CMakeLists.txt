
add_subdirectory(LArReco)
art_make_library(
  LIBRARIES
                fhiclcpp::fhiclcpp
                art::Persistency_Provenance canvas::canvas
                sbncode_GeoWrappers
)
set(  MODULE_LIBRARIES
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::X3d
        ROOT::Core
        messagefacility::MF_MessageLogger
        ROOT::EG
        sbn_LArReco
        sbncode_LArRecoProducer
        sbnobj::Common_CRT
        sbnobj::Common_Reco
        lardataobj::RecoBase_AssnsDicts_dict
        lardataobj::RecoBase_dict
        lardataobj::RecoBase
)

cet_build_plugin(CRTHitConverter art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(MCSFitAllPID art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTFlashTriggerFilter art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(PMTFlashTriggerMaker art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(RangePAllPID art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ShowerCosmicDistance art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(ShowerSelectionVars art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TrackScatterClosestApproach art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TrackStoppingChi2Fitter art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TransferPFParticleFlashMatch art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TransferPFParticleT0 art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(TransferTrackT0 art::module LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

