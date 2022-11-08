art_make_library(
    LIBRARIES
        larsim::PhotonPropagation
        larsim::PhotonPropagation_PhotonVisibilityService_service
        larsim::LegacyLArG4
        #sbncode_OpT0Finder_flashmatch_GeoAlgo
        sbncode_OpT0Finder_flashmatch_Base
        larcorealg::GeoAlgo
        art::Framework_Core
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        ROOT::Minuit
)

install_headers()
install_source()
