art_make_library( 
    LIB_LIBRARIES
        #sbncode_OpT0Finder_flashmatch_GeoAlgo
        larcorealg::GeoAlgo
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardataalg::DetectorInfo
        larsim::PhotonPropagation
        lardata::DetectorPropertiesService
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        art::Utilities canvas::canvas
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        CLHEP::CLHEP
        Boost::system
        cetlib::cetlib cetlib_except::cetlib_except
)

install_headers()
install_source()
