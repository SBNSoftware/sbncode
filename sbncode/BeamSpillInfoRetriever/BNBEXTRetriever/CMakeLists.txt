
cet_build_plugin(BNBEXTRetriever art::module
    LIBRARIES
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        sbnobj::Common_POTAccounting
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        artdaq_core::artdaq-core_Utilities
        lardata::Utilities
        larcore::Geometry_AuxDetGeometry_service
)

install_headers()
install_fhicl()
install_source()

