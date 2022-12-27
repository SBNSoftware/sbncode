find_package(ifbeam)
find_package(ifdh_art)

cet_build_plugin(NuMIRetriever art::module
    LIBRARIES
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        sbnobj::Common_POTAccounting
        ifbeam::ifbeam
        ifdh_art::IFBeam_service
        lardata::Utilities
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        artdaq_core::artdaq-core_Utilities
)

install_headers()
install_fhicl()
install_source()
