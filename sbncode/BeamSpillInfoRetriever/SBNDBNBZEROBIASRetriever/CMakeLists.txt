find_package(ifbeam)
find_package(ifdh_art)

cet_build_plugin(SBNDBNBZEROBIASRetriever art::module
    LIBRARIES
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::X3d
        Boost::system
        messagefacility::MF_MessageLogger
        ifbeam::ifbeam
        ifdh_art::IFBeam_service
	SQLite::SQLite3
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
        artdaq_core::artdaq-core_Utilities
        sbn_MWRData
        sbnobj::Common_POTAccounting
        larcorealg::CoreUtils
)

install_headers()
install_fhicl()
install_source()

