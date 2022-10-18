find_package(ifbeam)
find_package(ifbeam_art)


cet_make_library(LIBRARIES Boost::system
        LIBRARY_NAME sbn_BNBSpillInfoRetriever_MWRData
        SOURCE MWRData.cpp
)


simple_plugin(BNBRetriever module
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::X3d
        Boost::system
        messagefacility::MF_MessageLogger
        ifbeam::ifbeam
        ifdh_art::IFBeam_service
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        artdaq_core::artdaq-core_Utilities
        sbn_BNBSpillInfoRetriever_MWRData
        sbnobj::Common_POTAccounting
        larcorealg::CoreUtils
)

install_headers()
install_fhicl()
install_source()

