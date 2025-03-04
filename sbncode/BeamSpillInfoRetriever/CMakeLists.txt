add_subdirectory(BNBRetriever)
add_subdirectory(SBNDBNBRetriever)
add_subdirectory(SBNDBNBZEROBIASRetriever)
add_subdirectory(SBNDBNBEXTRetriever)
add_subdirectory(NuMIRetriever)
add_subdirectory(BNBEXTRetriever)
add_subdirectory(NuMIEXTRetriever)
add_subdirectory(job)

art_make_library(LIBRARIES Boost::system
        LIBRARY_NAME sbn_MWRData
        SOURCE MWRData.cpp
)

art_make_library(
    LIBRARIES
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        ROOT::X3d
        Boost::system
        messagefacility::MF_MessageLogger
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
        artdaq_core::artdaq-core_Utilities

        LIBRARY_NAME sbn_SBNDPOTTools
        SOURCE SBNDPOTTools.cpp
)

install_headers()
install_fhicl()
install_source()
