add_subdirectory(ICARUSBNBRetriever)
add_subdirectory(ICARUSBNBEXTRetriever)
add_subdirectory(SBNDBNBRetriever)
add_subdirectory(SBNDBNBEXTRetriever)
add_subdirectory(SBNDBNBZEROBIASRetriever)
add_subdirectory(NuMIRetriever)
add_subdirectory(NuMIEXTRetriever)
add_subdirectory(job)

find_package(ifbeam)
find_package(ifdh_art)

art_make_library(LIBRARIES Boost::system
        LIBRARY_NAME sbn_MWRData
        SOURCE MWRData.cpp
)

art_make_library(
    LIBRARIES
        art::Persistency_Common
        art::Utilities 
        ifbeam::ifbeam
        ifdh_art::IFBeam_service
        messagefacility::MF_MessageLogger
        sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
        sbnobj::Common_POTAccounting
        sbn_MWRData
        larcorealg::CoreUtils

        LIBRARY_NAME sbn_POTTools
        SOURCE POTTools.cpp
)

install_headers()
install_fhicl()
install_source()
