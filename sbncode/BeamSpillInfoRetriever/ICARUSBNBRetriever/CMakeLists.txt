find_package(ifbeam)
find_package(ifdh_art)

cet_build_plugin(ICARUSBNBRetriever art::module
    LIBRARIES
        sbn_SBNDPOTTools
        SQLite::SQLite3
)

install_headers()
install_fhicl()
install_source()
