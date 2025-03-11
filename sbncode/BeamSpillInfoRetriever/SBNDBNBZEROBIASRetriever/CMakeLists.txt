find_package(ifbeam)
find_package(ifdh_art)

cet_build_plugin(SBNDBNBZEROBIASRetriever art::module
    LIBRARIES
        sbn_SBNDPOTTools
)

install_headers()
install_fhicl()
install_source()

