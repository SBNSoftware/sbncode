find_package(ifbeam)
find_package(ifbeam_art)
include_directories( $ENV{IFBEAM_FQ_DIR}/include )
include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFDH_ART_INC} )
cet_find_library( IFBEAMSERVICE        NAMES IFBeam_service    PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAM               NAMES ifbeam            PATHS ENV IFBEAM_LIB   NO_DEFAULT_PATH )



simple_plugin(NuMIRetriever module
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib cetlib_except
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
