include_directories( $ENV{IFDHC_FQ_DIR}/inc )

# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaCore
         LIB_LIBRARIES StandardRecordProxy
                       StandardRecord
                       OscLib
                       pthread
                       ${ROOT_BASIC_LIB_LIST}
                       ${IFDH}
        )

install_headers()
install_source()
