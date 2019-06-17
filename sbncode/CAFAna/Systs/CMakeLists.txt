# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaSysts
         LIB_LIBRARIES CAFAnaCore
                       StandardRecordProxy
        )

install_headers()
install_source()
