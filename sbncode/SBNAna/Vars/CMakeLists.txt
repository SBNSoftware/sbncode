# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME SBNAnaVars
         LIB_LIBRARIES CAFAnaCore
                       StandardRecordProxy
                       OscLib
        )

install_headers()
install_source()
