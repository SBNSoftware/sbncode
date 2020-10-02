# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME SBNAnaCuts
         LIB_LIBRARIES CAFAnaCore
                       SBNAnaVars
                       StandardRecordProxy
        )

install_headers()
install_source()
