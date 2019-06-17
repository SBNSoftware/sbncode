# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaVars
         LIB_LIBRARIES CAFAnaCore
                       StandardRecordProxy
                       OscLib
        )
