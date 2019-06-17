# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaCuts
         LIB_LIBRARIES CAFAnaCore
                       StandardRecordProxy
        )
