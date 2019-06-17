# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaExtrap
         LIB_LIBRARIES CAFAnaCore
        )
