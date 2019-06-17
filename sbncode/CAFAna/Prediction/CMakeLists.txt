# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaPrediction
         LIB_LIBRARIES CAFAnaCore CAFAnaSysts CAFAnaExtrap
        )
