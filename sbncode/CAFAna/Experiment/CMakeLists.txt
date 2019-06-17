# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaExperiment
         LIB_LIBRARIES CAFAnaCore CAFAnaVars CAFAnaPrediction
                       OscLib
        )
