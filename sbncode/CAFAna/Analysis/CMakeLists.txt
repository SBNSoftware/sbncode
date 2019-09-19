# Allow "StandardRecord/" to find the right directory
include_directories(..)

art_make(LIBRARY_NAME CAFAnaAnalysis
         LIB_LIBRARIES CAFAnaCore CAFAnaExperiment
                       ${ROOT_MINUIT2}
        )

install_headers()
install_source()
