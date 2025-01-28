add_subdirectory(BNBRetriever)
add_subdirectory(SBNDBNBRetriever)
add_subdirectory(SBNDBNBZEROBIASRetriever)
add_subdirectory(SBNDBNBEXTRetriever)
add_subdirectory(NuMIRetriever)
add_subdirectory(BNBEXTRetriever)
add_subdirectory(NuMIEXTRetriever)
add_subdirectory(SBNDBNBLightAnalyzer)
add_subdirectory(job)


art_make_library(LIBRARIES Boost::system
        LIBRARY_NAME sbn_MWRData
        SOURCE MWRData.cpp

)



install_headers()
install_fhicl()
install_source()
