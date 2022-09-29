cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAM NAMES ifbeam PATHS ENV IFBEAM_LIB NO_DEFAULT_PATH )
cet_find_library( WDA NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )

find_ups_product(sbnanaobj)

add_subdirectory(SinglePhotonAnalysis)
add_subdirectory(SBNEventWeight)
add_subdirectory(NuMuSelection)
add_subdirectory(OpDet)
add_subdirectory(OpT0Finder)

add_subdirectory(CAFMaker)
add_subdirectory(Metadata)

add_subdirectory(Calibration)
add_subdirectory(FlashMatch)
add_subdirectory(LArRecoProducer)
add_subdirectory(TPCReco)


add_subdirectory(BeamSpillInfoRetriever)
add_subdirectory(FluxReader)
add_subdirectory(EventGenerator)
add_subdirectory(PID)
add_subdirectory(GeometryTools)
add_subdirectory(CosmicID)

# Supera
if ( NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Supera/.git" )
  execute_process(COMMAND git submodule update --init --recursive
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# We need to run Supera/setup.sh first, otherwise it won't have a CMakeLists.txt.
execute_process(COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/Supera/setup.sh icarus
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Supera )
add_subdirectory(Supera)
