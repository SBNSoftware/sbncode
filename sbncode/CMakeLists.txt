cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

find_ups_product(sbnanaobj)

add_subdirectory(SBNEventWeight)
add_subdirectory(NuMuSelection)
add_subdirectory(OpDet)
add_subdirectory(OpT0Finder)

add_subdirectory(CAFMaker)

add_subdirectory(FlashMatch)
add_subdirectory(LArRecoProducer)

