cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

add_subdirectory(EventWeight)
add_subdirectory(OscLib)
add_subdirectory(CAFAna)
add_subdirectory(FlatMaker)

