# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.

include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})
include_directories($ENV{NUG4_INC})

file(RELATIVE_PATH RELATIVE_CURRENT_SOURCE_DIR ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Add in ROOT dictionary
FUTURE_ROOT_GENERATE_DICTIONARY(TrackReducer_dict 
  ${RELATIVE_CURRENT_SOURCE_DIR}/Data.h  
  LINKDEF linkdef.h)
  
add_library(sbnanalysis_TrackReducer_Data SHARED TrackReducer_dict.cxx)
set_target_properties(sbnanalysis_TrackReducer_Data PROPERTIES PUBLIC_HEADER "Data.h;")
add_custom_command(TARGET sbnanalysis_TrackReducer_Data POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libTrackReducer_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

add_library(sbnanalysis_TrackReducer SHARED TrackReducer.cxx)

target_link_libraries(sbnanalysis_TrackReducer
  sbnanalysis_SBNOscReco_RecoUtils
  sbnanalysis_TrackReducer_Data
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_TrackReducer_Data DESTINATION lib PUBLIC_HEADER DESTINATION ${SBN_INC_DIR}/${RELATIVE_CURRENT_SOURCE_DIR})
install(TARGETS sbnanalysis_TrackReducer DESTINATION lib)


