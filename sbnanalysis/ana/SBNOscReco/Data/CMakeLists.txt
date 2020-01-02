file(RELATIVE_PATH RELATIVE_CURRENT_SOURCE_DIR ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Add in ROOT dictionary
FUTURE_ROOT_GENERATE_DICTIONARY(SBNOscReco_dict 
  ${RELATIVE_CURRENT_SOURCE_DIR}/CRTMatch.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/FlashMatch.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/Mode.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/RecoEvent.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/RecoParticle.h  
  ${RELATIVE_CURRENT_SOURCE_DIR}/RecoTrack.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/TrueParticle.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/TruthMatch.h 
  ${RELATIVE_CURRENT_SOURCE_DIR}/DetInfo.h 
  ${RELATIVE_CURRENT_SOURCE_DIR}/MCType.h
  ${RELATIVE_CURRENT_SOURCE_DIR}/FlatInteraction.h
  LINKDEF linkdef.h)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#REFLEX_GENERATE_DICTIONARY(SBNOscReco_dict classes.h
#  SELECTION classes_def.xml)
  
add_library(sbnanalysis_SBNOscReco_classes SHARED SBNOscReco_dict.cxx)
set_target_properties(sbnanalysis_SBNOscReco_classes PROPERTIES PUBLIC_HEADER "CRTMatch.h;FlashMatch.h;Mode.h;RecoEvent.h;RecoParticle.h;RecoTrack.h;TruthMatch.h;DetInfo.h;MCType.h;TrueParticle.h;")
add_custom_command(TARGET sbnanalysis_SBNOscReco_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libSBNOscReco_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_classes DESTINATION lib PUBLIC_HEADER DESTINATION ${SBN_INC_DIR}/${RELATIVE_CURRENT_SOURCE_DIR})
