# Add in ROOT dictionary
ROOT_GENERATE_DICTIONARY(SBNOscReco_dict CRTMatch.h  FlashMatch.h  Mode.h  RecoEvent.h  RecoParticle.h  RecoTrack.h  TruthMatch.h LINKDEF linkdef.h)
add_library(sbnanalysis_SBNOscReco_classes SHARED SBNOscReco_dict.cxx)
add_custom_command(TARGET sbnanalysis_SBNOscReco_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libSBNOscReco_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_classes DESTINATION lib)


