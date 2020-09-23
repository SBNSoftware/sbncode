add_library(sbnanalysis_SBNOscReco_LArReco SHARED TrackMomentumCalculator.cxx TrajectoryMCSFitter.cxx)

target_link_libraries(sbnanalysis_SBNOscReco_LArReco
  lardataobj_RecoBase
  Minuit
  Minuit2
)

install(TARGETS sbnanalysis_SBNOscReco_LArReco DESTINATION lib)

