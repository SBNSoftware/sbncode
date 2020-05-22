add_library(sbnanalysis_SBNOscReco_RecoUtils SHARED RecoUtils.cc GeoUtil.cc)
include_directories($ENV{NUG4_INC})

target_link_libraries(sbnanalysis_SBNOscReco_RecoUtils
  sbnanalysis_ProviderManager
  lardataobj_RawData
  lardata_RecoObjects
  lardataobj_RecoBase
)

install(TARGETS sbnanalysis_SBNOscReco_RecoUtils DESTINATION lib)

