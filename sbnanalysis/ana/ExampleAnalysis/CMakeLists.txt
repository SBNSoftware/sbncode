project(ExampleAnalysis)

include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_library(sbnanalysis_ExampleAnalysis_ExampleTools SHARED ExampleTools.cxx)
add_library(sbnanalysis_ExampleAnalysis_ExampleSelection SHARED ExampleSelection.cxx)

add_executable(ExampleAnalysis ExampleAnalysis.cxx)

target_link_libraries(
  ExampleAnalysis
  sbnanalysis_ExampleAnalysis_ExampleTools
  sbnanalysis_ExampleAnalysis_ExampleSelection
  sbnanalysis_Processor
  canvas
  cetlib
  gallery
  nusimdata_SimulationBase
  larcoreobj_SummaryData
  lardataobj_RecoBase
)

install(TARGETS sbnanalysis_ExampleAnalysis_ExampleTools DESTINATION lib)
install(TARGETS sbnanalysis_ExampleAnalysis_ExampleSelection DESTINATION lib)

install(TARGETS ExampleAnalysis DESTINATION ${PROJECT_NAME}/bin)

