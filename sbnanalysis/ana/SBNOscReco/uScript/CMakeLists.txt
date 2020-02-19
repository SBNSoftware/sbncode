# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(uscript_lib chunk.cc value.cc vm.cc compile.cc scanner.cc tclass.cc common.cc)

target_link_libraries(uscript_lib
  ${ROOT_LIBRARIES}
)

add_executable(uscript main.cc)

target_link_libraries(uscript
  uscript_lib
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_Event
)

install(TARGETS uscript_lib DESTINATION lib)
install(TARGETS  uscript DESTINATION bin)
