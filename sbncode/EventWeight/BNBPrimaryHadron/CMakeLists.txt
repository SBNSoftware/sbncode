include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{NUTOOLS_INC})
include_directories($ENV{GENIE_INC}/GENIE)
file(GLOB WGHCALCFILES *WeightCalc*)

art_make_library(
  LIBRARY_NAME sbncode_BNBPrimaryHadronWeightCalculators
  SOURCE ${WGHCALCFILES}
  LIBRARIES
    sbncode_WeightCalculators
    lardata_Utilities
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${CLHEP}
) 

install_headers()
install_fhicl()
install_source()

