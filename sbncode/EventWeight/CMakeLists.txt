include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{NUTOOLS_INC})
include_directories($ENV{GENIE_INC}/GENIE)

file(GLOB WGHCALCFILES *WeightCalc* Factory*)

art_make_library(
  LIBRARY_NAME sbncode_WeightCalculators 
  SOURCE ${WGHCALCFILES}
  LIBRARIES
    lardata_Utilities
    nutools_NuReweight_art
    nutools_NuReweight
    larsim_Simulation
    nutools_RandomUtils_NuRandomService_service
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    art_Persistency_Common canvas canvas
    art_Persistency_Provenance canvas canvas
    art_Utilities canvas canvas
    cetlib cetlib_except
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
)

art_make(
  LIBRARY_NAME sbncode_EventWeight
  EXCLUDE ${WGHCALCFILES}
  LIB_LIBRARIES
    sbncode_BNBPrimaryHadronWeightCalculators
    sbncode_WeightCalculators
    ${Boost_SYSTEM_LIBRARY}
  MODULE_LIBRARIES
    sbncode_WeightCalculators
    sbncode_BNBPrimaryHadronWeightCalculators
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
add_subdirectory(jobs)

