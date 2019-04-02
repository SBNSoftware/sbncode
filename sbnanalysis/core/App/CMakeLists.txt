art_make(
  LIBRARY_NAME sbncode_SBNAnalysis
  LIB_LIBRARIES
    sbnanalysis_core_DataTypes_dict
    sbncode_sbnanalysis_core_Base_Processor
    ${Boost_SYSTEM_LIBRARY}
    # Processors
    sbncode_sbnanalysis_ana_ExampleAnalysis
  MODULE_LIBRARIES
    sbnanalysis_core_DataTypes_dict
    sbncode_sbnanalysis_core_Base_Processor
    ${MF_MESSAGELOGGER}
    # Processors
    sbncode_sbnanalysis_ana_ExampleAnalysis
)

install_headers()
install_fhicl()
install_source()

