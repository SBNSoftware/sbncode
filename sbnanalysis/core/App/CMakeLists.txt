art_make(
  LIBRARY_NAME sbncode_SBNAnalysis
  LIB_LIBRARIES
    sbnanalysis_core_DataTypes_dict
    sbncode_sbnanalysis_core_Base_Processor
    ${Boost_SYSTEM_LIBRARY}
  MODULE_LIBRARIES
    sbnanalysis_core_DataTypes_dict
    sbncode_sbnanalysis_core_Base_Processor
    ${MF_MESSAGELOGGER}
 
)

install_headers()
install_fhicl()
install_source()

