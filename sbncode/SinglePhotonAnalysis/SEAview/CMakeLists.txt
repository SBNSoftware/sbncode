link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)
find_package(larpandora)

art_make( 
	LIBRARY_NAME sbncode_SinglePhotonAnalysis_SEAview
    LIB_LIBRARIES
      larcorealg_Geometry
      larcore_Geometry_Geometry_service
      lardata_Utilities
      larevt_Filters
      lardataobj_RawData
      lardataobj_RecoBase
      lardataobj_MCBase
      larcoreobj_SummaryData
      larsim_MCCheater_BackTrackerService_service
      larsim_MCCheater_ParticleInventoryService_service
      larpandora_LArPandoraInterface
      LArPandoraContent 
      nusimdata_SimulationBase
      ${Boost_SERIALIZATION_LIBRARY}
      nurandom::RandomUtils_NuRandomService_service
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
      ${MF_MESSAGELOGGER}

      ${FHICLCPP}
      cetlib cetlib_except
      ${CLHEP}
      ${ROOT_BASIC_LIB_LIST}
      ${ROOT_GEOM} 
    MODULE_LIBRARIES
      larcorealg_Geometry
      larcore_Geometry_Geometry_service
      lardata_Utilities
      larevt_Filters
      lardataobj_RawData
      lardataobj_RecoBase
      lardataobj_MCBase
      larcoreobj_SummaryData	
      larsim_MCCheater_BackTrackerService_service
      larpandora_LArPandoraInterface
      SignalShapingServiceMicroBooNE_service
      nusimdata_SimulationBase
      ${Boost_SERIALIZATION_LIBRARY}
      nurandom::RandomUtils_NuRandomService_service
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
      ${MF_MESSAGELOGGER}
      ${FHICLCPP}
      cetlib cetlib_except
      ${CLHEP}
      ${ROOT_BASIC_LIB_LIST}
      ${ROOT_GEOM} )
