add_subdirectory(DataTypes)
add_subdirectory(fcl)

art_make(
          MODULE_LIBRARIES 
			   uboonecode_uboone_BasicTool_GeoAlgo


			   larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           #larsim_MCCheater_BackTracker_service
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_ArtDataHelper
                           lardata_RecoObjects larreco_Deprecated
                           lardataobj_RecoBase
                           lardata_RecoObjects larreco_Deprecated
                           lardataobj_AnalysisBase
                           lardata_Utilities
                           larpandora_LArPandoraInterface
                           larreco_Calorimetry
                           larreco_RecoAlg
                           larreco_RecoAlg_ClusterRecoUtil
                           larreco_ClusterFinder
                           nutools_RandomUtils_NuRandomService_service
                           nusimdata_SimulationBase
                           ${LARRECO_LIB}
                           ${LARDATA_LIB}
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
                           art_Persistency_Common canvas
                           art_Persistency_Provenance canvas
                           art_Utilities canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${CLHEP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_SPECTRUM}
                           ${ROOT_ROOFIT}
                           ${ROOT_ROOFITCORE}
        )
install_headers()
install_fhicl()
install_source()
