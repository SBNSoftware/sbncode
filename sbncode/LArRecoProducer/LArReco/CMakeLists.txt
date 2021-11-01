
art_make_library( LIBRARY_NAME sbn_LArReco
  SOURCE  TrackMomentumCalculator.cxx TrajectoryMCSFitter.cxx
                  LIBRARIES
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES}
        ${ART_FRAMEWORK_PRINCIPAL}
        art_Persistency_Common
        art_Utilities canvas
                    ${ROOT_BASIC_LIB_LIST}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${FHICLCPP}
                    cetlib cetlib_except
        lardataobj_RecoBase
        ${ROOT_MINUIT}
        ${ROOT_MINUIT2}
                )
