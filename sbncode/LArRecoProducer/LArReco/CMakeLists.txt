
art_make_library( LIBRARY_NAME sbn_LArReco
  SOURCE  TrackMomentumCalculator.cxx TrajectoryMCSFitter.cxx
                  LIBRARIES
        art::Framework_Core
        art::Framework_Services_Registry
        ${ART_FRAMEWORK_SERVICES}
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
                    ${ROOT_BASIC_LIB_LIST}
        messagefacility::MF_MessageLogger
        
        fhiclcpp::fhiclcpp
                    cetlib::cetlib cetlib_except
        lardataobj_RecoBase
        ${ROOT_MINUIT}
        ${ROOT_MINUIT2}
                )
