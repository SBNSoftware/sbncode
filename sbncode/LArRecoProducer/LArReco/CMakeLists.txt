art_make_library( NAME sbn_LArReco
                  LIBRARIES
                       art::Framework_Core
                       art::Framework_Services_Registry
                       art::Framework_Principal
                       art::Persistency_Common
                       art::Utilities canvas::canvas
                       cetlib::cetlib cetlib_except::cetlib_except
                       lardataobj::RecoBase larcorealg::CoreUtils lardata::RecoObjects
                       ROOT::Minuit
                       ROOT::Minuit2
                       ROOT::Graf3d
                )
