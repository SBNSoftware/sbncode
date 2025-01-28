art_make(  MODULE_LIBRARIES  larpandora::LArPandoraInterface
                             lardata::RecoObjects  
                             lardataobj::RecoBase
                             lardataobj::AnalysisBase
                             larreco::RecoAlg

                             art::Framework_Core
                             art::Framework_Principal
                             art::Framework_Services_Registry
                             art_root_io::tfile_support
                             art_root_io::TFileService_service
                             artdaq_core::artdaq-core_Utilities
                             sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
                             sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_SBND
                             art::Persistency_Common
                             art::Persistency_Provenance
                             art::Utilities

                             cetlib::cetlib
                             cetlib_except::cetlib_except
                             ${LARRECO_LIB}
                             ${LARDATA_LIB}
                             ${ART_FRAMEWORK_SERVICES_BASIC}
                             ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                             ${MF_UTILITIES}
                             ${FHICLCPP}
                             ${CLHEP}
                             sbnobj::Common_Reco


                             messagefacility::MF_MessageLogger
                             fhiclcpp::fhiclcpp
                             ROOT::Geom
                             ROOT::XMLIO
                             ROOT::Gdml
                             ROOT::Core

                             ${ROOT_BASIC_LIB_LIST}
                             ${ROOT_SPECTRUM}
                             ${ROOT_ROOFIT}
                             ${ROOT_ROOFITCORE}
        )
install_fhicl()
install_source()
