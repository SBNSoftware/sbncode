cet_enable_asserts()

GenerateRepoVersionSource(${CMAKE_PROJECT_NAME})

message(DEBUG "Generated source file '${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.h'")

art_make_library(
  LIBRARIES
    sbnobj::Common_Metadata
    larcorealg::headers
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
  )


set( sbn_meta_lib_list   art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         art::Framework_Services_System_TriggerNamesService_service
                         art::Persistency_Common canvas::canvas
                         art::Persistency_Provenance canvas::canvas
                         ROOT::Core
                         larcoreobj::SummaryData
    )

cet_build_plugin( FileCatalogMetadataSBN  art::service
               LIBRARIES ${sbn_meta_lib_list}
        )

cet_build_plugin( MetadataSBN art::service
               LIBRARIES
               ${sbn_meta_lib_list}
               art::Framework_IO
               art::Framework_IO_detail
        )

cet_build_plugin(SaveJobEnvironment art::ResultsProducer
  LIBRARIES
    sbncode::Metadata
    sbnobj::Common_Metadata
    art_root_io::TFileService_service
    art_root_io::tfile_support
    art::Framework_Services_Registry
  )

# cet_build_plugin(DumpJobEnvironment art::ResultsProducer
#   LIBRARIES
#     sbnobj::Common_Metadata
#   )

cet_build_plugin(DumpJobEnvironment art::Output
  LIBRARIES
    sbnobj::Common_Metadata
  )

cet_build_plugin(sbncodeRepositoryVersion art::tool
  LIBRARIES
    sbncode::Metadata
    sbnobj::Metadata
    sbnobj::Common_Metadata
  )


install_headers(LIST "${CMAKE_CURRENT_BINARY_DIR}/${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.h")
install_source(LIST "${CMAKE_CURRENT_BINARY_DIR}/${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.cxx")
install_headers()
install_fhicl()
install_source()
