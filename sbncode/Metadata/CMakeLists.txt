cet_enable_asserts()

set( sbn_meta_lib_list   ${ART_FRAMEWORK_CORE}
                         ${ART_FRAMEWORK_PRINCIPAL}
                         ${ART_FRAMEWORK_SERVICES_REGISTRY}
                         ${ART_FRAMEWORK_SERVICES_BASIC}
                         ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                         art_Framework_Services_System_TriggerNamesService_service
                         ${MF_MESSAGELOGGER}
                         ${ROOT_BASIC_LIB_LIST}
    )

simple_plugin( FileCatalogMetadataSBN  "service"
               ${sbn_meta_lib_list}
        )

simple_plugin( TFileMetadataSBN "service"
               ${sbn_meta_lib_list}
               art_Framework_IO
               art_Framework_IO_detail
        )

install_headers()
install_fhicl()
install_source()
