cet_enable_asserts()

set( sbn_meta_lib_list   art::Framework_Core
                         art::Framework_Principal
                         art::Framework_Services_Registry
                         art::Framework_Services_System_TriggerNamesService_service
                         art::Persistency_Common canvas::canvas
                         art::Persistency_Provenance canvas::canvas
                         ROOT::Core
                         larcoreobj::SummaryData
    )

simple_plugin( FileCatalogMetadataSBN  "service"
               ${sbn_meta_lib_list}
        )

simple_plugin( MetadataSBN "service"
               ${sbn_meta_lib_list}
               art::Framework_IO
               art::Framework_IO_detail
        )

install_headers()
install_fhicl()
install_source()
