cet_make_exec( extractCAFMetadata
               SOURCE extractCAFMetadata.cc
               LIBRARIES ${ROOT_BASIC_LIB_LIST}
               )

cet_script(diff_cafs)

install_headers()
install_source()
