cet_enable_asserts()

cet_build_plugin(PDMapAlgSimple art::tool LIBRARIES PRIVATE
            art::Utilities
)

install_headers()
install_source()
install_fhicl()
