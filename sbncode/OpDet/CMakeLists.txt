cet_enable_asserts()

cet_make_library(
        LIBRARY_NAME OpDet_PDMapAlg INTERFACE
        SOURCE
        PDMapAlg.h
)

cet_build_plugin(PDMapAlgSimple art::tool LIBRARIES PRIVATE
            art::Utilities
            sbncode::OpDet_PDMapAlg
)

install_headers()
install_source()
install_fhicl()
