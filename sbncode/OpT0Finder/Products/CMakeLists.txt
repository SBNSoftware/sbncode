# include_directories( $ENV{NUTOOLS_INC} )
art_make_library(
    LIBRARIES
        lardataobj::RecoBase
        ROOT::Core
        ROOT::EG
)

art_dictionary(
    DICTIONARY_LIBRARIES
        larcorealg::Geometry
        larcoreobj::SimpleTypesAndConstants
        lardataobj::RecoBase
        ROOT::Core
        ROOT::EG
)
install_headers()
install_fhicl()
install_source()
