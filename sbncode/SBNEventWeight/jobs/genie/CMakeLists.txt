#install_headers()
install_fhicl()
#install_source()
FILE(GLOB fcl_files *.fcl)

install_source( EXTRAS ${fcl_files} )


