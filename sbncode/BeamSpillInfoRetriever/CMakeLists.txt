add_subdirectory(BNBRetriever)
add_subdirectory(NuMIRetriever)
add_subdirectory(EXTRetriever)
add_subdirectory(job)

link_directories(
 $ENV{ROOTSYS}/lib
)

install_headers()
install_fhicl()
install_source()
