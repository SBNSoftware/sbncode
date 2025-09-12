include_directories ( $ENV{GEANT4REWEIGHT_INC} )
link_directories ( $ENV{GEANT4REWEIGHT_LIB} )

add_subdirectory(Base)
add_subdirectory(Calculators)
add_subdirectory(App)
add_subdirectory(jobs)

