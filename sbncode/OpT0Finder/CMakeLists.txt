# Add Geant4 directories as needed by flashmatch code
#include_directories($ENV{GEANT4_FQ_DIR}/include)

add_subdirectory(Products)
add_subdirectory(flashmatch)
add_subdirectory(job)
