add_custom_command(# Rebuild if Event.hh or the script itself change
                   DEPENDS $ENV{SBNCODE_DIR}/sbncode/StandardRecord/StandardRecord.h
                           ${CMAKE_CURRENT_SOURCE_DIR}/gen_flatrecord.py
                   OUTPUT FlatRecord.cxx
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_flatrecord.py
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/FlatMaker/
                   COMMAND cp Flat*.h $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/FlatMaker/
  )

art_make_library(LIBRARY_NAME FlatMaker
                 SOURCE       FlatRecord.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER} # CAFAnaCore
                 )

art_make_exec(NAME convert_to_flat
              SOURCE convert_to_flat.cc
              LIBRARIES FlatMaker
# sbnanalysis doesn't participate fully in the mrb build. Skip this for now
#                        sbnanalysis_Event
                        ${ROOT_BASIC_LIB_LIST})
