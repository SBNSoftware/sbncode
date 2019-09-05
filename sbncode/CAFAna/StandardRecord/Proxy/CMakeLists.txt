add_custom_command(# Rebuild if Event.hh or the script itself change
                   DEPENDS ../../../../sbnanalysis/core/Event.hh
                           parse_xml.py
                   OUTPUT SRProxy.cxx SRProxy.h
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/parse_xml.py
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/sbncode/CAFAna/StandardRecord/Proxy/
                   COMMAND cp SRProxy.h $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/sbncode/CAFAna/StandardRecord/Proxy/
  )

# Allow unqualified StandardRecord/ to find the right headers
include_directories(../..)

# This is a very picky error buried inside template instantiations
add_definitions(-Wno-int-in-bool-context)

art_make_library(LIBRARY_NAME StandardRecordProxy
                 SOURCE       SRProxy.cxx SRProxy.h BasicTypesProxy.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

#install(FILES SRProxy.h
#         DESTINATION ${${product}_inc_dir}/CAFAna/StandardRecord/Proxy/)
