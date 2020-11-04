FILE(GLOB SR_DEPENDENCIES $ENV{SBNCODE_DIR}/sbncode/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT SRProxy.cxx SRProxy.h FwdDeclare.h
                   COMMAND gen_srproxy -i sbncode/StandardRecord/StandardRecord.h -o SRProxy --target caf::StandardRecord --include-path $ENV{SBNCODE_DIR}:$ENV{ROOT_INC} -op sbncode/CAFAna/StandardRecord/Proxy/ --prolog $ENV{SBNCODE_DIR}/sbncode/CAFAna/StandardRecord/Proxy/Prolog.h --epilog-fwd $ENV{SBNCODE_DIR}/sbncode/CAFAna/StandardRecord/Proxy/EpilogFwd.h
  )

include_directories($ENV{SRPROXY_INC})

# This is a very picky error buried inside template instantiations
#add_definitions(-Wno-int-in-bool-context)

art_make_library(LIBRARY_NAME StandardRecordProxy
                 SOURCE       SRProxy.cxx Instantiations.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

#install(FILES SRProxy.h
#         DESTINATION ${${product}_inc_dir}/CAFAna/StandardRecord/Proxy/)
