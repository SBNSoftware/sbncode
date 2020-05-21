add_custom_command(# Rebuild if Event.hh or the script itself change
                   DEPENDS ../../../../sbnanalysis/core/Event.hh
                   OUTPUT SRProxy.cxx SRProxy.h FwdDeclare.h
                   COMMAND gen_srproxy -i sbnanalysis/core/Event.hh -o SRProxy --target event::Event --include-path $ENV{SBNCODE_DIR}:$ENV{ROOT_INC} -op sbncode/CAFAna/StandardRecord/Proxy/ --prolog $ENV{SBNCODE_DIR}/sbncode/CAFAna/StandardRecord/Proxy/Prolog.h --epilog-fwd $ENV{SBNCODE_DIR}/sbncode/CAFAna/StandardRecord/Proxy/EpilogFwd.h --extra-cflags ' -DGEN_FLATRECORD_CONTEXT'
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/sbncode/CAFAna/StandardRecord/Proxy/
                   COMMAND cp SRProxy.h FwdDeclare.h $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/sbncode/CAFAna/StandardRecord/Proxy/
  )

include_directories($ENV{SRPROXY_INC})

# This is a very picky error buried inside template instantiations
add_definitions(-Wno-int-in-bool-context)

art_make_library(LIBRARY_NAME StandardRecordProxy
                 SOURCE       SRProxy.cxx Instantiations.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

#install(FILES SRProxy.h
#         DESTINATION ${${product}_inc_dir}/CAFAna/StandardRecord/Proxy/)
