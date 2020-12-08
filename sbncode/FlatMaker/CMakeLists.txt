FILE(GLOB SR_DEPENDENCIES $ENV{SBNCODE_DIR}/sbncode/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT FlatRecord.cxx FlatRecord.h FwdDeclare.h
                   COMMAND gen_srproxy --flat -i sbncode/StandardRecord/StandardRecord.h -o FlatRecord --target caf::StandardRecord --include-path $ENV{SBNCODE_DIR}:$ENV{ROOT_INC} -op sbncode/FlatMaker --prolog $ENV{SBNCODE_DIR}/sbncode/FlatMaker/Prolog.h --extra-cflags ' -D_Float16=short'
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/FlatMaker/
                   COMMAND cp FlatRecord.h FwdDeclare.h $ENV{MRB_INSTALL}/sbncode/$ENV{SBNCODE_VERSION}/include/FlatMaker/
  )

include_directories($ENV{SRPROXY_INC})

art_make_library(LIBRARY_NAME FlatMaker
                 SOURCE       FlatRecord.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST}
                              ${ROOT_TREEPLAYER}
                              -L$ENV{CAFANACORE_LIB}
                              CAFAnaCoreExt
                 )

art_make_exec(NAME flatten_caf
              SOURCE flatten_caf.cc
              LIBRARIES FlatMaker CAFAnaCore
                        ${ROOT_BASIC_LIB_LIST})
