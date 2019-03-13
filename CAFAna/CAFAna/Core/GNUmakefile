include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhc.mk
include SoftRelTools/arch_spec_boost.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)Core
LIBCXXFILES := $(wildcard *.cxx)

include SoftRelTools/standard.mk

override LIBLIBS += $(LOADLIBES)  -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(IFDHC_FQ_DIR)/lib/ -lifdh -lOscLibFunc

