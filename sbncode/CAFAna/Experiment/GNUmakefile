include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)Experiment
LIBCXXFILES := $(wildcard *.cxx)

include SoftRelTools/standard.mk

override LIBLIBS += $(LOADLIBES)  -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)

