include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)Vars
LIBCXXFILES := $(wildcard *.cxx)

include SoftRelTools/standard.mk

# Need headers only for genie syst enum
override CPPFLAGS += -I$(NUTOOLS_INC) -I$(GENIE_INC)/GENIE/
override CPPFLAGS += -fno-var-tracking-assignments

override LIBLIBS += $(LOADLIBES)  -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) # -lMCReweightFunc

