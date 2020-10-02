
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)
# SUBDIRS     := scripts

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lart_Framework_Services_Registry -lPhysics -lart_Persistency_Common -lart_Framework_Core -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_nutools.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_ifdhc.mk
include SoftRelTools/arch_spec_ifdhart.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(NUTOOLS_LIB) -lStandardRecord 
