#
#

include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

# SUBDIRS     := Proxy NOvARwgt test

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

LIBLIBS     := -L$(ROOTSYS)/lib -lPhysics

# BINS := hadd_caf
# SCRIPTS := hadd_caf




########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
