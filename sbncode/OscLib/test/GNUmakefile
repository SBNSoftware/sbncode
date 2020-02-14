#
# $Id: GNUmakefile,v 1.4 2012-08-16 23:10:04 denis Exp $
#
BINS = $(basename $(wildcard *.cc))
SIMPLEBINS = $(BINS)
PACKAGE = OscLib
LIB_TYPE := 
ROOTCINT := 

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk

override CPPFLAGS  += -I$(GSL_INC)
override LDFLAGS   += -L$(GSL_LIB)
override LOADLIBES += -lgsl -lgslcblas
override LIBLIBS   += -L$(GSL_LIB) -lgsl -lgslcblas

BINLIBS += $(ROOTLIBS) -lOscLibFunc -L$(GSL_LIB) -lgsl -lgslcblas 

