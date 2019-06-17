PACKAGE     := CAFAna
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
SUBDIRS     := Core Vars Cuts Decomp Experiment Extrap Prediction Analysis Systs # Vars Cuts Systs Unfold Decomp Experiment Extrap FC Prediction Analysis XSec grid
LIBCXXFILES := $(wildcard *.cxx)


BINS := cafe
SCRIPTS := cafe

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhc.mk

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE) -lCAFAnaCore -lCAFAnaVars -lCAFAnaCuts -lCAFAnaSysts -lCAFAnaUnfold -lCAFAnaDecomp -lCAFAnaExperiment -lCAFAnaExtrap -lCAFAnaFC -lCAFAnaPrediction -lCAFAnaAnalysis -lCAFAnaXSec

# Explicitly link in root libraries. Why is this necessary? Shouldn't
# arch_spec_root.mk do that? Nonetheless, without this we get undefined symbol
# errors...
override LIBLIBS += `root-config --libs` -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(IFDHC_FQ_DIR)/lib/ -lifdh
