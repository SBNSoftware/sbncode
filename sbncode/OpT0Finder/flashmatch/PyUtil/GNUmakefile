#
# This is an example GNUmakefile for my packages
#
PACKAGE_NAME = PyUtil

# specific names for this package
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
HEADERS_DEST = flashmatch/core/$(PACKAGE_NAME)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell flashmatch-config --includes)
INCFLAGS += $(shell ${FMATCH_PYTHON_CONFIG} --includes)
INCFLAGS += -I$(shell ${FMATCH_PYTHON} -c "import numpy; print(numpy.get_include())")

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(FMATCH_BASEDIR)/Makefile/Makefile.${OSNAME}

LDFLAGS += -L$(shell ${FMATCH_PYTHON_CONFIG} --prefix)/lib/
LDFLAGS += $(shell ${FMATCH_PYTHON_CONFIG} --ldflags)

# call the common GNUmakefile
include $(FMATCH_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build: # nothing extra for this package
pkg_clean: # nothing extra for this package
