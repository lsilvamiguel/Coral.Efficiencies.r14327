top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CondDB
TEST_PROGRAMS=

APPL_SRCS = CondDbHandler.$(SRC_EXT)

include $(top_srcdir)/makefiles/General.mk

C_INCLUDES +=-I/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/HepUtilities/0.5.0.0/include
LIBS +=-L/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/HepUtilities/0.5.0.0/lib -lHepUtilities

