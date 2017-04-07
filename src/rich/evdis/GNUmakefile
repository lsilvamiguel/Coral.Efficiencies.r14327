top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsRCEvdis
TEST_PROGRAMS=
#VERBOSE=true

SRCS      = CsRCGHit.$(SRC_EXT)          \
	    CsRCGCathode.$(SRC_EXT)      \
	    CsRCGRing.$(SRC_EXT)         \
            CsRCGTrack.$(SRC_EXT)        \
            CsRCGRingTrack.$(SRC_EXT)    \
            CsRCEDisplay.$(SRC_EXT)      \
	    CsRCGraph.$(SRC_EXT)         \
            CsRCIFtoDisplay.$(SRC_EXT)	

SRCS_CINT = CsRichOneDisplay.$(SRC_EXT) CsRCDisplay.$(SRC_EXT)

SRCS_ROOT = CsRCEvdisDict.$(SRC_EXT)

APPL_SRCS = $(SRCS_CINT) $(SRCS_ROOT) $(SRCS) 

HEAD_EXT  = h
HEAD  = ${SRCS_CINT:%.$(SRC_EXT)=%.$(HEAD_EXT)}
#HEAD = CsRCDisplay.$(HEAD_EXT) 

CLEAN_FILES = *Dict.* *Dict_rdict.pcm

# Tell to general which extra libs.... and remove them from general...

#ifeq "$(OS)" "Linux"
#C_FLAGS += -fno-rtti -fPIC
#endif

# General... do not remove...

include $(top_srcdir)/makefiles/General.mk

# This is particular for ROOT...

$(SRCS_ROOT): $(HEAD) LinkDef.h
	@$(ROOTSYS)/bin/rootcint -v4 -f $(SRCS_ROOT) -c -I`root-config --incdir` -I$(top_srcdir)/include $(HEAD)  LinkDef.h

# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): $(SRCS_ROOT)

