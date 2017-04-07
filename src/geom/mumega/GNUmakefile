top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsMumega
TEST_PROGRAMS=

APPL_SRCS = CsMumegaCluster.$(SRC_EXT)          \
            CsMumegaHit.$(SRC_EXT)              \
            CsMumegaPlane.$(SRC_EXT)            \
            CsMumegaPlaneDict.$(SRC_EXT)        \
            CsMumegaTimeCal.$(SRC_EXT)          \
            CsPixelMumegaCluster.$(SRC_EXT)     \
            CsPixelMumegaPlane.$(SRC_EXT)       \
            CsPixelMumegaPlaneDict.$(SRC_EXT)   \
            CsRectPixelMumegaCluster.$(SRC_EXT)     \
            CsRectPixelMumegaPlane.$(SRC_EXT)       \
            CsRectPixelMumegaPlaneDict.$(SRC_EXT)

include $(top_srcdir)/makefiles/General.mk

CsMumegaPlaneDict.$(SRC_EXT): CsMumegaPlane.h CsMumegaPlaneLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c $^

CsPixelMumegaPlaneDict.$(SRC_EXT): CsPixelMumegaPlane.h CsPixelMumegaPlaneLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c $^

CsRectPixelMumegaPlaneDict.$(SRC_EXT): CsRectPixelMumegaPlane.h CsRectPixelMumegaPlaneLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c $^

CLEAN_FILES += *Dict.* *Dict_rdict.pcm

# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): CsMumegaPlaneDict.$(SRC_EXT) CsPixelMumegaPlaneDict.$(SRC_EXT) CsRectPixelMumegaPlaneDict.$(SRC_EXT)

