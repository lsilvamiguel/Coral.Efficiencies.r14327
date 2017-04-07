top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsGEM
TEST_PROGRAMS=

APPL_SRCS = CsGEMCluster.$(SRC_EXT)         \
            CsGEMHit.$(SRC_EXT)             \
            CsGEMPlane.$(SRC_EXT)           \
            CsGEMPlaneDict.$(SRC_EXT)       \
            CsGEMTimeCal.$(SRC_EXT)         \
            CsPixelGEMCluster.$(SRC_EXT)    \
            CsPixelGEMPlane.$(SRC_EXT)      \
            CsPixelGEMPlaneDict.$(SRC_EXT)

include $(top_srcdir)/makefiles/General.mk

CsGEMPlaneDict.$(SRC_EXT): CsGEMChan.h CsGEMCluster.h CsGEMHit.h CsGEMPlane.h CsGEMPlaneHeader.h CsGEMPlanePar.h CsGEMTimeCal.h CsGEMPlaneLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c $^

CsPixelGEMPlaneDict.$(SRC_EXT): CsPixelGEMCluster.h CsPixelGEMPlane.h CsPixelGEMPlanePar.h CsPixelGEMPlaneLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c $^

CLEAN_FILES += *Dict.* *Dict_rdict.pcm

# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): CsGEMPlaneDict.$(SRC_EXT) CsPixelGEMPlaneDict.$(SRC_EXT)
