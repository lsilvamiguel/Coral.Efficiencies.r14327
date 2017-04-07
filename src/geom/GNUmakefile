top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsGeom
TEST_PROGRAMS=

APPL_SRCS = CsDet.$(SRC_EXT)                        \
			CsDetector.$(SRC_EXT)                   \
			CsDWDetector.$(SRC_EXT)                 \
			CsCEDARDetector.$(SRC_EXT)              \
			CsStrawTubesDetector.$(SRC_EXT)         \
			CsStrawTubesDetector-calib.$(SRC_EXT)   \
			CsDriftChamberDetector.$(SRC_EXT)       \
			CsMWPCDetector.$(SRC_EXT)               \
			CsMW1Detector.$(SRC_EXT)                \
			CsGEMDetector.$(SRC_EXT)                \
			CsPixelGEMDetector.$(SRC_EXT)           \
			CsMicroMegaDetector.$(SRC_EXT)          \
                        CsPixelMumegaDetector.$(SRC_EXT)        \
			MMLibrary.$(SRC_EXT)                    \
			CsSiTrackerDetector.$(SRC_EXT)          \
			CsDFiberHodoDetector.$(SRC_EXT)         \
			CsJFiberHodoDetector.$(SRC_EXT)         \
			CsTriggerHodoDetector.$(SRC_EXT)        \
			CsBMSDetector.$(SRC_EXT)                \
			CsRICH1Detector.$(SRC_EXT)              \
			CsDriftTubeDetector.$(SRC_EXT)          \
			CsRichWallDetector.$(SRC_EXT)           \
			CsRwRecons.$(SRC_EXT)                   \
			CsRwChargeRecons.$(SRC_EXT)             \
			CsMiscDetector.$(SRC_EXT)               \
			CsRPDetector.$(SRC_EXT)                 \
			CsField.$(SRC_EXT)                      \
			CsZone.$(SRC_EXT)                       \
			CsGeom.$(SRC_EXT)                       \
			CsMaterialMap.$(SRC_EXT)                \
			CsRTRelation.$(SRC_EXT)                 \
			CsRICH1UpGrade.$(SRC_EXT)               \
			CathodeAPV.$(SRC_EXT)                   \
			CathodeMAPMT.$(SRC_EXT)                 \
			RayTrace.$(SRC_EXT)

include $(top_srcdir)/makefiles/General.mk

CLEAN_FILES += *Dict.*

ifneq ("$(COMPASS_SETUP)","")
C_FLAGS+=-DCOMPASS_SETUP=$(COMPASS_SETUP)
endif

