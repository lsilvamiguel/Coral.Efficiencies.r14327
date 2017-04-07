top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsGeomDB
TEST_PROGRAMS=main
TEST_PROGRAMS=testUpdater
TEST_PROGRAMS=testReader

SCHEMA_DDL = \
GdbPerStation.ddl	\

APPL_SRCS = \
CsGeoDbDetector.$(SRC_EXT) \
CsGeoDbUpdater.$(SRC_EXT) \
CsGeoDbReader.$(SRC_EXT) \
GdbDetector.$(SRC_EXT) \
GdbDBuilder.$(SRC_EXT) \
CsCOMGEANT.$(SRC_EXT) \
GdbStation.$(SRC_EXT) \
GdbPerStation_ddl.$(SRC_EXT) \
GdbPerStation.$(SRC_EXT) \
GdbPoint3D.$(SRC_EXT) \
FfrReader.$(SRC_EXT) \
CsRotation.$(SRC_EXT) \
GdbRotation.$(SRC_EXT) \
GdbMaterial.$(SRC_EXT) \
GdbSolidType.$(SRC_EXT) \
GdbVolume.$(SRC_EXT) \
GdbWireInfo.$(SRC_EXT) \
GdbDetectorType.$(SRC_EXT) \
GdbDetTableCont.$(SRC_EXT) \

# Tell to general which extra libs.... and remove them from general...
C_INCLUDES += -I$(top_srcdir)/src/geoutil/g3read

include $(top_srcdir)/makefiles/General.mk
