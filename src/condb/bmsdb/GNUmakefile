top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsBMSDB
TEST_PROGRAMS=main
TEST_PROGRAMS=testUpdater
TEST_PROGRAMS=testReader

##### Add *.ddl file #######
SCHEMA_DDL = \
CsBMSpConstants.ddl	\

APPL_SRCS = \
BmsContainerConf.$(SRC_EXT) \
CsBmsDbUpdater.$(SRC_EXT) \
CsBmsDbReader.$(SRC_EXT) \
CsBMSpConstants.$(SRC_EXT) \
CsBMSpConstants_ddl.$(SRC_EXT) \

include ${top_srcdir}/makefiles/General.mk
