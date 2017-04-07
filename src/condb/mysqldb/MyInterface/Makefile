top_srcdir = ../../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY = MySQLDB
TEST_PROGRAMS = Test

APPL_SRCS = MySQLDB.$(SRC_EXT)          \
            MySQLInterface.$(SRC_EXT)   \
            MySQLDBInterface.$(SRC_EXT) \
            MyInterface.$(SRC_EXT)

CLEAN_FILES     += GNUmakefile libMySQLDB.a libMySQLDBInterface.a

include $(top_srcdir)/makefiles/General.mk

all: local_libs

.PHONY: local_libs
local_libs: libMySQLDB.a libMySQLDBInterface.a

GNUmakefile:
	@ln -s Makefile GNUmakefile

libMySQLDB.a: $(OBJS)
		ar cr $@ $^

libMySQLDBInterface.a: $(OS)/MySQLInterface.o $(OS)/MyInterface.o $(OS)/MySQLDBInterface.o
		ar cr $@ $^

