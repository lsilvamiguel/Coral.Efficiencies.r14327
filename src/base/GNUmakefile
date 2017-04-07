top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsBase
TEST_PROGRAMS=

APPL_SRCS = Coral.$(SRC_EXT)      \
	    CsInit.$(SRC_EXT)     

F77_SRCS  = inithbook.$(F77_EXT)

include $(top_srcdir)/makefiles/General.mk

all:	CoralCompilationDate.h
CoralCompilationDate.h:	
	(echo -n "#define CORAL_COMPILATION_DATE "; date +%Y%m%d)     > CoralCompilationDate.h 
	(echo -n "#define CORAL_SVN_REVISION \""; svnversion -n ../../; echo "\"" )  >> CoralCompilationDate.h 

CLEAN_FILES += CoralCompilationDate.h
