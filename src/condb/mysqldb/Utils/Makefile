top_srcdir=../../../..
include $(top_srcdir)/makefiles/Include.mk

APPL_SRCS = loadOldFileInDB.$(SRC_EXT) \
            getDBFilePath.$(SRC_EXT) \
            addFileInDB.$(SRC_EXT) \
            get1stEvtPath.$(SRC_EXT) \
            getTgtCurrents.$(SRC_EXT) \
            Utilities.$(SRC_EXT) \
            SyncDB.$(SRC_EXT)

EXECS = loadOldFileInDB \
        getDBFilePath \
        addFileInDB \
        get1stEvtPath \
        getTgtCurrents

CLEAN_FILES += $(EXECS)

include $(top_srcdir)/makefiles/General.mk

all :  $(EXECS)

SyncDB :	$(OS)/SyncDB.o $(OS)/Utilities.o
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

loadOldFileInDB :	$(OS)/loadOldFileInDB.o
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

addFileInDB :	$(OS)/addFileInDB.o
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

getDBFilePath :	$(OS)/getDBFilePath.o ../MyInterface/libMySQLDBInterface.a
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

get1stEvtPath :	$(OS)/get1stEvtPath.o ../MyInterface/libMySQLDBInterface.a
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

getTgtCurrents :	$(OS)/getTgtCurrents.o ../MyInterface/libMySQLDBInterface.a
		$(CXX) $(C_FLAGS) $(C_INCLUDES) -o $@ $^ $(LIBS)

