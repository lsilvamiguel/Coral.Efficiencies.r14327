top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

C_INCLUDES += -I../ntuple/example

APPL_SRCS = main.$(SRC_EXT) \
            DstProd.$(SRC_EXT) \
            DstRead.$(SRC_EXT) \
            deleteDst.$(SRC_EXT) \
            RawRead.$(SRC_EXT) \
            histInfo.$(SRC_EXT) \
            mDstInfo.$(SRC_EXT) \
            RunTest.$(SRC_EXT) \
            getInfo.$(SRC_EXT) \
            WriteRawEvents.$(SRC_EXT) \
            PrintEventList.$(SRC_EXT) \
            prod_control.$(SRC_EXT) \
            removeTBnames.$(SRC_EXT)

ifeq "$(USE_RFIO)" "yes"
APPL_SRCS += oradbTest.$(SRC_EXT)
endif

# not compiling
#            castorFile.$(SRC_EXT) \

include $(top_srcdir)/makefiles/General.mk

all: DstProd  \
     DstRead \
     DstDelete \
     RawRead \
     histInfo \
     mDstInfo \
     runInfo \
     eventInfo \
     WriteRawEvents \
     PrintEventList \
     prodControl \
     removeTBnames

ifeq "$(USE_RFIO)" "yes"
all: oradbTest
endif

# not compiling
#     castorFile \

.PHONY: DstProd
DstProd:        $(OS)/ $(OS)/DstProd
.PHONY: DstRead
DstRead:        $(OS)/ $(OS)/DstRead
.PHONY: DstDelete
DstDelete:      $(OS)/ $(OS)/DstDelete
.PHONY: RawRead
RawRead:        $(OS)/ $(OS)/RawRead
.PHONY: histInfo
histInfo:       $(OS)/ $(OS)/histInfo
.PHONY: mDstInfo
mDstInfo:       $(OS)/ $(OS)/mDstInfo
.PHONY: runInfo
runInfo:        $(OS)/ $(OS)/runInfo
.PHONY: oradbTest
oradbTest:      $(OS)/ $(OS)/oradbTest
.PHONY: eventInfo
eventInfo:      $(OS)/ $(OS)/eventInfo
.PHONY: castorFile
castorFile:     $(OS)/ $(OS)/castorFile
.PHONY: WriteRawEvents
WriteRawEvents: $(OS)/ $(OS)/WriteRawEvents
.PHONY: PrintEventList
PrintEventList: $(OS)/ $(OS)/PrintEventList
.PHONY: prodControl
prodControl:    $(OS)/ $(OS)/prodControl
.PHONY: removeTBnames
removeTBnames:  $(OS)/ $(OS)/removeTBnames

$(OS)/DstProd: $(OS)/main.$(OBJ_EXT) $(OS)/DstProd.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/DstRead: $(OS)/main.$(OBJ_EXT) $(OS)/DstRead.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/DstDelete: $(OS)/deleteDst.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/RawRead: $(OS)/main.$(OBJ_EXT) $(OS)/RawRead.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/histInfo: $(OS)/histInfo.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/mDstInfo: $(OS)/mDstInfo.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/runInfo: $(OS)/RunTest.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/oradbTest: $(OS)/oradbTest.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/eventInfo: $(OS)/getInfo.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/castorFile: $(OS)/castorFile.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif
#	$(C++)  $(C_FLAGS) $(C_INCLUDES) -I$(ORA_SRC) -o  $@  castorFile.cpp \
#	$(ORA_SRC)/$(OS)/RunInfo.$(OBJ_EXT) $(ORA_SRC)/$(OS)/CsOraSession.$(OBJ_EXT)  $(ORA_SRC)/$(OS)/CsStoreMisc.$(OBJ_EXT) \
#	-L${DIR_ORACLE}/lib -locci9 -lclntsh -lshift
##	@scp $(CASTOR_FILE) objsrvvy@lxplus:public/castorFile.exe
##	@scp castorFile.pl objsrvvy@lxplus:public/castorFile.pl

$(OS)/WriteRawEvents: $(OS)/WriteRawEvents.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/PrintEventList: $(OS)/PrintEventList.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/prodControl: $(OS)/prod_control.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

$(OS)/removeTBnames: $(OS)/removeTBnames.$(OBJ_EXT)
	@echo "($(C++)) creating $@"
ifndef VERBOSE
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
else
	$(C++) $(C_FLAGS) $^ $(LIBS) -o $@
endif

