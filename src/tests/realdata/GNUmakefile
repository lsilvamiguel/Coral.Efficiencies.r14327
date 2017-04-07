top_srcdir  = ../../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY = HManager

APPL_SRCS = HManager.$(SRC_EXT) \
            EffManager.$(SRC_EXT) \
            ResManager.$(SRC_EXT) \
            RTManager.$(SRC_EXT) \
            Obj3D.$(SRC_EXT) \
            Point.$(SRC_EXT) \
            Utils.$(SRC_EXT) \
            Fit.$(SRC_EXT) \
            TH2Fit.$(SRC_EXT) \
            DetFileManager.$(SRC_EXT) \
            DetectorInfo.$(SRC_EXT) \
            DeadZoneInfo.$(SRC_EXT) \
            MagnetInfo.$(SRC_EXT) \
            RTInfo.$(SRC_EXT) \
            Dict.$(SRC_EXT) 

CLEAN_FILES = *Dict.* *Dict_rdict.pcm

LIB_TYPE = shared

include $(top_srcdir)/makefiles/General.mk

### all 
all: spoint \
     traf \
     DetOff \
     ManyChunk \
     ReloadAll
spoint:    $(OS)/ $(OS)/spoint 
traf:      $(OS)/ $(OS)/traf
DetOff:    $(OS)/ $(OS)/DetOff
ManyChunk: $(OS)/ $(OS)/ManyChunk
ReloadAll: $(OS)/ $(OS)/ReloadAll

# executables
$(OS)/DetOff: $(OS)/CsDetOff.o $(OS)/Opt.o $(APPL_OBJS)
	@echo creating $@
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/ManyChunk: $(OS)/CsManyChunk.o $(OS)/Opt.o $(APPL_OBJS)
	@echo creating $@
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/ReloadAll: $(OS)/CsReloadAll.o $(APPL_OBJS)
	@echo creating $@
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/spoint: $(OS)/Main.spoint.o $(APPL_OBJS)
	@echo creating $@
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/traf: $(OS)/Main.traf.o $(APPL_OBJS)
	@echo creating $@
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

Dict.$(SRC_EXT): HManager.h \
         		 EffManager.h \
		         ResManager.h \
		         RTManager.h \
        		 DetFileManager.h \
		         DetectorInfo.h \
        		 DeadZoneInfo.h \
		         MagnetInfo.h \
        		 RTInfo.h \
		         Utils.h \
        		 Fit.h \
		         TH2Fit.h \
        		 LinkDef.h
	@echo creating ROOT Dictionary $@
	@$(ROOTSYS)/bin/rootcint -f $@ -c $^

#
# C++ Compilation
#
$(OS)/%.$(OBJ_EXT): %.$(SRC_EXT) $(SCHEMA_HDRS) | $(ODIRS)
	@echo "(C++) compiling $*"
ifndef VERBOSE
	@$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
else
	$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
endif

# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): Dict.$(SRC_EXT)

