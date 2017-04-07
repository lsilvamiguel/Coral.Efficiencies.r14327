top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

APPL_SRCS = Align.$(SRC_EXT) \
            CheckTracks.$(SRC_EXT) \
            DrawTracks3D.$(SRC_EXT) \
            DrawTracks3DFromTree.$(SRC_EXT) \
            Tracks.$(SRC_EXT) \
            DetFileManager.$(SRC_EXT) \
            DetectorInfo.$(SRC_EXT) \
            DeadZoneInfo.$(SRC_EXT) \
            MagnetInfo.$(SRC_EXT) \
            UpdateFromAlign.$(SRC_EXT) \
            CleanDetFile.$(SRC_EXT) \
            Utils.$(SRC_EXT) \
            Macro.$(SRC_EXT) \
            AddOffset.$(SRC_EXT) \
            Point.$(SRC_EXT) \
            Obj3D.$(SRC_EXT) \
            Opt.$(SRC_EXT) \
            Fit.$(SRC_EXT) \
            TH2Fit.$(SRC_EXT) \
            Dict.$(SRC_EXT) 

F77_SRCS = millepede.$(F77_EXT)

EXTRA_SRCS = Main.align.$(SRC_EXT)           \
             Main.UpdateFromAlign.$(SRC_EXT) \
             Main.CleanDetFile.$(SRC_EXT)    \
             Main.AddOffset.$(SRC_EXT)       \
             Main.CheckTracks.$(SRC_EXT)     \
             Main.traf.$(SRC_EXT)            \
             Main.eff.$(SRC_EXT)             \
             Main.display3D.$(SRC_EXT)       \
             CsManyChunk.$(SRC_EXT)          \
             CsReloadAll.$(SRC_EXT)

OPTIONS_FILES = trafdic.2002.bare.opt \
                trafdic.2003.bare.opt \
                trafdic.2004.bare.opt \
                trafdic.2006.bare.opt \
                trafdic.2007.bare.opt \
                trafdic.2008.minus.bare.opt trafdic.2008.plus.bare.opt \
                trafdic.2009.bare.opt \
                trafdic.h2012.bare.opt

CLEAN_FILES += Dict.* Dict_rdict.pcm
DISTCLEAN_FILES += $(OPTIONS_FILES)

include $(top_srcdir)/makefiles/General.mk

all: $(SHLIB_PREFIX)align.$(SHLIB_EXT) \
     traf \
     display3D \
     align  \
     updateFromAlign  \
     cleanDetFile  \
     checkTracks \
     addOffset \
     ManyChunk \
     eff \
     options_files

$(SHLIB_PREFIX)align.$(SHLIB_EXT):  $(OS)/ $(OS)/$(SHLIB_PREFIX)align.$(SHLIB_EXT)
traf:                               $(OS)/ $(OS)/traf
display3D:                          $(OS)/ $(OS)/display3D 
align:                              $(OS)/ $(OS)/align 
updateFromAlign:                    $(OS)/ $(OS)/updateFromAlign 
cleanDetFile:                       $(OS)/ $(OS)/cleanDetFile 
checkTracks:                        $(OS)/ $(OS)/checkTracks 
addOffset:                          $(OS)/ $(OS)/addOffset 
ManyChunk:                          $(OS)/ $(OS)/ManyChunk
ReloadAll:                          $(OS)/ $(OS)/ReloadAll
eff:                                $(OS)/ $(OS)/eff

##### Options files: C2OPT  (i.e. (.opt.)c file to .opt file)
# Create ``bare'' options files from "../user/trafdic.*.opt.c", to serve as
# a basis (via include statements) for alignment-oriented options files.
options_files: $(OPTIONS_FILES)

# Define some var's used in the cpp compilation, depending upon local system
C2OPT_HEADER:="s%C2OPT_HEADER%// .bare.opt file created from .opt.c counterpart on $(shell date)%"
# Which computer system? => set "C2OPT_COMPUTERFARM"
ifneq (,$(findstring ccali,$(HOST)))    # @ Lyon 
    C2OPT_COMPUTERFARM:=3
endif
ifneq (,$(findstring compass,$(HOST)))  # @ gridKa
    C2OPT_COMPUTERFARM:=2
endif
ifneq (,$(findstring lxplus,$(HOST)))  # @ CERN
    C2OPT_COMPUTERFARM:=1
endif
ifndef C2OPT_COMPUTERFARM
    C2OPT_COMPUTERFARM:=1
endif
C2OPT_FLAG:=0     # To be overriden on command line to enable special settings

# Use sed to:
# - Prepend header
# - Remove #1 directives
# - Change other #line directives to blank lines (This allows to clear up
#  the text)
# - Insert TAB before "//" introducing ``in-line'' comment
%.bare.opt	:  ../user/%.opt.c GNUmakefile
	@echo $(CC) $< ">" $@
	@$(CC) -E -ffreestanding -traditional-cpp -C \
		-DC2OPT_COMPUTERFARM=$(C2OPT_COMPUTERFARM) \
		-DC2OPT_FLAG=$(C2OPT_FLAG) \
		-DC2OPT_BARE_FILE \
		$< | \
	sed -e "s%C2OPT_HEADER%// \"$@\" file created from \"$<\", on $(shell date)%" \
	    -e "/^# 1/ d" \
	    -e "/^#/ s/.*//" \
	    -e "/^\/\//! s%\([^/ ] \)//%\1	//%" \
	> $@

Dict.$(SRC_EXT): Align.h \
                 Tracks.h \
                 DetFileManager.h \
                 DetectorInfo.h \
                 DeadZoneInfo.h \
                 MagnetInfo.h \
                 CheckTracks.h \
                 DrawTracks3D.h \
                 Macro.h \
                 Fit.h \
                 TH2Fit.h \
                 Utils.h \
                 LinkDef.h
	@echo creating ROOT Dictionary $@
	@$(ROOTSYS)/bin/rootcint -f $@ -c $^

### Tree and Detector Managment library
$(OS)/$(SHLIB_PREFIX)align.$(SHLIB_EXT): $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) -shared -Wl,-symbolic -o $@ $^ $(C_FLAGS) $(LIBS)

### Executable    
$(OS)/align: $(OS)/Main.align.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/updateFromAlign: $(OS)/Main.UpdateFromAlign.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/cleanDetFile: $(OS)/Main.CleanDetFile.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/addOffset: $(OS)/Main.AddOffset.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/checkTracks: $(OS)/Main.CheckTracks.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/traf: $(OS)/Main.traf.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/eff: $(OS)/Main.eff.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/display3D: $(OS)/Main.display3D.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

### Utilities
$(OS)/ManyChunk: $(OS)/CsManyChunk.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/ReloadAll: $(OS)/CsReloadAll.o $(APPL_OBJS)
	@echo "($(C++)) creating $@" 
	@$(C++) $(C_FLAGS) $^ $(LIBS) -o $@

$(OS)/%.$(OBJ_EXT): %.$(SRC_EXT) $(SCHEMA_HDRS) | $(ODIRS)
	@echo "(C++) compiling $*"
ifndef VERBOSE
	@$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
else
	$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
endif


# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): Dict.$(SRC_EXT)

