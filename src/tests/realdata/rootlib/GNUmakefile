DEPEND_FILE = depend.mk
SRCS = $(wildcard *.cc)
OS = Linux
C++ = gcc
C_INCLUDES = -I$(ROOTSYS)/include -I$(CORAL)/include
C_FLAGS = `$(ROOTSYS)/bin/root-config --libs` -lMinuit

all: $(OS) $(OS)/DBShift \
           $(OS)/DBScale \
           libHManager.so
        
###
$(OS)/DBShift:   $(OS)/DBShift.o
	@echo creating $@
	@$(C++)  $(C_FLAGS) $(C_INCLUDES) $^ -o $@

###
$(OS)/DBScale:   $(OS)/DBScale.o
	@echo creating $@
	@$(C++)  $(C_FLAGS) $(C_INCLUDES) $^ -o $@

### HMANAGER To handle coral root files ###
libHManager.so: $(OS)/DetInfo.o \
                $(OS)/DetFileInfo.o \
                $(OS)/RTInfo.o \
                $(OS)/HManager.o \
                $(OS)/RTManager.o \
                $(OS)/ResManager.o \
                $(OS)/EffManager.o \
                $(OS)/TH2Fitter.o \
                $(OS)/Fit.o \
                $(OS)/Macro.o \
                $(OS)/Utils.o \
                $(OS)/HDict.o 
	@echo creating library $@
	@$(C++) -shared $(C_FLAGS) $(C_INCLUDES) $^ -o $@
  
HDict.h:        HDict.cc
HDict.cc:       HManager.h \
                RTManager.h \
                ResManager.h \
                EffManager.h \
                DetInfo.h \
                DetFileInfo.h \
                RTInfo.h \
                TH2Fitter.h \
                Fit.h \
                Utils.h \
                Macro.h\
                $(CORAL)/include/CsRTGridPoint.h\
                LinkDef.h
	@echo creating ROOT Dictionary $@
	@$(ROOTSYS)/bin/rootcint -f HDict.cc -c $(CFLAGS) $(C_INCLUDES) -I$(CORAL)/include/ $^


depend $(DEPEND_FILE) : 
	@echo creating dependencies
	@$(C++) $(C_INCLUDES) -MM $(SRCS) \
  | sed -e 's|^\(.*\.o:\)|$(OS)/\1|'>$(DEPEND_FILE)

clean: 
	@echo cleaning
	@rm -rf $(OS) libHManager.so \
                libDFManager.so \
                HDict.cc \
                HDict.h  \
                DFDict.cc \
                DFDict.h 

# rules
$(OS):
	@echo "creating platform directory for system $(OS)"
	@mkdir $(OS)
  
$(OS)/%.o: %.cc 
	@echo "($(C++)) compiling $*"
	@$(C++) -g $(C_INCLUDES) -c $*.cc -o $(OS)/$*.o
  
$(OS)/%.o: %.C 
	@echo "($(C++)) compiling $*"
	@$(C++) -g $(C_INCLUDES) -c $*.C -o $(OS)/$*.o

include $(DEPEND_FILE)
