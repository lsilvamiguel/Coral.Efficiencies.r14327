# Using conditional include here because "make distclean" should work, even if
# Makefile.inc doesn't exist (yet).
top_dir         = ../..
-include $(top_dir)/Makefile.inc


CXXFLAGS       += -I$(top_dir)/src


ifeq "$(LIB_TYPE)" "shared"
libDaqDataDecoding = $(top_dir)/src/libDaqDataDecoding.$(ver).so
MYLIBS         += -Wl,-rpath,$(top_dir)/src $(libDaqDataDecoding)
else
libDaqDataDecoding = $(top_dir)/src/libDaqDataDecoding.$(ver).a
MYLIBS         += $(libDaqDataDecoding)
endif


MYLIBS         += $(LIB_EXPAT)


CXXFLAGS       += `root-config --cflags`
MYLIBS         += `root-config --glibs`


CXXFLAGS       += -I$(DIR_DATE_LIB)/monitoring
ifeq "$(USE_RFIO)" "yes"
MYLIBS         += $(LIB_DATE_LIB) $(LIB_RFIO)
else
MYLIBS         += $(DIR_DATE_LIB)/monitoring/Linux/libmonitorstdalone.a
endif


.PHONY: all
all:            MurphyTV errorDumpAll

MurphyTV:       Main.o CollectErrs.o MurphyTV.o $(libDaqDataDecoding)
	$(CXX) -o MurphyTV Main.o CollectErrs.o MurphyTV.o $(MYLIBS) -lpthread
	@echo "MurphyTV done"
errorDumpAll:   errorDumpAll.o CollectErrs.o $(libDaqDataDecoding)
	$(CXX) -o errorDumpAll errorDumpAll.o CollectErrs.o $(MYLIBS) -lpthread
	@echo "errorDumpAll done"

errorDumpAll.o: errorDumpAll.cxx CollectErrs.h
	$(CXX) $(CXXFLAGS) -c errorDumpAll.cxx

Main.o:         Main.cxx CollectErrs.h MurphyTV.h
	$(CXX) $(CXXFLAGS) -c Main.cxx

MurphyTV.o:     MurphyTV.cxx MurphyTV.h CollectErrs.h helpgui.h
	$(CXX) $(CXXFLAGS) -c MurphyTV.cxx

CollectErrs.o:  CollectErrs.cxx CollectErrs.h
	$(CXX) $(CXXFLAGS) -c CollectErrs.cxx


.PHONY: clean
clean:
	@rm -f *.o MurphyTV errorDumpAll .depend


.PHONY: distclean
distclean: clean


.PHONY: depend
depend:
	@rm -f .depend
	@$(MAKE) .depend


.depend:
	@$(CXX) $(CXXFLAGS) -I. -MM *.cxx > .depend


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include .depend
endif
endif

