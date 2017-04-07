# Using conditional include here because "make distclean" should work, even if
# Makefile.inc doesn't exist (yet).
top_dir         = ../..
-include        $(top_dir)/Makefile.inc

ifeq "$(LIB_TYPE)" "shared"
libDaqDataDecoding = $(top_dir)/src/libDaqDataDecoding.$(ver).so
DLIBS           = -Wl,-rpath,$(top_dir)/src
else
libDaqDataDecoding = $(top_dir)/src/libDaqDataDecoding.$(ver).a
DLIBS           = 
endif
DLIBS          += $(libDaqDataDecoding) $(LIB_EXPAT)

LIBSnoROOT:=$(LIBS)

CXXFLAGS        += -I$(top_dir)/src

CXXFLAGS += `root-config --cflags`
LIBS     += `root-config --libs`

decode_chips    = decode_ChipF1.o decode_ChipADC.o decode_ChipAPV.o decode_ChipGassiplex.o

.PHONY: all
all:            ddd maps plugin.so

ddd:	        ddd.o $(libDaqDataDecoding)
		$(CXX) $(CXXFLAGS) -o $@ ddd.o $(DLIBS) $(LIBS) -lpopt

plugin.so:      plugin.cc
		$(CXX) $(CXXFLAGS)  -shared -o $@ $+ `root-config --cflags --libs` $(DLIBS) $(LIBS)

maps:	        maps.o $(libDaqDataDecoding)
		echo $(LIBSnoROOT)
		$(CXX) $(CXXFLAGS) -o $@ maps.o $(DLIBS) $(LIBSnoROOT) -lpopt

feor:	feor.o $(libDaqDataDecoding)
		$(CXX) $(CXXFLAGS) -o $@ feor.o $(DLIBS) $(LIBS)

fake_event: fake_event.o
		$(CXX) $(CXXFLAGS) -o $@ fake_event.o $(DLIBS) $(LIBS)


.PHONY: clean
clean:
		@rm -rf *.so *.a *.o *.root *.ps *.log log *.o doc Run_* *.err .depend \
		ddd fake_event maps feor


.PHONY: distclean
distclean: clean
		@rm -fr Makefile.inc config.h config.status config.cache

.PHONY: depend
depend:
	@rm -f .depend
	@$(MAKE) .depend

.depend:
	@$(CXX) $(CXXFLAGS) -I. -MM ddd.cc plugin.cc maps.cc feor.cc fake_event.cc > .depend


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include        .depend
endif
endif

