top_dir       = ..
-include $(top_dir)/Makefile.inc

INCLUDES += -I../src

ifeq ($(USE_CERNLIB),yes)
SRC = db.cc gui.cc testreco.cc testRCecorr.cc testRCpcorr.cc ToyCalorimeter.cc
else
SRC = db.cc
endif

OBJ = $(SRC:.cc=.o)

##########################################################################################

.PHONY: all
#all: 		r rec calibration online gui
ifeq ($(USE_CERNLIB),yes)
all: db gui testreco testRCecorr testRCpcorr testEdepCorr testTiSdepCorr
else
all: db
endif

db: db.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

r: r.o CalorimeterGAMS.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

rec: main.o ftest.o CalorimeterGAMS.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

gui: gui.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

calibration: calibration.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

decoding_dir = ../../../DaqDataDecoding
date_dir     = /afs/cern.ch/compass/online/daq/dateV37
online: $(top_dir)/src/$(libReco) CalorimeterGAMS.o online.cc
	$(CXX) $(INCLUDES) $(CXXFLAGS) -Dlong32=int -I$(decoding_dir) -I$(date_dir) \
	-o online online.cc CalorimeterGAMS.o \
	-Xlinker -rpath -Xlinker $(decoding_dir) $(decoding_dir)/libDaqDataDecoding.a \
	../src/libReco.a $(date_dir)/monitoring/Linux/libmonitor.a $(CERN_LIB) $(LIBS) -lg2c -lnsl

testreco: testreco.o ToyCalorimeter.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

testRCecorr: testRCecorr.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

testRCpcorr: testRCpcorr.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

testEdepCorr: testEdepCorr.o ToyCalorimeter.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

testTiSdepCorr: testTiSdepCorr.o ToyCalorimeter.o $(top_dir)/src/$(libReco)
ifneq "$(VERBOSE)" ""
	$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
else
	@echo Linking $@
	@$(CXX) -o $@ $+ $(LD_OPT) $(CXXFLAGS)
endif

.PHONY: clean
clean:
	@rm -fr $(clean_mask) $(OBJ) t1 t2 rec db gui reco.db gui-DB calibration online testreco testEdepCorr testRCecorr testRCpcorr testTiSdepCorr .depend

.PHONY: distclean
distclean: clean

.PHONY: depend
depend:
	@rm -f .depend
	@$(MAKE) .depend

.depend:
	@$(CXX) -MM $(CXXFLAGS) $(INCLUDES) $(SRC) > .depend

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include .depend
endif
endif

