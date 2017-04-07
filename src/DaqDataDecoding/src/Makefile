# Using conditional include here because "make distclean" should work, even if
# Makefile.inc doesn't exist (yet).
-include        ../Makefile.inc

OBJS            = ObjectXML.o DateEquipment.o Exception.o DaqOption.o DaqError.o DaqMap.o \
                  DaqEvent.o Event1Run.o CatchInfo.o DaqEventsManager.o SLink.o \
                  Chip.o ChipF1.o ChipSinica.o ChipCol.o ChipADC.o ChipSADC.o ChipSADC_huffman.o \
                  ChipAPV.o ChipAPVRICH.o ChipGandalf.o \
                  ChipHotGeSiCA.o ChipGassiplex.o Scaler.o PixelMMDecoding.o \
                  Trigger.o TriggerTime.o OnlineFilter.o utils.o Stat.o

ifeq "$(LIB_TYPE)" "shared"
  all:          libDaqDataDecoding.$(ver).so
  libDaqDataDecoding.$(ver).so:      $(OBJS)
		$(CXX) -shared -o $@ $(OBJS)
		ln -fs libDaqDataDecoding.$(ver).so libDaqDataDecoding.so
else
  all:          libDaqDataDecoding.$(ver).a
  libDaqDataDecoding.$(ver).a:  $(OBJS)
		ar cr $@ $(OBJS)
		ln -fs libDaqDataDecoding.$(ver).a libDaqDataDecoding.a
endif

.PHONY: clean
clean:
		@rm -rf *.so *.a *.o *.root *.ps *.log log doc .depend

.PHONY: distclean
distclean: clean
		@rm -fr config.h


.PHONY: depend
depend:
	@rm -f .depend
	@$(MAKE) .depend


.depend:
	@$(CXX) $(CXXFLAGS) -I. -MM *.cc > .depend


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include        .depend
endif
endif
