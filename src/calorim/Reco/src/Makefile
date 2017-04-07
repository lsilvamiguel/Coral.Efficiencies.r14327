-include ../Makefile.inc

SRC +=                                      \
Calibration.cc                              \
Calorimeter.cc                              \
CalorimeterAlignmentInTime.cc               \
CalorimeterCellsMatrix.cc                   \
CalorimeterDB.cc                            \
CalorimeterOptions.cc                       \
CalorimeterParticle.cc                      \
Cell.cc                                     \
CellDataRaw.cc                              \
CellType.cc                                 \
Cluster.cc                                  \
DataBase.cc                                 \
DataMapping.cc                              \
GUICalorimeter.cc                           \
GUICellsMatrix.cc                           \
MCConstruction.cc                           \
myROOT_utils.cc                             \
Noise.cc                                    \
OneParticleResponse.cc                      \
Reconstruction.cc                           \
ReconstructionCombined.cc                   \
ReconstructionCombinedCorrection.cc         \
ReconstructionCombinedFitInput.cc           \
ReconstructionKolosov.cc                    \
ReconstructionKolosovFindClusters.cc        \
ReconstructionKolosovFindGamma.cc           \
ReconstructionKolosovFit.cc                 \
ReconstructionLednev.cc                     \
ReconstructionTest.cc                       \
Shower.cc                                   \
ShowerProfile.cc                            \
ShowerProfileLednev.cc                      \
StatInfo.cc                                 \
SubSet.cc

ifeq "$(USE_Qt)" "yes"
CREATE_FILES += GUIBaseCalorimeter_moc.cc   \
                GUICalorimeter_moc.cc       \
                GUIBaseCalorimeter_uic.h    \
                GUICalorimeterMove_uic.h
SRC += $(filter %.cc, $(CREATE_FILES))      \
       GUIBaseCalorimeter.cc                \
       GUICalorimeterMove.cc
endif

OBJS += $(SRC:.cc=.o)

FINC += -Ifortran_reco/commonblocks/ -Ifortran_reco/functions/

FSRC  = fortran_reco/separa.F fortran_reco/badcheck.F fortran_reco/ecal2reco.F \
        fortran_reco/fastreco.F fortran_reco/Ecell.F fortran_reco/showerfit.F

OBJS += $(FSRC:.F=.o)

clean_mask    = .depend tmp_dep *.o fortran_reco/*.o *.a *.so \
                $(CREATE_FILES) \
                *.bck *.bak *~ *.root *.log *.dvi *.aux fort.* *Cint.[ch] *Cint.cxx \
                *.out *.gprof gdb-transcript core doc html

distclean_mask = Reco_config.h


##########################################################################################

.PHONY: all
ifeq "$(LIB_TYPE)" "shared"
  all:          libReco.$(ver).so

  libReco.$(ver).so: $(OBJS)
		$(CXX) $(CXXFLAGS) -shared -o $@ $(OBJS)
		ln -fs libReco.$(ver).so libReco.so
else
  all:          libReco.$(ver).a

  libReco.$(ver).a: $(OBJS)
		ar cr $@ $(OBJS)
		ln -fs libReco.$(ver).a libReco.a
endif

Cint.cc: $(Cint)
		rootcint -f Cint.cxx -c `root-config --cflags` $(Cint)
		@mv Cint.cxx Cint.cc

ifeq "$(USE_Qt)" "yes"
GUIBaseCalorimeter_uic.h: GUIBaseCalorimeter.ui
		$(UIC) -o GUIBaseCalorimeter_uic.h GUIBaseCalorimeter.ui

GUICalorimeterMove_uic.h: GUICalorimeterMove.ui
		$(UIC) -o GUICalorimeterMove_uic.h GUICalorimeterMove.ui

GUIBaseCalorimeter_moc.cc: GUIBaseCalorimeter.h
		$(MOC) -o GUIBaseCalorimeter_moc.cc GUIBaseCalorimeter.h

GUICalorimeter_moc.cc: GUICalorimeter.h
		$(MOC) -o GUICalorimeter_moc.cc GUICalorimeter.h
endif

.PHONY:
clean:
		@rm -fr $(clean_mask)

.PHONY:
distclean: clean
		@rm -fr $(distclean_mask)

.PHONY: depend
depend:
	rm -f .depend
	$(MAKE) .depend

.depend: $(CREATE_FILES)
	@i=0; $(CXX) -MM $(CXXFLAGS) $(INCLUDES) $(SRC) > tmp_dep || i=1 ; $(CXX) -MM $(FOPT) $(FINC) $(FSRC) >> tmp_dep || i=1 ; if test $$i -eq 1 ; then rm -f tmp_dep ; exit 1 ; fi
# post-process the depend file, entries for fortran files are created
# with file.o instead of fortran_reco/file.o
	@for i in $(FSRC:.F=.o) ; do sed -i -e "s*^`basename $$i`*$$i*g" tmp_dep ; done
# avoid stale .depend file in case of error
	@mv tmp_dep .depend

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include .depend
endif
endif
