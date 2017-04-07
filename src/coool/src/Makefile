include         ../Makefile.inc

OBJS            = 	$(OS)/StartFrame.o $(OS)/MainFrame.o \
			$(OS)/VariousSettings.o $(OS)/Geometry.o\
			$(OS)/PlanePanel.o $(OS)/Plane1VPanel.o \
			$(OS)/Plane2VPanel.o  $(OS)/PlaneRPD_SADC_Panel.o \
			$(OS)/PlaneDriftChamberPanel.o $(OS)/PlaneRPD_F1_Panel.o \
			$(OS)/MonitorPanel.o \
			$(OS)/Monitor.o $(OS)/Variable.o \
			$(OS)/GroupPanel.o $(OS)/GroupScifiJPanel.o \
			$(OS)/GroupMumegaPanel.o \
			$(OS)/GroupMuonWallAPanel.o \
			$(OS)/Plane.o $(OS)/Plane1V.o $(OS)/PlaneMumega.o \
			$(OS)/PlanePMumega.o $(OS)/PlaneDriftChamber2V.o \
			$(OS)/PlaneDriftChamber.o $(OS)/PlaneScifiJ.o \
			$(OS)/PlaneTrigHodo.o $(OS)/PlaneMISC.o \
			$(OS)/PlaneCEDAR.o \
			$(OS)/PlanePrimakoffHodo.o \
			$(OS)/PlaneBeamKiller.o \
			$(OS)/PlaneVBHodo.o \
			$(OS)/PlaneVBOX.o \
			$(OS)/PlaneVeto.o \
			$(OS)/PlaneMwpc.o $(OS)/PlaneScifiG.o \
			$(OS)/PlaneBMS.o $(OS)/PlaneRICH.o \
			$(OS)/PlaneAPV.o $(OS)/PlaneGEM.o \
                        $(OS)/PlanePGEM.o $(OS)/PlaneSili.o \
			$(OS)/PlaneRiAPV.o $(OS)/PlaneRiAPVPanel.o \
			$(OS)/PlaneHCAL1.o $(OS)/PlaneHCAL2.o \
			$(OS)/PlaneHCALT.o \
                        $(OS)/PlaneECAL0.o \
			$(OS)/PlaneHCAL2Sum.o \
			$(OS)/PlaneECAL1.o $(OS)/PlaneECAL1Panel.o \
			$(OS)/PlaneFEM.o $(OS)/PlaneFEMPanel.o \
			$(OS)/PlaneECAL2.o $(OS)/PlaneECAL2Panel.o \
                        $(OS)/PlaneMuonWallA.o \
			$(OS)/PlaneMuonWallB.o $(OS)/PlaneScaler.o \
			$(OS)/Plane2V.o $(OS)/PlaneRCath.o $(OS)/PlaneRICH_MAPMT.o\
			$(OS)/PlaneRPD_SADC.o $(OS)/PlaneRPD_F1.o $(OS)/GroupRPD.o \
			$(OS)/PlaneSandwich.o \
			$(OS)/PlaneSandwichPanel.o \
			$(OS)/GeomPlane.o $(OS)/GeomPlaneMumega.o $(OS)/GeomPlanePMumega.o \
			$(OS)/PlaneTrigger_SADC_F1.o \
			$(OS)/GeomPlaneMwpc.o  $(OS)/GeomPlaneGEM.o $(OS)/GeomPlanePGEM.o \
			$(OS)/GeomPlaneSili.o $(OS)/GeomPlaneScifiJ.o \
			$(OS)/GeomPlaneScifiG.o $(OS)/GeomPlaneMuonWallA.o \
			$(OS)/PlaneStrawTubes.o  $(OS)/GroupStraws.o $(OS)/Group.o \
			$(OS)/GroupScifiJ.o $(OS)/GroupTrigHodo.o\
			$(OS)/GroupDAQ.o $(OS)/GroupMwpc.o \
			$(OS)/GroupScifiG.o $(OS)/GroupMumega.o \
			$(OS)/GroupGEM.o $(OS)/GroupSili.o \
			$(OS)/GroupPGEM.o \
			$(OS)/GroupMuonWallA.o $(OS)/BitButton.o \
			$(OS)/GroupBeamKiller.o \
			$(OS)/GroupCEDAR.o \
			$(OS)/GroupRiAPV.o \
			$(OS)/GroupECAL1.o \
			$(OS)/GroupRICH_MAPMT.o \
			$(OS)/Reference.o $(OS)/MwpcReconst.o \
			$(OS)/MwpcEventDisplay.o \
			$(OS)/PlaneRW.o \
			$(OS)/PlaneTcsPhase.o \
			$(OS)/GroupRW.o \
			$(OS)/PlaneCamera.o $(OS)/GroupCamera.o \
			$(OS)/PlaneGandalf.o \
			$(OS)/Dict.o

dict_headers     =  	Plane.h Plane1V.h PlaneMumega.h  \
			PlanePMumega.h PlaneDriftChamber2V.h \
			PlaneDriftChamber.h PlaneScifiJ.h \
			PlaneTrigHodo.h PlaneMISC.h PlaneVeto.h \
			PlaneCEDAR.h \
			PlanePrimakoffHodo.h \
			PlaneBeamKiller.h \
			PlaneVBHodo.h \
			PlaneVBOX.h \
			PlaneMwpc.h PlaneScifiG.h PlaneBMS.h \
			PlaneRICH.h PlaneRICH_MAPMT.h PlaneAPV.h \
                        PlaneGEM.h PlanePGEM.h PlaneSili.h \
			PlaneRiAPV.h PlaneRiAPVPanel.h \
			PlaneHCAL1.h PlaneHCAL2.h PlaneHCALT.h \
                        PlaneECAL0.h \
			PlaneHCAL2Sum.h \
                        PlaneECAL1.h PlaneECAL1Panel.h \
                        PlaneFEM.h PlaneFEMPanel.h \
                        PlaneECAL2.h PlaneECAL2Panel.h \
                        PlaneMuonWallA.h PlaneMuonWallB.h \
			PlaneScaler.h\
			Plane2V.h PlaneRCath.h \
			PlaneRPD_SADC.h PlaneRPD_F1.h \
			PlaneRPD_SADC.h \
			PlaneSandwich.h \
			PlaneSandwichPanel.h \
			PlaneTrigger_SADC_F1.h \
			PlaneStrawTubes.h \
			PlaneCamera.h GroupCamera.h \
			PlaneGandalf.h \
			MainFrame.h StartFrame.h \
			VariousSettings.h Geometry.h PlaneRPD_F1_Panel.h \
			PlaneDriftChamberPanel.h PlaneRPD_SADC_Panel.h\
			PlanePanel.h Plane1VPanel.h Plane2VPanel.h \
			Monitor.h Variable.h\
			MonitorPanel.h \
			Group.h GroupScifiJ.h GroupTrigHodo.h \
			GroupPanel.h GroupScifiJPanel.h GroupMumegaPanel.h \
			GroupMuonWallAPanel.h \
			GroupDAQ.h GroupMwpc.h GroupScifiG.h GroupMumega.h\
			GroupGEM.h GroupSili.h GroupMuonWallA.h GroupStraws.h \
			GroupPGEM.h \
			GroupBeamKiller.h \
			GroupCEDAR.h \
			GroupRiAPV.h \
			GroupRICH_MAPMT.h \
			GroupRPD.h \
			PlaneRW.h \
			PlaneTcsPhase.h \
			GroupRW.h GroupECAL1.h \
			BitButton.h Reference.h

MYLIBS =	-L../expat/xmlparse -lexpat
CXXFLAGS += -g -O3 -Wuninitialized 

#ifeq "$(LIB_TYPE)" "shared"
#  all:          libDaqDataDecoding.$(ver).so
#  libDaqDataDecoding.$(ver).so: $(OS)  $(OBJS)
#		$(CXX) -shared -o $@ $(OBJS)
#		ln -fs libDaqDataDecoding.$(ver).so libDaqDataDecoding.so
#else
#  all:          libDaqDataDecoding.$(ver).a
#  libDaqDataDecoding.$(ver).a: $(OS)  $(OBJS)
#		ar cr $@ $(OBJS)
#		ln -fs libDaqDataDecoding.$(ver).a libDaqDataDecoding.a
#endif

all :	$(OS)	coool bcoool

show:	
	@echo $(CXXFLAGS)

$(OS):
	@echo "creating platform directory for system $(OS)"
	@mkdir $(OS)

coool:	 $(OS) $(OBJS) $(OS)/coool.o
		$(CXX) -v
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(OS)/coool.o $(LIBS) $(MYLIBS)

cooolGJ:	$(OS) $(OBJS) $(OS)/coool.o
		$(CXX) -v
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(OS)/coool.o $(LIBS) $(MYLIBS)

bcoool:	 $(OS) $(OBJS)	$(OS)/bcoool.o
		$(CXX) -v
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(OS)/bcoool.o $(LIBS) $(MYLIBS)

Debug:		 $(OS)/Debug.o
		$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(MYLIBS)


Dict.h :  Dict.cc
Dict.cc :	$(dict_headers) LinkDef.h
			$(ROOTSYS)/bin/rootcint -f Dict.cc -c                 \
			$(CXXFLAGS) $(dict_headers) LinkDef.h

$(OS)/VariousSettings.o: VariousSettings.cc
		$(CXX) -I. -I$(top_dir)/expat/xmlparse $(CXXFLAGS) -c VariousSettings.cc -o $@

.PHONY: clean
clean:
		@rm -rf $(OS) *.o *.so *.a *.root *.ps *.log log doc coool bcoool \
		Dict.*

.PHONY: distclean
distclean: clean
		@rm -fr config.h .depend

.depend:
		-$(CXX) $(CXXFLAGS) -I. -MM *.cc | sed -e 's|^\(.*\.o:\)|$(OS)/\1|' > .depend

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include        .depend
endif
endif

