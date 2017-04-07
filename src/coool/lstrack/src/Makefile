include         ../Makefile.inc


OBJS  =  	Structures.o FitCTime.o \
		TrackMaker.o Tracker.o Station.o Doublet.o 

DICT_HDRS    := Structures.h \
                TrackMaker.h Tracker.h Station.h Doublet.h \

DICTM_HDRS   := $(DICT_HDRS) MMan_2triplets.h 


all :  liblstrack.a libMMan.so TestGeom

liblstrack.a:  $(OBJS) LsDict.o
	ar cr $@ $^
	cp $@ ../lib/

libMMan.so:  $(OBJS) MMan_2triplets.o MManGeometry.o Dictm.o
	$(CXX) -shared -o $@ $^
	cp $@ ../lib/

TestGeom: libMMan.so TestGeom.o 
	$(CXX) -o $@ $^ $(LIBS)
	cp $@ ../bin/

LsDict.cc LsDict.h : $(DICT_HDRS)
	$(ROOTSYS)/bin/rootcint -f LsDict.cc -c $(CXXFLAGS) $(DICT_HDRS) \
	LinkDef.h 

Dictm.cc Dictm.h : $(DICTM_HDRS)
	$(ROOTSYS)/bin/rootcint -f Dictm.cc -c $(CXXFLAGS) $(DICTM_HDRS) \
	LinkDefM.h

clean:
	@rm -rf *.so *.a *.o .depend *.root *.ps *.log log doc *.o \
	*Dict* $(execs) 
distclean:
		$(MAKE) clean
		@rm -fr .depend Makefile.inc config.h config.status config.cache doc 

depend:
		$(CXX) $(CXXFLAGS) -I. -MM *.cc > .depend

tags:           *.cc *.h
		ctags -C -d --globals --members -w *.cc

-include        .depend



