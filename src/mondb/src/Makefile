include         ../Makefile.inc

#MONITOR          =  /afs/cern.ch/compass/detector/monitor
#MYSQL            =  $(MONITOR)/mysql

#CXXFLAGS        += -Wall -I$(top_dir)/src -I$(top_dir)/expat/xmlparse
#CXXFLAGS        += -I$(MYSQL)/include
#CXXFLAGS        += -g -Wall

#LD               = g++
#SOFLAGS          = -shared 

# LIBS             = -L/usr/lib/mysql -lmysqlclient -lz 
#LIBS            += -lThread  -L$(MYSQL)/lib -lmysqlclient -lz

OBJS             =  Utilities.o monDB.o Add2DB.o SyncDB.o LsDB.o CatDB.o




all :  SyncDB Add2DB LsDB CatDB libmonDB.a Test TestSpeed

Add2DB :	Add2DB.o Utilities.o
		$(CXX)  -o $@ $^ $(LIBS)

SyncDB :	SyncDB.o Utilities.o
		$(CXX)  -o $@ $^ $(LIBS)

LsDB :		LsDB.o Utilities.o
		$(CXX)  -o $@ $^ $(LIBS)

CatDB :		CatDB.o Utilities.o
		$(CXX)  -o $@ $^ $(LIBS)

libmonDB.a:  	monDB.o
		ar cr $@ $^
		cp $@ ../lib/

Test :		Test.o libmonDB.a
		$(CXX)  -o $@ Test.o ../lib/libmonDB.a $(LIBS)

TestSpeed :	TestSpeed.o libmonDB.a
		$(CXX)  -o $@ TestSpeed.o ../lib/libmonDB.a $(LIBS)

clean:
		@rm -rf *.so *.a *.o .depend *.root *.ps *.log log doc *.o $(execs) SyncDB Add2DB LsDB CatDB
distclean:
		$(MAKE) clean
		@rm -fr .depend Makefile.inc config.h config.status config.cache doc 

depend:
		$(CXX) $(CXXFLAGS) -I. -MM *.cc > .depend

tags:           *.cc *.h
		ctags -C -d --globals --members -w *.cc

-include        .depend



