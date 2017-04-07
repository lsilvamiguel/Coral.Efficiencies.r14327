
MYSQL            =  /afs/cern.ch/compass/detector/monitor/mysql

CXXFLAGS        += -Wall -I$(top_dir)/src -I$(top_dir)/expat/xmlparse
CXXFLAGS        += -I$(MYSQL)/include
CXXFLAGS        += -g -Wall

LD               = g++ -g
#SOFLAGS          = -shared 

LIBS             =  -L$(MYSQL)/lib -lmysqlclient -lz

OBJS             =  MyInterface.o 

all :  libMyInterface.a Test


libMyInterface.a: $(OBJS) 	
		ar cr $@ $^

Test :		Test.o libMyInterface.a
		$(CXX)  -o $@ $^ $(LIBS)

clean:
		@rm -rf *.so *.a *.o .depend  

distclean:
		$(MAKE) clean
		@rm -fr .depend 

depend:
		$(CXX) $(CXXFLAGS) -I. -MM *.cc > .depend

tags:           *.cc *.h
		ctags -C -d --globals --members -w *.cc

-include        .depend



