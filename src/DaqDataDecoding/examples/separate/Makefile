# You must change the line with
#   DaqDataDecoding = <path>
# The name 'DaqDataDecoding' should be path to installed
# version of DaqDataDecoding library. It should NOT point to DaqDataDecoding
# source code!
DaqDataDecoding = 

include         $(DaqDataDecoding)/etc/Makefile.inc
CXXFLAGS        += -I$(DaqDataDecoding)/include

simple:		simple.o
		$(CXX) -o simple simple.o $(DaqDataDecoding)/lib/libDaqDataDecoding.a $(LIBS)

clean:
		@rm -rf *.o simple
