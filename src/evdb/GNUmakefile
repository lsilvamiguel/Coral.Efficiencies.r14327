top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsObjStore

OPTIM_LEVEL = debug

APPL_SRCS = CsDateStore.$(SRC_EXT) CsDDDStore.$(SRC_EXT)

include $(top_srcdir)/makefiles/General.mk


C_INCLUDES += -I$(top_srcdir)/src/evdis

production: prod util

prod: $(OS)/eventout

util: $(OS)/dateReader $(OS)/eventout $(OS)/sizeout

$(OS)/dateReader: dateReader.cc
	g++ dateReader.cc CsDateReader.cc -O4 -o $(OS)/dateReader $(C_INCLUDES) \
	-L/$(CORAL)/lib/$(OS) \
	$(CORAL)/src/DaqDataDecoding/expat/xmlparse/libexpat.a \
	$(CORAL)/src/DaqDataDecoding/src/libDaqDataDecoding.a \
	$(CORAL)/src/DaqDataDecoding/expat/xmlparse/libexpat.a \
	-pthread -L/usr/local/lib -lshift -lnsl -ldl

$(OS)/sizeout: sizeout.cc
	g++ sizeout.cc CsDateReader.cc -O4 -o $(OS)/sizeout $(C_INCLUDES) \
	-L/$(CORAL)/lib/$(OS) \
	$(CORAL)/src/DaqDataDecoding/expat/xmlparse/libexpat.a \
	$(CORAL)/src/DaqDataDecoding/src/libDaqDataDecoding.a \
	$(CORAL)/src/DaqDataDecoding/expat/xmlparse/libexpat.a \
	-pthread -L/usr/local/lib -lshift -lnsl -ldl

$(OS)/eventout: eventout.cc
	g++ eventout.cc -o $(OS)/eventout -L/usr/local/lib -lshift -lnsl \
	$(C_INCLUDES) -I$(CORAL)/include/DaqDataDecoding \
	-L/$(CORAL)/lib/$(OS) -lDaqDataDecoding -lexpat \
	-L/usr/local/lib -lshift -lnsl
