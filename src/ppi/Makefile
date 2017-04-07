top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY = CsPPI

APPL_SRCS = src/CsPP.cc \
            src/CsPPI.cc \
            src/CsPPI_EC02time.cc

C_INCLUDES += -Iinclude/

include $(top_srcdir)/makefiles/General.mk

all: $(MYLIB) tools

.PHONY: tools
tools: 
	@$(MAKE) -C tools all

.PHONY: clean
clean: tools-clean

.PHONY: tools-clean
tools-clean:
	@$(MAKE) -C tools clean

