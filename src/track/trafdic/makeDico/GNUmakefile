top_srcdir = ../../../..
include $(top_srcdir)/makefiles/Include.mk

# $Id: GNUmakefile,v 1.72 2010/05/03 10:43:18 ybedfer Exp $
#
# Makefile to build the makeDico helper.
#

ifndef OS
  OS=Linux
  ERROR:=$(shell echo WARNING: environment variable OS is not defined\! [ will be = Linux] >& 2 )
endif

CDIR=src
C_INCLUDES += -I../src -I../src/Aux

# define extra sources to create a dependency file
EXTRA_SRCS = $(CDIR)/makeDico.cc \
             $(CDIR)/CsField.cc

########## DEFINES ENABLED from COMMAND LINE
ifdef DEBUG
	C_FLAGS+=-DDEBUG
endif

########## INCLUDE coral MAKEFILES
include $(top_srcdir)/makefiles/General.mk

.DEFAULT_GOAL:=$(OS)/makeDico

########## RULE TO MAKE makeDico
$(OS)/makeDico:	$(OS)/$(CDIR)/makeDico.o $(OS)/$(CDIR)/CsField.o
	@echo "Linking: $(OS)/makeDico"
	$(C++) -o $(OS)/makeDico $(C_FLAGS) $(OS)/$(CDIR)/makeDico.o $(OS)/$(CDIR)/CsField.o $(LIBS)

########## makeDico.o
$(OS)/$(CDIR)/makeDico.o: ./$(CDIR)/makeDico.cc | $(OS)/$(CDIR)
	@echo "LATTICE compilation: makeDico.cc"
	@$(C++) -c $(C_INCLUDES) $(C_FLAGS) -o $(OS)/$(CDIR)/makeDico.o $(CDIR)/makeDico.cc

########## CsField.o, i.e. special version of "./src/geom/CsField.cc" for the
########## generation of the dico.
$(OS)/$(CDIR)/CsField.o: ./$(CDIR)/CsField.cc | $(OS)/$(CDIR)
	@echo "LATTICE compilation: CsField.cc"
	@$(C++) -c $(C_INCLUDES) $(C_FLAGS) -o $(OS)/$(CDIR)/CsField.o $(CDIR)/CsField.cc

########## directory for object files
$(OS)/$(CDIR):
	@echo creating platform directory for system "$@"
	@mkdir -p $@
