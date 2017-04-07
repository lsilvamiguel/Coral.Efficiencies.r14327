top_srcdir = ../..
LIB_TYPE = static
include $(top_srcdir)/makefiles/Include.mk

##### Program name
MYPROGRAM = myprogram

##### Source files
APPL_SRCS = main.cc CoralUser.cc

##### Options files
OPTIONS_FILES = trafdic.2002.opt \
                trafdic.2003.opt \
                trafdic.2004.opt \
                trafdic.2006.opt \
                trafdic.2007.opt trafdic.2007.transv.opt \
                trafdic.2008.minus.opt trafdic.2008.plus.opt \
                trafdic.2009.opt trafdic.2009.dvcs.opt trafdic.2009.primakoff.opt CsPPI.2009.opt CsPPI.2009.dvcs.opt \
                trafdic.2010.opt \
                trafdic.2011.opt \
                trafdic.h2012.opt trafdic.m2012.opt trafdic.2012.DY.opt \
		trafdic.2014.opt trafdic.2015.opt

DISTCLEAN_FILES += $(OPTIONS_FILES)

##### what's to do...
all : $(OS)/ $(OS)/$(MYPROGRAM) options_files

include $(top_srcdir)/makefiles/General.mk

$(OS)/$(MYPROGRAM): $(APPL_OBJS) $(top_srcdir)/lib/$(OS)/*.a
	$(C++) $(C_FLAGS) -o $(OS)/$(MYPROGRAM) $(APPL_OBJS) $(LIBS) 

##### Options files: C2OPT  (i.e. (.opt.)c file to .opt file)
# Some of the options files are to be compiled by cpp from a template ".opt.c"
# counterpart to yield a file valid on local computer system.
options_files: $(OPTIONS_FILES)

# if $(HOST) is empty set it to the current hostname
ifeq ("$(HOST)", "")
HOST := $(shell hostname)
endif

# Define some var's used in the cpp compilation, depending upon local system
C2OPT_HEADER:="s%C2OPT_HEADER%// .opt file created from .opt.c counterpart on $(shell date)%"
# Which computer system? => set "C2OPT_COMPUTERFARM"
ifneq (,$(findstring pcco,$(HOST)))     # On-line
C2OPT_COMPUTERFARM:=4
endif
ifneq (,$(findstring ccage,$(HOST)))    # @ Lyon 
C2OPT_COMPUTERFARM:=3
endif
ifneq (,$(findstring compass,$(HOST)))  # @ gridKa
C2OPT_COMPUTERFARM:=2
endif
ifneq (,$(findstring lxplus,$(HOST)))   # @ CERN
C2OPT_COMPUTERFARM:=1
endif
ifndef C2OPT_COMPUTERFARM
C2OPT_COMPUTERFARM:=1
endif
C2OPT_FLAG:=0     # To be overriden on command line to enable special settings

# C2OPT is based on cpp (invoked by gcc -E)
# Since cpp action fluctuates from one gcc version to the next, it is
# supplemented by some sed commands (A very bad situation, I agree... Although
# these fuctuations affect only comments.)
# Use sed to:
# - Prepend header
# - Remove #1 directives
# - Change other #line directives to blank lines (This allows to clear up
#  the text)
# - Insert TAB before "//" introducing ``in-line'' comment
%.opt: %.opt.c GNUmakefile
	@echo $(CC) $< ">" $@
	@$(CC) -E -ffreestanding -traditional-cpp -C \
		-DC2OPT_COMPUTERFARM=$(C2OPT_COMPUTERFARM) \
		-DC2OPT_FLAG=$(C2OPT_FLAG) \
		$< | \
	sed -e "s%C2OPT_HEADER%// \"$@\" file created from \"$<\", on $(shell date)%" \
	    -e "/^# 1/ d" \
	    -e "/^#/ s/.*//" \
	    -e "/^\/\//! s%\([^/ ] \)//%\1	//%" \
	> $@

%.transv.opt: %.opt.c GNUmakefile
	@echo $(CC) $< ">" $@
	@$(CC) -E -ffreestanding -traditional-cpp -C \
		-DC2OPT_COMPUTERFARM=$(C2OPT_COMPUTERFARM) \
		-DC2OPT_FLAG=$(C2OPT_FLAG) -DC2OPT_TRANSVERSE=1 \
		$< | \
	sed -e "s%C2OPT_HEADER%// \"$@\" file created from \"$<\", on $(shell date)%" \
	    -e "/^# 1/ d" \
	    -e "/^#/ s/.*//" \
	    -e "/^\/\//! s%\([^/ ] \)//%\1	//%" \
	> $@

%.longitud.opt: %.opt.c GNUmakefile
	@echo $(CC) $< ">" $@
	@$(CC) -E -ffreestanding -traditional-cpp -C \
		-DC2OPT_COMPUTERFARM=$(C2OPT_COMPUTERFARM) \
		-DC2OPT_FLAG=$(C2OPT_FLAG) -DC2OPT_LONGITUDINAL=1 \
		$< | \
	sed -e "s%C2OPT_HEADER%// \"$@\" file created from \"$<\", on $(shell date)%" \
	    -e "/^# 1/ d" \
	    -e "/^#/ s/.*//" \
	    -e "/^\/\//! s%\([^/ ] \)//%\1	//%" \
	> $@
