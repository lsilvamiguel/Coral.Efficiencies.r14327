# $Id: General.mk,v 1.81 2010/11/03 16:45:46 suhl Exp $

ifndef top_srcdir
  MSG:=$(shell echo ERROR: top_srcdir is not defined\! >& 2 )
  exit
else
  include $(top_srcdir)/makefiles/Include.mk
endif

ifeq "$(ObjectFilesDeletion)" "yes"
  DELETEOBJS = rm $(OS)/*.o 	
else
  DELETEOBJS = echo -n
endif

CONFIG = $(top_srcdir)

C_OBJS    = $(addprefix $(OS)/, $(APPL_SRCS:.$(SRC_EXT)=.$(OBJ_EXT)))
F77_OBJS  = $(addprefix $(OS)/, $(F77_SRCS:.$(F77_EXT)=.$(OBJ_EXT)))
APPL_OBJS = $(F77_OBJS) $(C_OBJS)

ODIRS     = $(sort $(OS)/ $(dir $(APPL_OBJS)))

# This works only with one value for TEST_PROGRAMS
# In order to make it work with more files, the
# rules must be modified
TEST_MAIN = $(addsuffix .$(SRC_EXT), $(TEST_PROGRAMS))
TEST_OS   = $(addsuffix .$(OS), $(TEST_PROGRAMS))
TEST_OBJ  = $(addprefix $(OS)/, $(TEST_MAIN:.$(SRC_EXT)=.$(OBJ_EXT)))

OBJS = $(APPL_OBJS)

ifdef VERBOSE
MSG:=$(shell echo ----- >& 2 )
MSG:=$(shell echo C_FLAGS =  $(C_FLAGS) >& 2 )
MSG:=$(shell echo C_INCLUDES =  $(C_INCLUDES) >& 2 )
MSG:=$(shell echo LIBS =  $(LIBS) >& 2 )
MSG:=$(shell echo ----- >& 2 )
endif

.PHONY: all
all: lib install

.PHONY: default
default: all

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(DEPEND_FILE)
endif
endif


#
# C++ Compilation
#
$(C_OBJS) $(TEST_OBJ): $(OS)/%.$(OBJ_EXT): %.$(SRC_EXT) $(SCHEMA_HDRS) | $(ODIRS)
	@echo "(C++) compiling $*"
ifdef VERBOSE
	$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
else
	@$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) ${C_FLAGS} $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
endif


#
# Fortran Compilation
#
$(F77_OBJS): $(OS)/%.$(OBJ_EXT): %.$(F77_EXT) | $(ODIRS)
	@echo "(Fortran) compiling $*"
ifdef VERBOSE
	$(F77) -c $*.$(F77_EXT) $(FFLAGS) $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
else
	@$(F77) -c $*.$(F77_EXT) $(FFLAGS) $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
endif


#
# System specific C++ settings
#
ifdef LIBRARY
  ifeq "$(LIB_TYPE)" "shared"
    MYLIB   = $(OS)/lib$(LIBRARY).$(SHLIB_EXT)
#    MAKELIB = $(MKSHLIB) $(MYLIB) $(SHLIB_FLAGS) $(C_FLAGS)  $(APPL_OBJS)
    MAKELIB = $(MKSHLIB) $(MYLIB) $(C_FLAGS)  $(OBJS)
    C_FLAGS += -fPIC
  else
    MYLIB = $(OS)/lib$(LIBRARY).$(ARCH_EXT)
#    MAKELIB = ar -r  $(MYLIB) $(APPL_OBJS) 
    MAKELIB = ar -cr  $(MYLIB) $(OBJS) 
  endif
  INSTALLCPP = $(top_srcdir)/lib/$(MYLIB)
endif

$(ODIRS):
	@echo creating platform directory for system "$@"
	@mkdir -p $@

.PHONY: lib
lib: $(MYLIB) 

.PHONY: install
install: $(INSTALLCPP) 

$(INSTALLCPP): $(MYLIB)
	@mkdir -p $(top_srcdir)/lib/$(OS)
	@if ! test -L $(INSTALLCPP); then echo install $(MYLIB) in the development area; ln -fs `pwd`/$(MYLIB) $(INSTALLCPP); fi

$(MYLIB): $(OBJS) | $(OS)/
	@echo creating library $(MYLIB)
	@echo $(OS) $(OBJS)
	@$(MAKELIB)
	@$(DELETEOBJS)

.PHONY: test
test: $(TEST_OS)

$(TEST_OS): $(TEST_OBJ) GNUmakefile | $(OS)/
	@echo preparing the executable of $(TEST_PROGRAMS)
	$(C++) $(C_FLAGS) -o $(OS)/$(TEST_PROGRAMS) $(TEST_OBJ) $(LIBS)

ifndef depend_used
.PHONY: depend
depend:
	@rm -f $(DEPEND_FILE)
	@$(MAKE) $(DEPEND_FILE)

$(DEPEND_FILE):
	@error=0 ; \
	count=0 ; \
	for i in $(APPL_SRCS) $(TEST_MAIN) $(filter %.$(SRC_EXT), $(EXTRA_SRCS)) ; do \
		$(C++) $(C_INCLUDES) $(C_FLAGS) -MM -MT `echo $(OS)/$$i | sed -e 's|\(.*\.\)$(SRC_EXT)$$|\1$(OBJ_EXT)|g'` $$i >> tmp_dep || error=1 ; \
		if test $$error -eq 1 ; then \
			rm -f tmp_dep ; \
			exit 1 ; \
		fi ; \
		count=$$((count+1)) ; \
	done ; \
	for i in $(F77_SRCS) $(filter %.(F77_EXT), $(EXTRA_SRCS)) ; do \
		$(C++) $(FFLAGS) -MM -MT `echo $(OS)/$$i | sed -e 's|\(.*\.\)$(F77_EXT)$$|\1$(OBJ_EXT)|g'` $$i >> tmp_dep || error=1 ; \
		if test $$error -eq 1 ; then \
			rm -f tmp_dep ; \
			exit 1 ; \
		fi ; \
		count=$$((count+1)) ; \
	done ; \
	if test $$count -gt 0 ; then \
		mv tmp_dep $(DEPEND_FILE) ; \
	fi 
endif

.PHONY: clean
clean: 
	@rm -rf core *~ \#*\# vc50.pdb vc50.pch .inslog2 tca.map tca.log \
	$(OS) $(TEST_OS) $(CLEAN_FILES) $(DEPEND_FILE)

.PHONY: distclean
distclean: clean
	@rm -rf $(DISTCLEAN_FILES)

.PHONY: show
show:
	@echo "OS:            $(OS)"
	@echo
	@echo "LIBS:          $(LIBS)"
	@echo
	@echo "INCDIRS:       $(INCDIRS)"
	@echo
	@echo "C_FLAGS:       $(C_FLAGS)"
	@echo
	@echo "C_INCLUDES:    $(C_INCLUDES)"
	@echo
	@echo "LIBRARY:       $(LIBRARY)"
	@echo
	@echo "MYLIB:         $(MYLIB)"
	@echo
	@echo "MAKELIB:       $(MAKELIB)"
	@echo
	@echo "INSTALLCPP:    $(INSTALLCPP)"
	@echo
	@echo "APPL_SRCS:     $(APPL_SRCS)"
	@echo
	@echo "F77_SRCS:      $(F77_SRCS)"
	@echo
	@echo "ODIRS:         $(ODIRS)"
	@echo
	@echo "C_OBJS:        $(C_OBJS)"
	@echo
	@echo "APPL_OBJS:     $(APPL_OBJS)"
	@echo
	@echo "TEST_OBJ:      $(TEST_OBJ)"
	@echo
	@echo "TEST_PROGRAMS: $(TEST_PROGRAMS)"
	@echo
	@echo "DEPEND_FILE:   $(DEPEND_FILE)"
	@echo
