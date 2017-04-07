top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

APPL_SRCS = CsPPI_EC02tools.cc

CLEAN_FILES += merger dumper

include $(top_srcdir)/makefiles/General.mk

all:  merger dumper

merger: $(OS)/merger.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ -o $@

dumper: $(OS)/dumper.o $(APPL_OBJS)
	@echo "($(C++)) creating $@"
	@$(C++) $(C_FLAGS) $^ -o $@

$(OS)/%.$(OBJ_EXT): %.$(SRC_EXT) $(SCHEMA_HDRS) | $(ODIRS)
	@echo "(C++) compiling $*"
ifndef VERBOSE
	@$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) $(C_FLAGS) $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
else
	$(C++) -c $*.$(SRC_EXT) $(C_INCLUDES) $(C_FLAGS) $(C_OUTPUT_TO)"$(OS)/$*.$(OBJ_EXT)"
endif


