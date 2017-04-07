top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

# Program name
MYPROGRAM = richone

# Source files
APPL_SRCS = main.$(SRC_EXT) CoralUser.cc

include $(top_srcdir)/makefiles/General.mk

# what's to do...
.PHONY: all
all: $(OS)/$(MYPROGRAM) | $(OS)/

$(OS)/$(MYPROGRAM): $(APPL_OBJS)
	$(C++) $(C_FLAGS) $(C_INCLUDES) -o $(OS)/$(MYPROGRAM) $(APPL_OBJS) $(LIBS)


