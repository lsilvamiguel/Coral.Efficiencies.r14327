top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk
include $(top_srcdir)/makefiles/General.mk

# Program name
MYPROGRAM = myprogram

# Source files
APPL_SRCS = main.$(SRC_EXT) CoralUser.cc

# what's to do...
all : $(OS) $(OS)/$(MYPROGRAM) 

$(OS)/$(MYPROGRAM): $(APPL_OBJS) 
	$(C++) $(C_FLAGS) $(C_INCLUDES) -o $(OS)/$(MYPROGRAM) $(APPL_OBJS) $(LIBS) 

