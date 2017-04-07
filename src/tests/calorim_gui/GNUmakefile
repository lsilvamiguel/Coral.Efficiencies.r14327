top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

##### Program name
MYPROGRAM = calorim_gui

##### what's to do...
all : $(OS)/ $(OS)/$(MYPROGRAM)

##### Source files
APPL_SRCS = main.cc

include $(top_srcdir)/makefiles/General.mk

$(OS)/$(MYPROGRAM): $(APPL_OBJS) $(top_srcdir)/lib/$(OS)/*.a
	$(C++) $(C_FLAGS) -o $(OS)/$(MYPROGRAM) $(APPL_OBJS) $(LIBS) 
