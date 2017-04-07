top_srcdir = ../../..
LIB_TYPE = static
include $(top_srcdir)/makefiles/Include.mk

##### Program name
MYPROGRAM = led_laser

##### what's to do...
all : $(OS)/ $(OS)/$(MYPROGRAM)

##### Source files
APPL_SRCS = main.cc CoralUser.cc RunDataBase.cc StoreLED.cc DevMonitorFEM.cc StoreDeltaADC.cc

include $(top_srcdir)/makefiles/General.mk

$(OS)/$(MYPROGRAM): $(APPL_OBJS) $(top_srcdir)/lib/$(OS)/*.a
	$(C++) $(C_FLAGS) -o $(OS)/$(MYPROGRAM) $(APPL_OBJS) $(LIBS) 

