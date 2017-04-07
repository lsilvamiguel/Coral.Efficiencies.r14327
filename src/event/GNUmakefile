top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

C_FLAGS += -Wno-trigraphs

LIBRARY=CsEvent

APPL_SRCS = CsBuildParticles.cc  \
            CsCluster.cc         \
	    CsDigit.cc           \
	    CsEvent.cc           \
            CsEventUtils.cc      \
            CsMCDigit.cc         \
            CsMW1Pid.cc          \
	    CsParticle.cc        \
            CsRecoEvent.cc       \
	    PIDdoBeamID.cc       \
	    PIDdoMuonID.cc


# Tell to general which extra libs.... and remove them from general...

include $(top_srcdir)/makefiles/General.mk

C_INCLUDES += -I$(top_srcdir)/src/evdb -I$(DIR_ROOT)/include

ifeq "$(USE_ORACLE)" "yes"
  C_INCLUDES += -I$(top_srcdir)/src/oracle 
endif
