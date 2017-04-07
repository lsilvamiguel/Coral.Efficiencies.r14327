top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY = CsCalorim

APPL_SRCS =                 \
CsCalorimeter.cc            \
CsCalorimeterDecode.cc      \
CsDigitizerSADC.cc          \
DigitizerSADCBase.cc          \
DigitizerSADCN.cc          \
DigitizerSADCF.cc          \
CsECAL0.cc                  \
CsECAL1.cc                  \
CsECAL2.cc                  \
CsHCAL1.cc                  \
CsHCAL2.cc                  \
TrigGroupTest.cc            \
CsCalorimeterGUI.cc         \
CsTimingInfo.cc

ifeq "$(USE_Qt)" "yes"

APPL_SRCS   += CsCalorimeterGUI_moc.cc
CLEAN_FILES += CsCalorimeterGUI_moc.cc

endif

# too many warnings in ROOT headers
# C_FLAGS += -Wundef

# too many warnings in DDD and CLHEP headers
# C_FLAGS += -Wshadow

include $(top_srcdir)/makefiles/General.mk


ifeq "$(USE_Qt)" "yes"

$(DEPEND_FILE): CsCalorimeterGUI_moc.cc

CsCalorimeterGUI_moc.cc: CsCalorimeterGUI.h
	$(Qt_MOC) -o CsCalorimeterGUI_moc.cc CsCalorimeterGUI.h

endif
