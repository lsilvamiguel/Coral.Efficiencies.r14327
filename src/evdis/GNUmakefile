top_srcdir = ../..
include $(top_srcdir)/makefiles/Include.mk

LIBRARY=CsEvdis
TEST_PROGRAMS=

APPL_SRCS = CsEvdis.$(SRC_EXT)          \
            CsEvdisEvent.$(SRC_EXT)     \
            CsEvdisDict.$(SRC_EXT)      \
            CsEvdisMagField.$(SRC_EXT)  \
            CsEvdisTrack.$(SRC_EXT)

DICT_SRCS=$(filter-out CsEvdisDict.$(SRC_EXT), $(APPL_SRCS))
DICT_HDRS=$(DICT_SRCS:.$(SRC_EXT)=.h)

include $(top_srcdir)/makefiles/General.mk

CsEvdisDict.$(SRC_EXT): $(DICT_HDRS) CsEvdisLinkDef.h
	@echo "Generating ROOT dictionary $@..."
	@rootcint -f $@ -c -I$(top_srcdir)/include $^

CLEAN_FILES += *Dict.* *Dict_rdict.pcm

# ROOT dictionaries have to be created before the depandany file is created
$(DEPEND_FILE): CsEvdisDict.$(SRC_EXT)
