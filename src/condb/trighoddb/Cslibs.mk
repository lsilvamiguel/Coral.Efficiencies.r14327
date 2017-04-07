# $Id: Cslibs.mk,v 1.3 2009/04/08 07:59:37 suhl Exp $
CS_LIB = -L$(top_srcdir)/lib/$(OS) \
		-lCsCond	\
		-lCsTransMgr	\
		-lCsBase     \
		-lCsHist     \
		-lCsEvdis     \
		-lCsDecod    \
		-lCsEvent    \
		-lCsTrack    \
		-lCsEvmc     \
		-lCsObjStore  \
		-lCsGeom     \
		-lCsGEM      \
		-lCsRaw      \
		-lCsUtils    
#	 -lTraffic    

#C_INCLUDES += -I$(top_srcdir)/include/ 
# Again a trouble fro ROOT. This is a temporary solution...
C_INCLUDES += -I$(top_srcdir)/include/ -I$(top_srcdir)/src/evdis/ -I$(top_srcdir)/src/hist/
