top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

# $Id: GNUmakefile,v 1.72 2010/05/03 10:43:18 ybedfer Exp $
#
# Makefile to build "lattice" track fitting library for CORAL
# It is to be called from CORAL makefiles  with targets "clean" "lib" and
# "install"
# 	
#   Other targets are for building "lattice" utilities:
#  - "makeDico"
#  - an executable, "trafdic" 
#


ifndef OS
  OS=Linux
  ERROR:=$(shell echo WARNING: environment variable OS is not defined\! [ will be = Linux] >& 2 )
endif

CDIR=src
C_INCLUDES += -I./$(CDIR) -I./$(CDIR)/Aux

########## DEFINES

C_FLAGS+=-DCURVATURE_TOLERANCE
C_FLAGS+=-DDRIFTS_SUBDUED
#C_FLAGS+=-DQUADRATIC_PROJ
C_FLAGS+=-DCUT_KF_UPON_0
C_FLAGS+=-DSCIFIS_ENHANCED=2
C_FLAGS+=-DFRINGE_FIELD_KF
#C_FLAGS+=-DSKIP_XWCOMB
C_FLAGS+=-DPR_LOW_P_FORWARD
C_FLAGS+=-DPR_2ZONE_ZPROJ
C_FLAGS+=-DInActive_CHECK_URANGE
C_FLAGS+=-DPR_USE_GEM_CORR
C_FLAGS+=-DPR_GP_WHICHSIDE_INFO
C_FLAGS+=-DPR_ST_WHICHSIDE_INFO
ifdef PR_HADRON_2004    # This is specific to the 2004 hadron pilot data
C_FLAGS+=-DPR_ENHANCE_CENTRAL_GEMs
endif

C_FLAGS+=-DQNewton_ENABLED
#C_FLAGS+=-DTraFDIc_HISTO_STAT    # Extra ``define'' to be enabled on option

C_FLAGS+=-DEXTENSIVE_DICO

########## DEFINES ENABLED from COMMAND LINE
ifdef DEBUG
	C_FLAGS+=-DDEBUG
endif
ifdef ASSOCIATE_MIRRORS
	C_FLAGS+=-DASSOCIATE_MIRRORS
endif

########## SOURCES
LIBRARY = Trafdic

APPL_SRCS := $(wildcard $(CDIR)/*.cc) 
# Fortran sources
F77_SRCS  := $(wildcard src/Fortran/*.F)


########## INCLUDE coral MAKEFILES
include $(top_srcdir)/makefiles/General.mk

##########
Doc: 
	@echo
	@echo "Generate documentation ..."
	@mkdir -p doxygen
	@cd doxygen; rm -rf html; rm -rf latex; rm -rf man; doxygen ../Doxyfile; cd ..
	@echo
