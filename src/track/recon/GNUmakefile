top_srcdir = ../../..
include $(top_srcdir)/makefiles/Include.mk

#
# Makefile to build COMPASS 
# track reconstruction program  in the first spectrometer (RECON)
# library
# 	

LIBRARY = Recon

APPL_SRCS = sources/CsReconBridging.$(SRC_EXT) \
            sources/CsReconFitting.$(SRC_EXT) \
            sources/CsReconPrepattern.$(SRC_EXT) \
            sources/GetClusters.$(SRC_EXT) \
            sources/RecCluster.$(SRC_EXT) \
            sources/RecIni.$(SRC_EXT) \
            sources/Recon.$(SRC_EXT) \
            sources/RecOpt.$(SRC_EXT)
F77_SRCS = sources/Fortran/assoc.$(F77_EXT) \
           sources/Fortran/check_hit.$(F77_EXT) \
           sources/Fortran/clean_up.$(F77_EXT) \
           sources/Fortran/extrap_sm2.$(F77_EXT) \
           sources/Fortran/find_best.$(F77_EXT) \
           sources/Fortran/follow.$(F77_EXT) \
           sources/Fortran/gluehits.$(F77_EXT) \
           sources/Fortran/mag_field.$(F77_EXT) \
           sources/Fortran/mc_ini.$(F77_EXT) \
           sources/Fortran/multi_coulomb.$(F77_EXT) \
           sources/Fortran/ntuple.$(F77_EXT) \
           sources/Fortran/planes1.$(F77_EXT) \
           sources/Fortran/planes2.$(F77_EXT) \
           sources/Fortran/planes3.$(F77_EXT) \
           sources/Fortran/proj.$(F77_EXT) \
           sources/Fortran/proj_fringe.$(F77_EXT) \
           sources/Fortran/rec_det1_mgmg.$(F77_EXT) \
           sources/Fortran/rec_det1_new.$(F77_EXT) \
           sources/Fortran/rec_det1_old.$(F77_EXT) \
           sources/Fortran/rec_det2.$(F77_EXT) \
           sources/Fortran/recon_end.$(F77_EXT) \
           sources/Fortran/recon.$(F77_EXT) \
           sources/Fortran/recon_ini.$(F77_EXT) \
           sources/Fortran/recon_out.$(F77_EXT) \
           sources/Fortran/remdup.$(F77_EXT) \
           sources/Fortran/remover_hits.$(F77_EXT) \
           sources/Fortran/setup.$(F77_EXT) \
           sources/Fortran/tkfit_frext.$(F77_EXT) \
           sources/Fortran/tkfit_frfl.$(F77_EXT) \
           sources/Fortran/tkfit_pr.$(F77_EXT) \
           sources/Fortran/util.$(F77_EXT)

C_INCLUDES += -Iincludes/ -Iincludes/INCL/

FFLAGS += -Iincludes/INCL/

ifndef OS
  OS=Linux
  ERROR:=$(shell echo WARNING: environment variable OS is not defined\! [ will be = Linux] >& 2 )
endif

include $(top_srcdir)/makefiles/General.mk

exe :	lib install
	rm -f $(top_srcdir)/src/hardcoral/$(OS)/hardcoral
	$(MAKE) -C $(top_srcdir)/src/hardcoral/
	ln -f -s $(top_srcdir)/src/hardcoral/$(OS)/hardcoral $(EXE)
	@echo
	@echo $(EXE) "Done"

Doc: 
	@echo
	@echo "Generate documentation ..."
	@cd doxygen; rm -rf html; rm -rf latex; rm -rf man; doxygen; cd ..
	@echo

All:	clean Depend lib install exe Doc






