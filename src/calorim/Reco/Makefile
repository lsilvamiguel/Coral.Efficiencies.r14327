-include Makefile.inc

SUBDIRS       = src test

.PHONY: default
default:		
		@$(MAKE) -C src


.PHONY: all
all: default test	


.PHONY: test
test: default
		@$(MAKE) -C test


.PHONY: clean
clean:
		@rm -rf $(clean_mask) doc
		@for i in $(SUBDIRS); do $(MAKE) -C $$i clean; done


.PHONY: distclean
distclean:
		@rm -rf $(clean_mask) doc
		@for i in $(SUBDIRS); do $(MAKE) -C $$i distclean; done
		@rm -rf config.cache config.log config.status Makefile.inc


.PHONY: depend
depend:
		@for i in $(SUBDIRS); do $(MAKE) -C $$i depend; done


.PHONY: doc
doc:
		doxygen

new_Doxyfile    := perl -e '{@a=`cat Doxyfile` ; for (@a) { /PROJECT_NUMBER.*(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'
#new_version_cc  := perl -e '{@a=`cat src/version.cc`; $$new = "$(ver)"; $$new=~s/\./,/g; for (@a) { /Reco_version.*(\d+,\d+,\d+)/ && (s/$$1/$$new/)  ; print $$_}}'
new_main_doc    := perl -e '{@a=`cat main.doc` ; for (@a) { /Reco\-(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'

.PHONY: release
release:
		rm -rf Reco-*
		$(new_Doxyfile)   > Doxyfile.new
		mv Doxyfile.new   Doxyfile
		$(new_main_doc)   > main.doc.new
		mv main.doc.new   main.doc
		cp -a $(PWD) Reco-$(ver)
		(cd Reco-$(ver)                                      ;\
		echo ``>Makefile.inc                                            ;\
		$(MAKE) distclean                                               ;\
		echo Version date is $(ver_date)                                ;\
		cd ..)
		tar czf Reco-$(ver).tgz Reco-$(ver)
