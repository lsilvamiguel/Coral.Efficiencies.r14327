include         Makefile.inc

.PHONY:		all clean distclean

all:
		$(MAKE) -C expat
		$(MAKE) -C src

doc:
		doxygen

#		-umask 022;$(MKDIRHIER) $(stylesheets)
#		set -e; for i in dsssl/*.dsl dsssl/*.cat; do \
#		$(INSTALL_DATA) $$i $(stylesheets); \
#		done;

install:
		$(MKDIRHIER) $(includedir)
		set -e; for i in src/*.h; do \
		$(INSTALL_DATA) $$i $(includedir);\
		done;
		$(MKDIRHIER) $(libdir)
		$(INSTALL_DATA) src/libDaqDataDecoding.a $(libdir)/libDaqDataDecoding.a
		$(MKDIRHIER) $(bindir)
		$(INSTALL_PROGRAM) examples/rtree/gui $(bindir)
		@echo $(MKDIRHIER) $(docdir)/DaqDataDecoding
		$(MKDIRHIER) $(docdir)
		-set -e; for i in doc/html/*; do \
		$(INSTALL_DATA) $$i $(docdir);\
		done;
		$(MKDIRHIER) $(sysconfdir)
		$(INSTALL_DATA) Makefile.inc $(sysconfdir)

#		$(INSTALL_DATA) $$i $(docdir)/DaqDataDecoding/$$i; \
#ok		echo $(subst doc/html/,,doc/html/aa); \

new_Doxyfile    := perl -e '{@a=`cat Doxyfile` ; for (@a) { /PROJECT_NUMBER.*(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'
#new_version_cc  := perl -e '{@a=`cat src/version.cc`; $$new = "$(ver)"; $$new=~s/\./,/g; for (@a) { /DaqDataDecoding_version.*(\d+,\d+,\d+)/ && (s/$$1/$$new/)  ; print $$_}}'
new_main_doc    := perl -e '{@a=`cat main.doc` ; for (@a) { /DaqDataDecoding\-(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'

release:
		rm -rf DaqDataDecoding-*
		$(new_Doxyfile)   > Doxyfile.new
		mv Doxyfile.new   Doxyfile
		$(new_main_doc)   > main.doc.new
		mv main.doc.new   main.doc
		cp -a $(PWD) DaqDataDecoding-$(ver)
		(cd DaqDataDecoding-$(ver)                                      ;\
		echo ``>Makefile.inc                                            ;\
		$(MAKE) distclean                                               ;\
		echo Version date is $(ver_date)                                ;\
		cd ..)
		tar czf DaqDataDecoding-$(ver).tgz DaqDataDecoding-$(ver)

az_www:
		rm -rf /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding
		mv DaqDataDecoding-$(ver)/doc/html /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding
		mv DaqDataDecoding-$(ver).tgz /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding

check:
		@$(MAKE) -C examples/how-to check

clean:
		@if [ -e src ];             then $(MAKE) -C src             clean; fi
		@if [ -e expat ];           then $(MAKE) -C expat           clean; fi

distclean: clean
		@if [ -e src ];           then $(MAKE) -C src               distclean; fi
		@rm -fr .depend Makefile.inc config.h config.status config.cache config.log doc *.tgz DaqDataDecoding-*

depend:
		@$(MAKE) -C src             .depend

errors_list:
		@find -regex ".*\.cc\|.*\.h" | xargs sed -n --expression 's/\(.*DaqError("\)\(.*\)\(".*\)/\2/p' | sort | uniq
# DO NOT DELETE
