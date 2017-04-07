# Using conditional include here because "make distclean" should work, even if
# Makefile.inc doesn't exist (yet).
-include        Makefile.inc

.PHONY: default
default:
		@$(MAKE) -C src

.PHONY: all
all: default examples

.PHONY: examples
examples: default
		@$(MAKE) -C examples/how-to
ifeq "$(USE_DATE_LIB)" "yes"
		@$(MAKE) -C examples/MurphyTV
endif

.PHONY: doc
doc:
		doxygen

#		-umask 022;$(MKDIRHIER) $(stylesheets)
#		set -e; for i in dsssl/*.dsl dsssl/*.cat; do \
#		$(INSTALL_DATA) $$i $(stylesheets); \
#		done;

.PHONY: install
install:        default doc examples
		$(MKDIRHIER) $(includedir)
		set -e; for i in src/*.h; do \
		$(INSTALL_DATA) $$i $(includedir);\
		done;
		$(MKDIRHIER) $(libdir)
		-$(INSTALL_DATA) src/libDaqDataDecoding.a  $(libdir)/libDaqDataDecoding.a
		-$(INSTALL_DATA) src/libDaqDataDecoding.so $(libdir)/libDaqDataDecoding.so
		$(MKDIRHIER) $(bindir)
		@echo $(MKDIRHIER) $(docdir)/DaqDataDecoding
		$(MKDIRHIER) $(docdir)
		-set -e; for i in doc/html/*; do \
		$(INSTALL_DATA) $$i $(docdir);\
		done;
		$(MKDIRHIER) $(sysconfdir)
		$(INSTALL_DATA) Makefile.inc $(sysconfdir)
		$(INSTALL_PROGRAM) examples/how-to/maps $(bindir)
		$(INSTALL_PROGRAM) examples/how-to/ddd $(bindir)

#		$(INSTALL_DATA) $$i $(docdir)/DaqDataDecoding/$$i; \
#ok		echo $(subst doc/html/,,doc/html/aa); \

new_Doxyfile    := perl -e '{@a=`cat Doxyfile` ; for (@a) { /PROJECT_NUMBER.*(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'
#new_version_cc  := perl -e '{@a=`cat src/version.cc`; $$new = "$(ver)"; $$new=~s/\./,/g; for (@a) { /DaqDataDecoding_version.*(\d+,\d+,\d+)/ && (s/$$1/$$new/)  ; print $$_}}'
new_main_doc    := perl -e '{@a=`cat main.doc` ; for (@a) { /DaqDataDecoding\-(\d+\.\d+\.\d+)/ && (s/$$1/$(ver)/)  ; print $$_}}'

.PHONY: release
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

.PHONY: az_www
az_www:
		rm -rf /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding
		mv DaqDataDecoding-$(ver)/doc/html /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding
		mv DaqDataDecoding-$(ver).tgz /afs/cern.ch/user/z/zvyagin/www/DaqDataDecoding

.PHONY: check
check:
		@$(MAKE) -C examples/how-to check

.PHONY: clean
clean:
		@$(MAKE) -C src               clean
		@$(MAKE) -C examples/how-to   clean
ifeq "$(USE_DATE_LIB)" "yes"
		@$(MAKE) -C examples/MurphyTV clean
endif

.PHONY: distclean
distclean:
		@$(MAKE) -C src               distclean
		@$(MAKE) -C examples/how-to   distclean
ifeq "$(USE_DATE_LIB)" "yes"
		@$(MAKE) -C examples/MurphyTV distclean
endif
		@rm -fr Makefile.inc config.h config.status config.cache config.log config.log doc *.tgz DaqDataDecoding-*

.PHONY: depend
depend:
		@$(MAKE) -C src               depend
		@$(MAKE) -C examples/how-to   depend
ifeq "$(USE_DATE_LIB)" "yes"
		@$(MAKE) -C examples/MurphyTV depend
endif

.PHONY: errors_list
errors_list:
		@find -regex ".*\.cc\|.*\.h" | xargs sed -n --expression 's/\(.*DaqError("\)\(.*\)\(".*\)/\2/p' | sort | uniq
