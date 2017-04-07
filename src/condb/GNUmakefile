
all : bin lib

bin:
	mkdir -p bin

depend:
	$(MAKE) -C handler -f GNUmakefile depend

lib : Handler CDB_tools

Handler:
ifndef VERBOSE	
	@$(MAKE) -C handler -f GNUmakefile
else
	$(MAKE) -C handler -f GNUmakefile
endif

CDB_tools:
ifndef VERBOSE	
	@$(MAKE) -C cdb_tools -f Makefile	
else
	$(MAKE) -C cdb_tools -f Makefile	
endif

install:
ifndef VERBOSE	
	@$(MAKE) -C handler -f GNUmakefile install
else
	$(MAKE) -C handler -f GNUmakefile install
endif

clean:
	$(MAKE) -C handler -f GNUmakefile clean
#	$(MAKE) -C cdb_tools -f GNUmakefile clean
