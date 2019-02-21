
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_init

ifndef USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test)
else
all: test_gs2_init 
	$(call RUNCOMMAND,test_gs2_init,test_gs2_init.in)
endif

test_gs2_init.o:  ${TEST_MODS} 

test_gs2_init: test_gs2_init.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

