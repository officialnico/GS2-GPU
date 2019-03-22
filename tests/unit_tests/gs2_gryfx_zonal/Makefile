
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_gryfx_zonal

ifdef TESTNORUN
all:
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
all: test_gs2_gryfx_zonal 
	$(call RUNCOMMAND,test_gs2_gryfx_zonal,test_gs2_gryfx_zonal.in)
endif

test_gs2_gryfx_zonal.o:  ${TEST_MODS} 

test_gs2_gryfx_zonal: test_gs2_gryfx_zonal.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

