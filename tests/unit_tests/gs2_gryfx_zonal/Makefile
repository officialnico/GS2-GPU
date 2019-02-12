
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_gryfx_zonal

all: test_gs2_gryfx_zonal 
	$(call RUNCOMMAND,test_gs2_gryfx_zonal,test_gs2_gryfx_zonal.in)

test_gs2_gryfx_zonal.o:  ${TEST_MODS} 

test_gs2_gryfx_zonal: test_gs2_gryfx_zonal.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

