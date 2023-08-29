
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_gs2_init
clean_individual_test: LOCAL_TEST_FILES:=${TEST_BASE_NAME}
build_test_exe: ${TEST_BASE_NAME}

ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
ifndef USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test)
else
all: build_test_exe
	$(call RUNCOMMAND,test_gs2_init,test_gs2_init.in)
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

