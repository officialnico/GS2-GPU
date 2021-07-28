
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_timestep_change
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} *.vspace_*

ifdef TESTNORUN
all: ${TEST_BASE_NAME}
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: ${TEST_BASE_NAME}
ifndef USE_NETCDF
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test on ${TEST_BASE_NAME}.in)
else
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}.in)
# Expect the following to fail due to changes in the fields
#	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_maxwell_reinit.in)
# Expect the following to fail due to losing the explicit source terms
#	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_no_explicit_in_restart.in)
endif
# Expect the following to fail until we merge PR #490 or similar
#	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_in_memory.in)
endif
endif
${TEST_BASE_NAME}.o:  ${TEST_MODS}

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
