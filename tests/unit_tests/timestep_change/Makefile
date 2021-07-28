
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_timestep_change
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME}

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
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}.in)

endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS}

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
