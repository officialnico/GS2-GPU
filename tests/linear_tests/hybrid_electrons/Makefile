
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=hybrid_electrons
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} *~ *.vspace_*
build_test_exe: ${TEST_BASE_NAME}
ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: build_test_exe
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_he_low_collisions.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_ke_low_collisions.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_ae_low_collisions.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_he_high_collisions.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_ke_high_collisions.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_ae_high_collisions.in)
endif
endif
${TEST_BASE_NAME}.o:  ${TEST_MODS}

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
