sinclude ../Makefile.tests_common
TEST_BASE_NAME:=cyclone_itg_collisions
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
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_low.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_low_ae.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_split_low.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_split_low_ae.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_high.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_high_ae.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_split_high.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_split_high_ae.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_high_nonradau_old.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_regular_high_nonradau_new.in)
endif
endif
${TEST_BASE_NAME}.o:  ${TEST_MODS}

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

update_golden_answers: ${TEST_BASE_NAME}
	$(call messageOp,"Updating golden answers for","\\t:",${TEST_BASE_NAME},"","")
	$(QUIETSYM)$(GK_HEAD_DIR)/scripts/update_golden_answers.sh ${TEST_BASE_NAME}
