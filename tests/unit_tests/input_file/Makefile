
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=gs2
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} fort.52
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
	$(call RUNCOMMAND,${TEST_BASE_NAME},empty.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},invalid_input_name.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},invalid_input_name_indexed.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},invalid_via_include_pre.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},invalid_via_include_post.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},doesnt_exist.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},circular_dependency.in)
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},nested_include.in)
# Intentionally no input file specified
	$(call RUNCOMMAND_EXPECT_FAIL,${TEST_BASE_NAME},)
endif
endif

${TEST_BASE_NAME}: ${BINDIR}/gs2
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)ln -sf $^ $@
