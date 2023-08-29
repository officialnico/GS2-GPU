
sinclude ../Makefile.tests_common

TEST_BASE_NAME:=slab_itg_s_hat_0_low_res
clean_individual_test: LOCAL_TEST_FILES:=${TEST_BASE_NAME} gs2_slab_itg
build_test_exe: ${TEST_BASE_NAME}
COPY:=cp ../../../bin/gs2 gs2_slab_itg

ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: build_test_exe
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_hybrid.in)
	$(call COPY)
	$(call RUNCOMMAND,gs2_slab_itg,${TEST_BASE_NAME}.in)
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

