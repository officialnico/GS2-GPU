
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_nonlinear_terms
clean_individual_test: LOCAL_TEST_FILES:=${TEST_BASE_NAME}
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
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
ifndef USE_FFT
all:
	$(call SKIPCOMMAND,compiled without FFT library; compile with USE_FFT=on to run this test)
else
all: build_test_exe
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}.in)
endif
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

