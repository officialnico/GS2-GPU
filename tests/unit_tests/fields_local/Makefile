
sinclude ../Makefile.tests_common
TEST_BASE_NAME=test_fields_local
clean_individual_test: LOCAL_TEST_FILES:=${TEST_BASE_NAME} ${TEST_BASE_NAME}_tuning

build_test_exe: ${TEST_BASE_NAME} ${TEST_BASE_NAME}_tuning
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
	$(call RUNCOMMAND,${TEST_BASE_NAME}_tuning,${TEST_BASE_NAME}_tuning.in)
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

${TEST_BASE_NAME}_tuning.o:  ${TEST_MODS}

${TEST_BASE_NAME}_tuning: ${TEST_BASE_NAME}_tuning.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

