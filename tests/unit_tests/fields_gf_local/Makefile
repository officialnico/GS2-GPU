
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_fields_gf_local
clean_individual_test: LOCAL_TEST_FILES:=test_fields_gf_local
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
ifeq ($(shell expr ${NTESTPROCS} \< 6; echo $?),1)
all:
	$(call SKIPCOMMAND,current NTESTPROCS=${NTESTPROCS}; set NTESTPROCS>=6 to run this test)
else
all: build_test_exe
	$(call RUNCOMMAND,test_fields_gf_local,test_fields_gf_local.in)
endif
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

