
sinclude ../Makefile.tests_common

TEST_BASE_NAME:=test_gs2_diagnostics_new
clean_individual_test: LOCAL_TEST_FILES:=${TEST_BASE_NAME} old_diagnostics/.test_* old_diagnostics/*.{o,error,out*,vres,*mom*,eigenfunc,fields,g,heat*,jext,lpc,parity,phase,exit_reason} averages.dat tmpdata.dat
build_test_exe: ${TEST_BASE_NAME}

ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifdef USE_NEW_DIAG
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else

# Here we are essentially running two tests. First we run the old
# diagnostics. Then we run the new diagnostics and check the results
# are the same. The we run the new diagnostics for half the number of
# steps. Then we restart and run for another half, copying the file
# from the first half and appending to it, and check that the result
# is still unchanged.
all: build_test_exe
	$(call RUNCOMMAND,${TEST_BASE_NAME},old_diagnostics/${TEST_BASE_NAME}.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}.in)
	rm -f ${TEST_BASE_NAME}_start.out.nc
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_start.in)
	cp ${TEST_BASE_NAME}_start.out.nc ${TEST_BASE_NAME}_append.out.nc; echo
	$(call RUNCOMMAND,${TEST_BASE_NAME},${TEST_BASE_NAME}_append.in)

endif
else
all:
	$(call SKIPCOMMAND,compiled without new diagnostics; compile with USE_NEW_DIAG=on to run this test)
endif
endif

${TEST_BASE_NAME}.o:  ${TEST_MODS} 

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

