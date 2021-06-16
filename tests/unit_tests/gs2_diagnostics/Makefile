
sinclude ../Makefile.tests_common

clean_local_test:
	-$(QUIETSYM)rm -f test_gs2_diagnostics
clean_individual_test: clean_local_test

ifdef TESTNORUN
all: test_gs2_diagnostics 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: test_gs2_diagnostics 
	$(call RUNCOMMAND,test_gs2_diagnostics,test_gs2_diagnostics.in)
	$(call RUNCOMMAND,test_gs2_diagnostics,test_gs2_diagnostics_nspec_2.in)
endif
endif

test_gs2_diagnostics.o:  ${TEST_MODS} 

test_gs2_diagnostics: test_gs2_diagnostics.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

