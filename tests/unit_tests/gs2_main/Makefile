
sinclude ../Makefile.tests_common

clean_local_test:
	-$(QUIETSYM)rm -f test_gs2_main
clean_individual_test: clean_local_test

ifdef TESTNORUN
all: test_gs2_main 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
all: test_gs2_main 
	$(call RUNCOMMAND,test_gs2_main,test_gs2_main.in)
endif

test_gs2_main.o:  ${TEST_MODS} 

test_gs2_main: test_gs2_main.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

