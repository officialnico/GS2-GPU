
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_main

ifdef TESTNORUN
all:
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
all: test_gs2_main 
	$(call RUNCOMMAND,test_gs2_main,test_gs2_main.in)
endif

test_gs2_main.o:  ${TEST_MODS} 

test_gs2_main: test_gs2_main.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

