
sinclude ../Makefile.tests_common

clean_individual_test: LOCAL_TEST_FILES:=test_indexing

ifdef TESTNORUN
all: test_indexing
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_indexing
	$(call RUNCOMMAND,test_indexing,test_indexing_lexys.in)
	$(call RUNCOMMAND,test_indexing,test_indexing_xyles.in)
	$(call RUNCOMMAND,test_indexing,test_indexing_yxles.in)
	$(call RUNCOMMAND,test_indexing,test_indexing_lexys_lz.in)
	$(call RUNCOMMAND,test_indexing,test_indexing_xyles_lz.in)
	$(call RUNCOMMAND,test_indexing,test_indexing_yxles_lz.in)
endif
endif

test_indexing.o:  ${TEST_MODS} 

test_indexing: test_indexing.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

