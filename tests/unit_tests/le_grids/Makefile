
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_le_grids

ifdef TESTNORUN
all: test_le_grids 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_le_grids 
	$(call RUNCOMMAND,test_le_grids,test_le_grids.in)
endif
endif

test_le_grids.o:  ${TEST_MODS} 

test_le_grids: test_le_grids.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

