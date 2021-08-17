
sinclude ../Makefile.tests_common

clean_individual_test: LOCAL_TEST_FILES:=test_le_grids

ifdef TESTNORUN
all: test_le_grids 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: test_le_grids 
	$(call RUNCOMMAND,test_le_grids,test_le_grids.in)
endif
endif

test_le_grids.o:  ${TEST_MODS} 

test_le_grids: test_le_grids.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

