
sinclude ../Makefile.tests_common

clean_local_test:
	-$(QUIETSYM)rm -f test_theta_grid
clean_individual_test: clean_local_test

ifdef TESTNORUN
all: test_theta_grid
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: test_theta_grid
	$(call RUNCOMMAND,test_theta_grid,)
endif
endif

test_theta_grid.o:  ${TEST_MODS} 

test_theta_grid: test_theta_grid.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

