
sinclude ../Makefile.tests_common

clean_local_test:
	-$(QUIETSYM)rm -f wstar_units
clean_individual_test: clean_local_test

ifdef TESTNORUN
all: wstar_units 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: wstar_units 
	$(call RUNCOMMAND,wstar_units,wstar_units.list)
endif
endif

wstar_units.o:  ${TEST_MODS} 

wstar_units: wstar_units.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

