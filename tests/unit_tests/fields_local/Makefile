
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_fields_local test_fields_local_tuning

ifdef TESTNORUN
all: test_fields_local test_fields_local_tuning 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: test_fields_local test_fields_local_tuning 
	$(call RUNCOMMAND,test_fields_local,test_fields_local.in)
	$(call RUNCOMMAND,test_fields_local_tuning,test_fields_local_tuning.in)
endif
endif

test_fields_local.o:  ${TEST_MODS} 

test_fields_local: test_fields_local.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

test_fields_local_tuning.o:  ${TEST_MODS}

test_fields_local_tuning: test_fields_local_tuning.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

