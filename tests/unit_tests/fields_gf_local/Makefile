
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_fields_gf_local

ifdef TESTNORUN
all: test_fields_gf_local 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
ifeq ($(shell expr ${NTESTPROCS} \< 6; echo $?),1)
all:
	$(call SKIPCOMMAND,current NTESTPROCS=${NTESTPROCS}; set NTESTPROCS>=6 to run this test)
else
all: test_fields_gf_local 
	$(call RUNCOMMAND,test_fields_gf_local,test_fields_gf_local.in)
endif
endif
endif

test_fields_gf_local.o:  ${TEST_MODS} 

test_fields_gf_local: test_fields_gf_local.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

