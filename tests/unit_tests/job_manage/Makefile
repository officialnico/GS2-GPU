
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_job_manage_stop_file
ifdef TESTNORUN
all: test_job_manage_stop_file
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: test_job_manage_stop_file
	-$(QUIETSYM)rm -f test_job_manage.exit_reason
	$(call RUNCOMMAND,test_job_manage_stop_file,test_job_manage.in)
	-$(QUIETSYM)rm -f test_job_manage.stop
endif
endif

test_job_manage_stop_file.o:  ${TEST_MODS} 

test_job_manage_stop_file: test_job_manage_stop_file.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
