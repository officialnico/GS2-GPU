
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_job_manage_stop_file test_job_manage_avail_cpu_time

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_job_manage_stop_file test_job_manage_avail_cpu_time
	$(call RUNCOMMAND,test_job_manage_stop_file,test_job_manage.in)
	$(call RUNCOMMAND,test_job_manage_avail_cpu_time,test_job_manage.in)
	-$(QUIETSYM)rm -f test_job_manage.stop
endif

test_job_manage_stop_file.o:  ${TEST_MODS} 

test_job_manage_stop_file: test_job_manage_stop_file.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

test_job_manage_avail_cpu_time.o:  ${TEST_MODS} 

test_job_manage_avail_cpu_time: test_job_manage_avail_cpu_time.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
