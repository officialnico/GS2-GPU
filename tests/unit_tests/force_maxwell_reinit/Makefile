
sinclude ../Makefile.tests_common

#export LOCAL_TEST_FILES

clean_individual_test: LOCAL_TEST_FILES:=HELLO *.{error,exit_reason,fort*,g,gridgen*,nc*,optim,out*} OUTPUT HELLO

PYTHONCOMPARE:=./compare.py initial.out.nc phi restart.out.nc phi

define LOCALFAILCOMMAND
echo "FAILED: `basename $$PWD` ($(PYTHONCOMPARE))" >> results_of_test.txt ; $(TEST_DIR)/write_test_xml.sh `basename $$PWD` 1 0 1 0 $(TEST_RESULT_DIR); cat OUTPUT
endef

define LOCALPASSCOMMAND
echo "passed: `basename $$PWD` ($(PYTHONCOMPARE))" >> results_of_test.txt ; $(TEST_DIR)/write_test_xml.sh `basename $$PWD` 1 0 0 0 $(TEST_RESULT_DIR)
endef

ifdef TESTNORUN
all: test_force_maxwell_reinit 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
ifndef USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else

all: test_force_maxwell_reinit 
	@$(TESTCOMMAND)../../../bin/gs2 initial.in > OUTPUT
	@$(TESTCOMMAND)../../../bin/gs2 restart.in > OUTPUT
	@$(PYTHONCOMPARE) ; rc=$$? ; if test "$$rc" != 0; then $(LOCALFAILCOMMAND) ; else $(LOCALPASSCOMMAND) ;  fi
endif
endif
endif

test_force_maxwell_reinit: $(TEST_MODS)
