sinclude ../Makefile.tests_common

#export LOCAL_TEST_FILES

clean_individual_test: LOCAL_TEST_FILES:=*.{error,exit_reason,fort*,g,gridgen*,nc*,optim,out*,cres,db,yxdist,dump*,antot,gs,dc2,check*,OLD,NEW}

PYTHONCOMPARE:=./compare.py initial.out.nc phi restart.out.nc phi

define LOCALFAILCOMMAND
echo "FAILED: `basename $$PWD` ($(PYTHONCOMPARE))" >> results_of_test.txt ; $(TEST_DIR)/write_test_xml.sh `basename $$PWD` 1 0 1 0 $(TEST_RESULT_DIR); cat OUTPUT
endef

define LOCALPASSCOMMAND
echo "passed: `basename $$PWD` ($(PYTHONCOMPARE))" >> results_of_test.txt ; $(TEST_DIR)/write_test_xml.sh `basename $$PWD` 1 0 0 0 $(TEST_RESULT_DIR)
endef

ifdef TESTNORUN
all: test_cyclone_all_diagnostics
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

all: test_cyclone_all_diagnostics
	@$(TESTCOMMAND)../../../bin/gs2 cyclone_itg_low_res_old.in > OUTPUT.OLD
	@$(TESTCOMMAND)../../../bin/gs2 cyclone_itg_low_res_new.in > OUTPUT.NEW
endif
endif
endif

test_cyclone_all_diagnostics: $(TEST_MODS)
