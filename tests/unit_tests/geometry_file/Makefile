sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_geometry_file
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} grid.file
build_test_exe:
ifdef TESTNORUN
all:
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifdef USE_NETCDF
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
TESTCOMMAND =
all:
	$(call RUNCOMMAND,./test_geometry_file.py)
endif
else
all:
	$(call SKIPCOMMAND,compiled without NETCDF; compile with USE_NETCDF to run this test)
endif
endif
