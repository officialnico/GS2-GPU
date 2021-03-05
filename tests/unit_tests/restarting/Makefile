sinclude ../Makefile.tests_common

CLEANCOMMAND = ./runtest.py --clean

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
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else

# We're calling a Python script and don't need mpirun, so override TESTCOMMAND
TESTCOMMAND =
all:
	$(call RUNCOMMAND,./test_restarting)

.PHONY: all
endif
else
all:
	$(call SKIPCOMMAND,compiled without netCDF; compile with USE_NETCDF=on to run this test)
endif
endif
