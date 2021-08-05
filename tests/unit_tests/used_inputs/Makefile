sinclude ../Makefile.tests_common

PYTEST := $(shell which pytest)
GK_HAS_PYTEST := $(if $(PYTEST),yes,no)

ifdef TESTNORUN
all:
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})

ifeq ($(GK_HAS_PYTEST),yes)
ifdef USE_NETCDF
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else # TEST_LEVEL

# We want to run 'pytest' without MPI, so undefine TESTCOMMAND
TESTCOMMAND =

all:
	$(call RUNCOMMAND,pytest)

.PHONY: all

endif # TEST_LEVEL
else # USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without netCDF; compile with USE_NETCDF=on to run this test)
endif # USE_NETCDF
else # PYTEST
all:
	$(call SKIPCOMMAND,Could not find 'pytest'; please install with 'make create-test-virtualenv')
endif # PYTEST
endif # TESTNORUN
