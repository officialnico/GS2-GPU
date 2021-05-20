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

# This test is a little different from most of the others: the main
# test executable is a Python script that calls `gs2`. What we want is
# to call `gs2` with `mpirun` (or equivalent) and not the Python
# script. The usual `RUNCOMMAND` test comand calls its first argument
# with `TESTCOMMAND`, which in turn calls `TESTEXEC` (which is usually
# `mpirun`). This would end up calling the Python script with
# `mpirun`, rather than `gs2`. We could do:
#
#   RESTART_EXEC := $(TESTCOMMAND)
#   TESTCOMMAND =
#
# which allows to pass `TESTCOMMAND` to the `--gs2` test script
# argument. However, `TESTCOMMAND` ends with `./` and `BINDIR` is an
# absolute path. All of this is to say, we need to undefine
# `TESTCOMMAND` and reimplement it in a separate variable.
TESTCOMMAND =
ifdef TESTEXEC
RESTART_EXEC=$(call TESTEXEC,$(NTESTPROCS))
else
ifdef USE_MPI
RESTART_EXEC=mpirun -np $(NTESTPROCS)
else
RESTART_EXEC=./
endif
endif

all:
	$(call RUNCOMMAND,./test_restarting,--gs2="$(RESTART_EXEC) $(BINDIR)/gs2")

.PHONY: all
endif
else
all:
	$(call SKIPCOMMAND,compiled without netCDF; compile with USE_NETCDF=on to run this test)
endif
endif
