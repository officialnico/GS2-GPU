sinclude ../Makefile.tests_common

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
# `mpirun`). So let's use a copy of `TESTCOMMAND` to run `gs2`, and
# then we can undefine `TESTCOMMAND`
RESTART_EXEC := $(TESTCOMMAND)
TESTCOMMAND =

# Make a symlink of gs2 in the current directory, needed because
# `TESTCOMMAND` ends with `./`, and `BINDIR` is an absolute path, so
# we can't just do `$(RESTART_EXEC)$(BINDIR)/gs2`
clean_individual_test: TEST_NAME = restarting
link_gs2: clean_individual_test
	@ln -sfv $(BINDIR)/gs2 ./gs2

# Also remove the gs2 symlink
CLEANCOMMAND += ./gs2

all: link_gs2
	$(call RUNCOMMAND,./test_restarting,--gs2="$(RESTART_EXEC)gs2")

.PHONY: all
endif
else
all:
	$(call SKIPCOMMAND,compiled without netCDF; compile with USE_NETCDF=on to run this test)
endif
endif
