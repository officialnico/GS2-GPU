
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=test_kinetic_energy_transfer
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} gs2
build_test_exe: ${BINDIR}/gs2
ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
define TEST_KINETIC_ENERGY
$(call messageOp,"Running ","\\t:",$(1))
@./run_tests.py --file $(2) --verbosity 0 >> OUTPUT ; if test "$$?" != 0; then $(FAILCOMMAND) ; else $(PASSCOMMAND) ;  fi
endef

ifndef USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test)
else
all: build_test_exe
	@ln -sf ${BINDIR}/gs2 ./${TEST_BASE_NAME}
	$(call messageOp,"Running ","\\t:",${TESTCOMMAND}${TEST_BASE_NAME} ${TEST_BASE_NAME}.in)
	@${TESTCOMMAND}${TEST_BASE_NAME} ${TEST_BASE_NAME}.in > OUTPUT
	$(call TEST_KINETIC_ENERGY,${TEST_BASE_NAME},${TEST_BASE_NAME}.out.nc)
endif
endif
endif
