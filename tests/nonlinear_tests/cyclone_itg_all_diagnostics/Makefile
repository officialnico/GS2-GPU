sinclude ../Makefile.tests_common

TEST_BASE_NAME:=test_cyclone_all_diagnostics
build_test_exe: ${TEST_BASE_NAME}
clean_individual_test: LOCAL_TEST_FILES:=*.{error,exit_reason,fort*,g,gridgen*,nc*,optim,out*,cres,db,yxdist,dump*,antot,gs,dc2,check*,OLD,NEW}

ifdef TESTNORUN
all: build_test_exe
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

all: build_test_exe
	$(call RUNCOMMAND,../../../bin/gs2,cyclone_itg_low_res_old.in)
	$(call RUNCOMMAND,../../../bin/gs2,cyclone_itg_low_res_new.in)
endif
endif
endif

${TEST_BASE_NAME}: $(TEST_MODS)
