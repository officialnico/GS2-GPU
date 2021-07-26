
sinclude ../Makefile.tests_common

clean_individual_test: LOCAL_TEST_FILES:=test_fields_response_dumping gs2 OUTPUT_TEST

ifdef TESTNORUN
all: test_fields_response_dumping
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else

# Used to check the values in the netcdf output agree
define TEST_FIELDS_RESPONSE_DUMPING_TEST
$(call messageOp,"Running ","\\t:",$(1))
@./compare.py $(2) $(3) > OUTPUT_TEST ; if test "$$?" != 0; then $(FAILCOMMAND) ; else $(PASSCOMMAND) ;  fi
endef

# Used to check that there were no complaints about not being able to read the file.
define TEST_FIELDS_RESPONSE_CHECK_READ
$(call messageOp,"Checking ","\\t:",$(1))
@grep "Could not find response file" $(1).error ; if test "$$?" == 0; then $(FAILCOMMAND) ; else $(PASSCOMMAND) ;  fi
endef

ifndef USE_NETCDF
all:
	$(call SKIPCOMMAND,compiled without NetCDF library; compile with USE_NETCDF=on to run this test)
else
all: ${BINDIR}/gs2
	@$(MKDIR) local
	@$(MKDIR) implicit
	@ln -sf ${BINDIR}/gs2 ./gs2
# Run the implicit case as the baseline
	@${TESTCOMMAND}gs2 test_fields_response_dumping_implicit_dump.in > OUTPUT
# Run with the local field solver
	@${TESTCOMMAND}gs2 test_fields_response_dumping_local_dump.in >> OUTPUT
	$(call TEST_FIELDS_RESPONSE_DUMPING_TEST,test_fields_response_dumping,test_fields_response_dumping_implicit_dump.out.nc,test_fields_response_dumping_local_dump.out.nc)
# Implicit read
	@${TESTCOMMAND}gs2 test_fields_response_dumping_implicit_read.in >> OUTPUT
	$(call TEST_FIELDS_RESPONSE_DUMPING_TEST,test_fields_response_dumping,test_fields_response_dumping_implicit_dump.out.nc,test_fields_response_dumping_implicit_read.out.nc)
	$(call TEST_FIELDS_RESPONSE_CHECK_READ,test_fields_response_dumping_implicit_read)
# Local read
	@${TESTCOMMAND}gs2 test_fields_response_dumping_local_read.in >> OUTPUT
	$(call TEST_FIELDS_RESPONSE_DUMPING_TEST,test_fields_response_dumping,test_fields_response_dumping_implicit_dump.out.nc,test_fields_response_dumping_local_read.out.nc)
	$(call TEST_FIELDS_RESPONSE_CHECK_READ,test_fields_response_dumping_local_read)
# Implicit read from local
	@${TESTCOMMAND}gs2 test_fields_response_dumping_implicit_read_from_local.in >> OUTPUT
	$(call TEST_FIELDS_RESPONSE_DUMPING_TEST,test_fields_response_dumping,test_fields_response_dumping_implicit_dump.out.nc,test_fields_response_dumping_implicit_read_from_local.out.nc)
	$(call TEST_FIELDS_RESPONSE_CHECK_READ,test_fields_response_dumping_implicit_read_from_local)
# Local read from implicit
	@${TESTCOMMAND}gs2 test_fields_response_dumping_local_read_from_implicit.in >> OUTPUT
	$(call TEST_FIELDS_RESPONSE_DUMPING_TEST,test_fields_response_dumping,test_fields_response_dumping_implicit_dump.out.nc,test_fields_response_dumping_local_read_from_implicit.out.nc)
	$(call TEST_FIELDS_RESPONSE_CHECK_READ,test_fields_response_dumping_local_read_from_implicit)
	@-rm -rf local/ implicit/
endif
endif
endif
