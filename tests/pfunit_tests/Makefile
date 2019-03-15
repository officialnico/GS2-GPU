PFUNIT_TEST_EXE=pfunit_test_suite
PFUNIT_TEST_DIR=$(GK_HEAD_DIR)/tests/pfunit_tests

PFUNIT_LIB=-L$(PFUNIT_DIR)/lib -lpfunit
PFUNIT_INC=-I$(PFUNIT_DIR)/mod

ifeq ($(GK_SYSTEM),bitbucket)
  PFUNIT_TEST_PATH=/opt/atlassian/pipelines/agent/build/pfunit_tests
else
  PFUNIT_TEST_PATH=.
endif

%.o : %.F90
	$(FC) $(F90FLAGS) -c $< $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@

%.F90: %.pf
	$(PFUNIT_DIR)/bin/pFUnitParser.py $< $@

GS2_OBJS = $(filter-out gs2.o, $(gs2_mod))

PFUNIT_TESTS= $(shell find $(PFUNIT_TEST_DIR) -name '*.pf')

PFUNIT_OBJS = $(PFUNIT_TESTS:.pf=.o)
PFUNIT_DRIVER = $(PFUNIT_DIR)/include/driver.F90
PFUNIT_TEST_LIST = $(PFUNIT_TEST_DIR)/testSuites.inc



make_test: $(PFUNIT_TEST_EXE)

$(PFUNIT_DRIVER): $(PFUNIT_TEST_LIST)
	@touch $(PFUNIT_DRIVER)

$(PFUNIT_TEST_EXE): $(PFUNIT_DRIVER) $(PFUNIT_OBJS) $(GS2_OBJS)
	$(call messageOp,"Compiling","\\t:","pFUnit tests")
	$(QUIETSYM)$(FC) -o $@ $^ -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_LIB) -I. $(SIMPLEDATAIO_LIB_ABS) $(LIBS) $(F90FLAGS) -I$(PFUNIT_TEST_DIR)

run_test: make_test
	$(call messageOp,"Running ","\\t:","pFUnit tests")
	$(QUIETSYM)cd $(PFUNIT_TEST_DIR) && mpirun -np 4 $(PFUNIT_TEST_PATH)/$(PFUNIT_TEST_EXE) -robust -xml $(GK_HEAD_DIR)/test-results/pfunit.xml -o pfunit.log > /dev/null 2> pfunit.error && cat pfunit.log || cat pfunit.error

run_gcov:
	gcov *.f90

pfunit_tests: check_pfunit_install make_test run_test

check_pfunit_install:
ifeq ("$(wildcard $(PFUNIT_DIR)/bin/pFUnitParser.py)","")
	$(error pFUnit files not found. Make sure pFUnit is installed and PFUNIT_DIR is set to its location)
endif


ifndef USE_GCOV
pfunit_coverage: coverage_gcov_fail
coverage_gcov_fail:
	$(error Cannot make pfunit_coverage as USE_GCOV undefined. Rebuild everything with USE_GCOV=on)
else
pfunit_coverage: gs2 make_test run_test run_gcov
	mkdir -p coverage-pfunit
	mv *.gcov coverage-pfunit/.
endif

pfunit_clean:
	$(call messageRm,"pFUnit tests")
	$(QUIETSYM)$(RM) $(PFUNIT_TEST_EXE) *.o *.mod *~ *.gcov *.gcda *.gcno *.log *.error */*.o */*.mod */*~ */*.gcov */*.gcda */*.gcno
