PFUNIT_TEST_EXE=pfunit_test_suite

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

PFUNIT_TESTS= $(shell find ./ -name '*.pf')

PFUNIT_OBJS = $(PFUNIT_TESTS:.pf=.o)
PFUNIT_DRIVER = $(PFUNIT_DIR)/include/driver.F90

make_test: $(GS2_OBJS)
	$(QUIETSYM)$(RM) $(PFUNIT_TEST_DIR)/$(PFUNIT_TEST_EXE)
	$(QUIETSYM)$(MAKE) -C $(PFUNIT_TEST_DIR) $(PFUNIT_TEST_EXE)

$(PFUNIT_TEST_EXE): $(PFUNIT_DRIVER) $(PFUNIT_OBJS) testSuites.inc
	$(call messageOp,"Compiling","\\t:","pFUnit tests")
	$(QUIETSYM)$(FC) -o $@ $(GS2_OBJS) -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_DRIVER) $(PFUNIT_OBJS) $(PFUNIT_LIB) -I. $(SIMPLEDATAIO_LIB_ABS) $(LIBS) $(F90FLAGS)

run_test: make_test
	$(call messageOp,"Running ","\\t:","pFUnit tests")
	$(QUIETSYM)cd $(PFUNIT_TEST_DIR) && mpirun -np 4 $(PFUNIT_TEST_PATH)/$(PFUNIT_TEST_EXE) -robust -xml $(GK_HEAD_DIR)/test-results/pfunit.xml -o pfunit.log > /dev/null 2> pfunit.error && cat pfunit.log || cat pfunit.error

run_gcov:
	gcov *.f90

ifndef USE_PFUNIT
pfunit_tests: test_fail_target
test_fail_target:
	$(error Cannot make pfunit_tests as USE_PFUNIT undefined. Rebuild everything with USE_PFUNIT=on)
else
pfunit_tests: make_test run_test
endif

ifndef USE_PFUNIT
pfunit_coverage: coverage_pfunit_fail
coverage_pfunit_fail:
	$(error Cannot make pfunit_coverage as USE_PFUNIT undefined. Rebuild everything with USE_PFUNIT=on USE_GCOV=on)
else
ifndef USE_GCOV
pfunit_coverage: coverage_gcov_fail
coverage_gcov_fail:
	$(error Cannot make pfunit_coverage as USE_GCOV undefined. Rebuild everything with USE_PFUNIT=on USE_GCOV=on)
else
pfunit_coverage: gs2 make_test run_test run_gcov
	mkdir -p coverage-pfunit
	mv *.gcov coverage-pfunit/.
endif
endif

pfunit_clean:
	$(call messageRm,"pFUnit tests")
	$(QUIETSYM)$(RM) $(PFUNIT_TEST_EXE) *.o *.mod *~ *.gcov *.gcda *.gcno *.log *.error */*.o */*.mod */*~ */*.gcov */*.gcda */*.gcno
