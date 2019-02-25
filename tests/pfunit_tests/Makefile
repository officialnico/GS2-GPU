.PHONY: tests clean

TEST_DEPS?=libgs2.a functional_tests.o benchmarks.o

%.o : %.F90
	$(FC) $(F90FLAGS) -c $< $(PFUNIT_INC) -I ../../.

%.F90: %.pf
	$(PFUNIT_DIR)/bin/pFUnitParser.py $< $@

all: pfunit_tests

GS2_TMP = $(patsubst %.o,../../%.o,$(gs2_mod))
GS2_OBJS = $(filter-out ../../gs2.o, $(GS2_TMP))

SRCS = $(wildcard *.F90)
TESTS = $(wildcard *.pf)

OBJS = $(SRCS:.F90=.o) $(TESTS:.pf=.o)
DRIVER = $(PFUNIT_DIR)/include/driver.F90 

make_test: gs2 #pftest_utils
	@cd $(PFUNIT_TEST_DIR) && $(RM) tests.x
	@cd $(PFUNIT_TEST_DIR) && make tests.x

tests.x: $(DRIVER) $(OBJS) testSuites.inc
	$(FC) -o $@ $(GS2_OBJS) -DUSE_MPI=YES $(PFUNIT_INC) $(DRIVER) $(OBJS) -L$(PFUNIT_LIB) -I. ../../libsimpledataio.a $(LIBS) $(F90FLAGS)

ifeq ($(GK_SYSTEM),bitbucket)
run_test: gs2
		@cd $(PFUNIT_TEST_DIR) && mpirun -np 4 /opt/atlassian/pipelines/agent/build/pfunit_tests/tests.x -robust -xml $(GK_HEAD_DIR)/test-results/pfunit.xml
else
run_test: gs2
		@cd $(PFUNIT_TEST_DIR) && mpirun -np 4 ./tests.x -robust -xml $(GK_HEAD_DIR)/test-results/pfunit.xml
endif

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
	@$(RM) *.gcov *.gcda *.gcno
	@cd $(PFUNIT_TEST_DIR) && $(RM) *.o *.mod *.x *~ *.gcov *.gcda *.gcno

