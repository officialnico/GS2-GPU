PFUNIT_TEST_EXE=pfunit_test_suite
PFUNIT_TEST_DIR=$(GK_HEAD_DIR)/tests/pfunit_tests
PFUNIT_NPROC=4

PFUNIT_VERSION_MAJOR ?= 3

MACRO_DEFS += -DPFUNIT_VERSION_MAJOR=$(PFUNIT_VERSION_MAJOR)

ifeq ($(PFUNIT_VERSION_MAJOR),4)
  include $(PFUNIT_DIR)/include/PFUNIT.mk

  PFUNIT_PREPROCESSOR ?= $(PFUNIT_BIN_DIR)/funitproc
  PFUNIT_INC := -I$(PFUNIT_INCLUDE_DIR)
  PFUNIT_LIB := -L$(PFUNIT_TOP_DIR)/lib -lpfunit -lfunit $(FARGPARSE_LIBRARIES)
  F90FLAGS += $(PFUNIT_EXTRA_FFLAGS)
else
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_DIR)/bin/pFUnitParser.py
  PFUNIT_INC := -I$(PFUNIT_DIR)/mod
  PFUNIT_LIB := -L$(PFUNIT_DIR)/lib -lpfunit
endif

ifeq ($(GK_SYSTEM),bitbucket)
  PFUNIT_TEST_PATH=/opt/atlassian/pipelines/agent/build/pfunit_tests
else
  PFUNIT_TEST_PATH=.
endif


GS2_OBJS = $(filter-out gs2.o, $(gs2_mod))

PFUNIT_TESTS= $(shell find $(PFUNIT_TEST_DIR) -name '*.pf')

PFUNIT_OBJS = $(PFUNIT_TESTS:.pf=.o)
PFUNIT_INTERMEDIATES = $(PFUNIT_TESTS:.pf=.F90)
PFUNIT_DRIVER = $(PFUNIT_DIR)/include/driver.F90
PFUNIT_TEST_LIST = $(PFUNIT_TEST_DIR)/testSuites.inc

.INTERMEDIATE: $(PFUNIT_INTERMEDIATES)

.PHONY: pfunit_tests make_test run_test check_pfunit_install pfunit_coverage run_gcov pfunit_clean

make_test: $(PFUNIT_TEST_EXE)

$(PFUNIT_DRIVER): $(PFUNIT_TEST_LIST)
	@touch $(PFUNIT_DRIVER)

$(PFUNIT_TEST_EXE): $(GS2_OBJS) $(PFUNIT_DRIVER) $(PFUNIT_OBJS)
	$(call messageOp,"Compiling","\\t:","pFUnit tests")
	$(QUIETSYM)$(FC) -o $@ $^ -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_LIB) -I. $(SIMPLEDATAIO_LIB_ABS) $(LIBS) $(F90FLAGS) -I$(PFUNIT_TEST_DIR) -D_TEST_SUITES='"$(PFUNIT_TEST_LIST)"'

run_test: make_test
	-@$(RM) pfunit.log pfunit.error
	$(call messageOp,"Running ","\\t:","pFUnit tests")
	$(QUIETSYM)(touch pfunit.log ; tail -f pfunit.log &  mpirun -np $(PFUNIT_NPROC) $(PFUNIT_TEST_PATH)/$(PFUNIT_TEST_EXE) -o pfunit.log > /dev/null 2> pfunit.error ) || (cat pfunit.error ; exit 1)

run_gcov:
	gcov *.f90

pfunit_tests: check_pfunit_install make_test run_test

check_pfunit_install:
ifeq ("$(wildcard $(PFUNIT_PREPROCESSOR))","")
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
	$(QUIETSYM)$(RM) $(PFUNIT_TEST_EXE) $(PFUNIT_TEST_DIR)/*.o $(PFUNIT_TEST_DIR)/*.mod \
      $(PFUNIT_TEST_DIR)/*~ $(PFUNIT_TEST_DIR)/*.gcov $(PFUNIT_TEST_DIR)/*.gcda \
      $(PFUNIT_TEST_DIR)/*.gcno $(PFUNIT_TEST_DIR)/*.log $(PFUNIT_TEST_DIR)/*.error \
      $(PFUNIT_OBJS) $(PFUNIT_INTERMEDIATES) ./test_*mod ./wraptest_*mod
