PFUNIT_TEST_EXE=pfunit_test_suite
PFUNIT_TEST_DIR=$(GK_HEAD_DIR)/tests/pfunit_tests
PFUNIT_NPROC=4

# Default major version of pFUnit
# Versions 3 and 4 require different flags, etc.
PFUNIT_VERSION_MAJOR ?= 3
# Installation directory of pFUnit
# Default is the directory where `scripts/build_pfunit` puts it
PFUNIT_DIR ?= $(GK_HEAD_DIR)/externals/pFUnit_build/install

PFUNIT_LOG:=pfunit.log
PFUNIT_ERROR:=pfunit.error

MACRO_DEFS += -DPFUNIT_VERSION_MAJOR=$(PFUNIT_VERSION_MAJOR)

ifeq ($(PFUNIT_VERSION_MAJOR),4)
  PFUNIT_VERSION_MINOR ?= 8
  sinclude $(PFUNIT_DIR)/PFUNIT-${PFUNIT_VERSION_MAJOR}.${PFUNIT_VERSION_MINOR}/include/PFUNIT.mk
  PFUNIT_MOD := pfunit
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_BIN_DIR)/funitproc
  PFUNIT_INC := -I$(PFUNIT_INCLUDE_DIR)
  PFUNIT_LIB := -L$(PFUNIT_TOP_DIR)/lib -lpfunit -lfunit $(FARGPARSE_LIBRARIES)
  F90FLAGS += $(PFUNIT_EXTRA_FFLAGS)
else
  PFUNIT_VERSION_MINOR ?= 3
  PFUNIT_TOP_DIR := $(PFUNIT_DIR)
  PFUNIT_MOD := pfunit_mod
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_DIR)/bin/pFUnitParser.py
  PFUNIT_INC := -I$(PFUNIT_DIR)/mod
  PFUNIT_LIB := -L$(PFUNIT_DIR)/lib -lpfunit
endif

MACRO_DEFS += -DPFUNIT_MOD=$(PFUNIT_MOD)

PFUNIT_TEST_PATH :=.

GS2_OBJS := $(filter-out $(OBJDIR)/gs2.o, $(gs2_mod))

PFUNIT_TESTS:= $(shell find $(PFUNIT_TEST_DIR) -name '[^.\#]*.pf')

PFUNIT_OBJS := $(PFUNIT_TESTS:.pf=.o)
PFUNIT_INTERMEDIATES := $(PFUNIT_TESTS:.pf=.F90)
PFUNIT_DRIVER := $(PFUNIT_TOP_DIR)/include/driver.F90
PFUNIT_TEST_LIST := $(PFUNIT_TEST_DIR)/testSuites.inc

.INTERMEDIATE: $(PFUNIT_INTERMEDIATES)

.PHONY: pfunit_tests make_test run_test check_pfunit_install pfunit_coverage run_gcov pfunit_clean
.PHONY: build_pfunit_library clean_pfunit_library



##################################################
# Handle flags to pass to scripts/build_pfunit
##################################################
# How many jobs to use
BUILD_PFUNIT_JOBS ?= 1
# Where do we install to
BUILD_PFUNIT_INSTALL_DIR ?= $(GK_HEAD_DIR)/externals/pFUnit_build/install
# Where do we build
BUILD_PFUNIT_BUILD_DIR ?=
# Do we ignore the git submodule
BUILD_PFUNIT_IGNORE_GIT ?=
# Do we force the install
BUILD_PFUNIT_FORCE_INSTALL ?=

# Construct full set of flags
BUILD_PFUNIT_FLAGS ?=
BUILD_PFUNIT_FLAGS += -j ${BUILD_PFUNIT_JOBS} -v ${PFUNIT_VERSION_MAJOR} -i ${BUILD_PFUNIT_INSTALL_DIR} -V ${PFUNIT_VERSION_MINOR}

ifndef USE_MPI
BUILD_PFUNIT_FLAGS += -m
endif

ifdef USE_OPENMP
BUILD_PFUNIT_FLAGS += -o
endif

ifdef BUILD_PFUNIT_BUILD_DIR
BUILD_PFUNIT_FLAGS += -b ${BUILD_PFUNIT_BUILD_DIR}
endif

ifdef BUILD_PFUNIT_IGNORE_GIT
BUILD_PFUNIT_FLAGS += -g
endif

ifdef BUILD_PFUNIT_FORCE_INSTALL
BUILD_PFUNIT_FLAGS += -f
endif

test_build_pfunit_library:
	@echo scripts/build_pfunit ${BUILD_PFUNIT_FLAGS}

build_pfunit_library:
	$(call messageInfo,"Building pFUnit in location: $(BUILD_PFUNIT_INSTALL_DIR)")
	$(QUIETSYM)scripts/build_pfunit ${BUILD_PFUNIT_FLAGS}

clean_pfunit_library:
	$(call messageInfo "Removing pFUnit library")
	$(QUIETSYM)scripts/build_pfunit -r ${BUILD_PFUNIT_FLAGS}

make_test: $(PFUNIT_TEST_EXE)

$(PFUNIT_DRIVER): $(PFUNIT_TEST_LIST)
	@touch $(PFUNIT_DRIVER)

$(PFUNIT_OBJS): $(GS2_OBJS)

$(PFUNIT_TEST_EXE): $(PFUNIT_OBJS) $(PFUNIT_DRIVER) $(GS2_OBJS)
	$(call messageOp,"Compiling","\\t:","pFUnit tests")
	$(QUIETSYM)$(FC) -o $@ $^ -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_LIB) -I. $(LIBS) $(F90FLAGS) -I$(PFUNIT_TEST_DIR) -D_TEST_SUITES='"$(PFUNIT_TEST_LIST)"'

run_test: make_test
	-@rm -f pfunit.log pfunit.error
	$(call messageOp,"Running ","\\t:","pFUnit tests")
	$(QUIETSYM) touch $(PFUNIT_LOG); \
      mpirun -np $(PFUNIT_NPROC) $(PFUNIT_TEST_PATH)/$(PFUNIT_TEST_EXE) -o $(PFUNIT_LOG) > /dev/null 2> $(PFUNIT_ERROR) \
      && cat $(PFUNIT_LOG) \
      || (cat $(PFUNIT_LOG) $(PFUNIT_ERROR); exit 1)

run_gcov: run_test
#Note we can't just pass $(OBJDIR)/*.o as not all objects are used/generate gcda
#files, which gcov now indicates as an error, causing make to fail.
	gcov -p -o $(OBJDIR) $(patsubst %.gcda, %.o, $(wildcard $(addsuffix /*.gcda, $(OBJDIR))))

pfunit_tests: check_pfunit_install make_test run_test

check_pfunit_install:
ifeq ("$(wildcard $(PFUNIT_PREPROCESSOR))","")
	$(error \
pFUnit files not found. Make sure pFUnit is installed and PFUNIT_DIR is set to its location.\
Consider running 'make build_pfunit_library' or use 'scripts/build_pfunit' \
to try to build a copy of pFUnit if you do not have one)
endif


ifndef USE_GCOV
pfunit_coverage: coverage_gcov_fail
coverage_gcov_fail:
	$(error Cannot make pfunit_coverage as USE_GCOV undefined. Rebuild everything with USE_GCOV=on)
else
pfunit_coverage: gs2 make_test run_test run_gcov
	mkdir -p coverage-pfunit
	mv *.gcov coverage-pfunit/. || echo "Couldn't move generated gcov files to coverage-pfunit. This is likely because gcov files not found in ${GK_HEAD_DIR}."
endif

distclean: pfunit_clean
pfunit_clean:
	$(call messageRm,"pFUnit tests")
	$(QUIETSYM)$(RM) $(PFUNIT_TEST_EXE) $(PFUNIT_TEST_DIR)/*.o $(PFUNIT_TEST_DIR)/*.mod \
      $(PFUNIT_TEST_DIR)/*~ $(PFUNIT_TEST_DIR)/*.gcov $(PFUNIT_TEST_DIR)/*.gcda \
      $(PFUNIT_TEST_DIR)/*.gcno $(PFUNIT_TEST_DIR)/*.log $(PFUNIT_TEST_DIR)/*.error \
      $(PFUNIT_OBJS) $(PFUNIT_INTERMEDIATES) ./test_*mod ./wraptest_*mod
