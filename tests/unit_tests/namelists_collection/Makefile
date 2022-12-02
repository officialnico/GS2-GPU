sinclude ../Makefile.tests_common
include $(GK_HEAD_DIR)/Makefile.target_gs2

PYTEST := $(shell which pytest)
GK_HAS_PYTEST := $(if $(PYTEST),yes,no)

ifdef TESTNORUN
all:
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})

ifeq ($(GK_HAS_PYTEST),yes)
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else # TEST_LEVEL

# We want to run 'pytest' without MPI, so undefine TESTCOMMAND
TESTCOMMAND =

all: $(BINDIR)/write_default_input_file $(BINDIR)/gs2
	$(call RUNCOMMAND,PYTHONPATH=$(GK_HEAD_DIR)/scripts pytest)

.PHONY: all

endif # TEST_LEVEL
else # PYTEST
all:
	$(call SKIPCOMMAND,Could not find 'pytest'; please install with 'make create-test-virtualenv')
endif # PYTEST
endif # TESTNORUN
