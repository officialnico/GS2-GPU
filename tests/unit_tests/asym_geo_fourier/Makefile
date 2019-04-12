
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_asym_geo_fourier

ifdef TESTNORUN
all: test_asym_geo_fourier
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_asym_geo_fourier
	$(call RUNCOMMAND,test_asym_geo_fourier,)
endif
endif

ogyropsi.dat EQDSK:
	tar -xjf data.tar.bz2

test_asym_geo_fourier.o:  ${TEST_MODS} 

test_asym_geo_fourier: test_asym_geo_fourier.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

