
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_nonlinear_terms

ifdef TESTNORUN
all: test_nonlinear_terms 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
ifndef USE_FFT
all:
	$(call SKIPCOMMAND,compiled without FFT library; compile with USE_FFT=on to run this test)
else
all: test_nonlinear_terms 
	$(call RUNCOMMAND,test_nonlinear_terms,test_nonlinear_terms.in)
endif
endif
endif

test_nonlinear_terms.o:  ${TEST_MODS} 

test_nonlinear_terms: test_nonlinear_terms.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

