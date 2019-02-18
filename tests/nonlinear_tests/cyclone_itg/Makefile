
sinclude ../Makefile.tests_common

CLEANCOMMAND+= cyclone_itg_low_res

#Define this tests level
MY_TEST_LEVEL=2
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
ifndef USE_FFT
all:
	$(call SKIPCOMMAND,compiled without FFT library; compile with USE_FFT=on to run this test)
else
all: cyclone_itg_low_res 
	$(call RUNCOMMAND,cyclone_itg_low_res,cyclone_itg_low_res.in)
endif
endif

cyclone_itg_low_res.o:  ${TEST_MODS} 

cyclone_itg_low_res: cyclone_itg_low_res.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

