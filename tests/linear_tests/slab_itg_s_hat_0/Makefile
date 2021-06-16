
sinclude ../Makefile.tests_common

clean_local_test:
	-$(QUIETSYM)rm -f slab_itg_s_hat_0_low_res gs2_slab_itg
clean_individual_test: clean_local_test

COPY=cp ../../../bin/gs2 gs2_slab_itg

ifdef TESTNORUN
all: slab_itg_s_hat_0_low_res 
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: slab_itg_s_hat_0_low_res 
	$(call RUNCOMMAND,slab_itg_s_hat_0_low_res,slab_itg_s_hat_0_low_res.in)
	$(call COPY)
	$(call RUNCOMMAND,gs2_slab_itg,slab_itg_s_hat_0_low_res.in)
endif
endif

slab_itg_s_hat_0_low_res.o:  ${TEST_MODS} 

slab_itg_s_hat_0_low_res: slab_itg_s_hat_0_low_res.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

