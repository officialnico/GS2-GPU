
sinclude ../Makefile.tests_common

CLEANCOMMAND+= slab_itg_s_hat_0_low_res

COPY=cp ../../../gs2 gs2_slab_itg

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: slab_itg_s_hat_0_low_res 
	$(call RUNCOMMAND,slab_itg_s_hat_0_low_res,slab_itg_s_hat_0_low_res.in)
	$(call COPY)
	$(call RUNCOMMAND,gs2_slab_itg,slab_itg_s_hat_0_low_res.in)
endif

slab_itg_s_hat_0_low_res.o:  ${TEST_MODS} 

slab_itg_s_hat_0_low_res: slab_itg_s_hat_0_low_res.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

