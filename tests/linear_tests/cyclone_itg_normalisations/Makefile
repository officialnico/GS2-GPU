
sinclude ../Makefile.tests_common
TEST_BASE_NAME:=cyclone_normalisations
clean_individual_test: LOCAL_TEST_FILES:= ${TEST_BASE_NAME} *.vspace_*
build_test_exe: ${TEST_BASE_NAME}
ifdef TESTNORUN
all: build_test_exe
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: build_test_exe
ifeq "$(shell expr ${TEST_LEVEL} \<= 2)" "1"
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_baseline.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_electron_mass.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_high_temperature.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_density.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_charge.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_ions_length.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_baseline.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_electron_mass.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_high_temperature.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_density.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_charge.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_two_species_length.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_baseline.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_electron_mass.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_high_temperature.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_density.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_charge.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisionless_beta_length.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_baseline.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_electron_mass.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_high_temperature.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_density.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_charge.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_length.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_all.in)
else
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_baseline.in)
	$(call RUNCOMMAND,${TEST_BASE_NAME},collisional_all.in)
endif
endif
endif
${TEST_BASE_NAME}.o:  ${TEST_MODS}

${TEST_BASE_NAME}: ${TEST_BASE_NAME}.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

update_golden_answers:
	$(call messageOp,"Updating golden answers for","\\t:",${TEST_BASE_NAME},"","")
	$(QUIETSYM)$(GK_HEAD_DIR)/scripts/update_golden_answers.sh ${TEST_BASE_NAME}
