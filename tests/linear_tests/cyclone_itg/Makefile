
sinclude ../Makefile.tests_common

CLEANCOMMAND+= cyclone_itg_low_res cyclone_itg_replay

COPY=cp fields_local_single_mode.out.nc fields_local_single_mode_replay.out.nc

ifdef TESTNORUN
all: cyclone_itg_low_res cyclone_itg_replay
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: cyclone_itg_low_res cyclone_itg_replay
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_implicit_single_mode.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_single_mode.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_single_mode_nonadiabatic_init.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_single_mode_npassing_even.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_single_mode_npassing_odd.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_range.in)
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_local_box.in)
ifeq ($(shell expr ${NTESTPROCS} \< 6; echo $?),1)
	$(call SKIPCOMMAND,current NTESTPROCS=${NTESTPROCS}; set NTESTPROCS>=6 to run this test)
else
	$(call RUNCOMMAND,cyclone_itg_low_res,fields_gf_local_single_mode.in)
endif
ifdef USE_NEW_DIAG
	$(call COPY)
	$(call RUNCOMMAND,cyclone_itg_replay,fields_local_single_mode_replay.in)
else
	$(call SKIPCOMMAND,compiled without new diagnostics; compile with USE_NEW_DIAG=on to run this test)
endif
endif
endif

cyclone_itg_low_res.o cyclone_itg_replay.o:  ${TEST_MODS} 

cyclone_itg_low_res: cyclone_itg_low_res.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

cyclone_itg_replay: cyclone_itg_replay.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

