sinclude ../Makefile.tests_common

CLEANCOMMAND += zonal_flow_residual

ifdef TESTNORUN
all: zonal_flow_residual
	$(call SKIPCOMMAND,variable TESTNORUN is defined; unset TESTNORUN to run this test)
else
MY_TEST_LEVEL=1
VAL := $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,current TEST_LEVEL=${TEST_LEVEL}; set TEST_LEVEL>=${MY_TEST_LEVEL} to run this test)
else
all: zonal_flow_residual
	$(call RUNCOMMAND,zonal_flow_residual,zonal_flow_residual.in)
endif
endif

zonal_flow_residual.o:  ${TEST_MODS}

zonal_flow_residual: zonal_flow_residual.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$(call processName, $^),"\\n\\t\-\-\-\-\-\-\>",$(call processName, $@))
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
