
sinclude ../Makefile.tests_common

INPUTFILE=time_timestep_level_1.in
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=time_timestep_level_2.in
endif

all: time_timestep 
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_timestep,$(INPUTFILE)) ; \
	done
	#$(TESTCOMMAND) time_timestep time_timestep.in

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_timestep.o:  ${TEST_MODS} 

time_timestep: time_timestep.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

