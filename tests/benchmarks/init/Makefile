
sinclude ../Makefile.tests_common

INPUTFILE=time_init_level_1.in
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=time_init_level_2.in
endif

all: time_init 
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_init,$(INPUTFILE)) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_init.o:  ${TEST_MODS} 

time_init: time_init.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

