
sinclude ../Makefile.tests_common

INPUTFILE=time_ffts_level_1.in
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=time_ffts_level_2.in
endif

all: time_ffts
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_ffts,$(INPUTFILE)) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_ffts.o:  ${TEST_MODS} 

time_ffts: time_ffts.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

