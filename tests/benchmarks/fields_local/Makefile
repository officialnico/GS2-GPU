
sinclude ../Makefile.tests_common

INPUTFILE=time_fields_local_level_1.in
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=time_fields_local_level_2.in
endif

all: time_fields_local
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_fields_local,$(INPUTFILE)) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_fields_local.o:  ${TEST_MODS} 

time_fields_local: time_fields_local.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

