
sinclude ../Makefile.tests_common

INPUTFILE=time_redist_xxf_to_yxf_level_1.in
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=time_redist_xxf_to_yxf_level_2.in
endif

all: time_redist_xxf_to_yxf
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_xxf_to_yxf,$(INPUTFILE)) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_redist_xxf_to_yxf.o:  ${TEST_MODS} 

time_redist_xxf_to_yxf: time_redist_xxf_to_yxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

