
sinclude ../Makefile.tests_common
all: time_indexing_xxf
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_xxf,time_indexing_xxf_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_xxf,time_indexing_xxf_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_xxf,time_indexing_xxf_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_xxf.o:  ${TEST_MODS} 

time_indexing_xxf: time_indexing_xxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

