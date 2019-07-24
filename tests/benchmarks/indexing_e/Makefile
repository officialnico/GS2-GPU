
sinclude ../Makefile.tests_common
all: time_indexing_e
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_e,time_indexing_e_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_e,time_indexing_e_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_e,time_indexing_e_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_e.o:  ${TEST_MODS} 

time_indexing_e: time_indexing_e.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

