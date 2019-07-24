
sinclude ../Makefile.tests_common
all: time_indexing
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing,time_indexing_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing,time_indexing_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing,time_indexing_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing.o:  ${TEST_MODS} 

time_indexing: time_indexing.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

