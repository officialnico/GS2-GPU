
sinclude ../Makefile.tests_common
all: time_indexing_le
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_le,time_indexing_le_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_le,time_indexing_le_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_le,time_indexing_le_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_le.o:  ${TEST_MODS} 

time_indexing_le: time_indexing_le.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

