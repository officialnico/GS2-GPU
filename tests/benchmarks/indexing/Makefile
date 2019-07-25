
sinclude ../Makefile.tests_common
all: time_indexing
	  echo "NOTE: indexing benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing,time_indexing_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing,time_indexing_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing,time_indexing_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing.o:  ${TEST_MODS} 

time_indexing: time_indexing.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

