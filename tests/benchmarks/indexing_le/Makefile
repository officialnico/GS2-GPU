
sinclude ../Makefile.tests_common
all: time_indexing_le
	  echo "NOTE: indexing_xxf benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_le,time_indexing_le_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_le,time_indexing_le_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_le,time_indexing_le_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_le.o:  ${TEST_MODS} 

time_indexing_le: time_indexing_le.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

