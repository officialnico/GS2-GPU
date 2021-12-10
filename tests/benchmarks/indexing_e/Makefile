
sinclude ../Makefile.tests_common
all: time_indexing_e
	  echo "NOTE: indexing_e benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_e,time_indexing_e_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_e,time_indexing_e_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_e,time_indexing_e_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_e.o:  ${TEST_MODS} 

time_indexing_e: time_indexing_e.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

