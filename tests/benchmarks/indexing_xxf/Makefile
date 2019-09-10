
sinclude ../Makefile.tests_common
all: time_indexing_xxf
	  echo "NOTE: indexing_xxf benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_xxf,time_indexing_xxf_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_xxf,time_indexing_xxf_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_xxf,time_indexing_xxf_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_xxf.o:  ${TEST_MODS} 

time_indexing_xxf: time_indexing_xxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

