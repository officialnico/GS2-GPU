
sinclude ../Makefile.tests_common
all: time_indexing_lz
	  echo "NOTE: indexing_lz benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_lz,time_indexing_lz_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_lz,time_indexing_lz_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_lz,time_indexing_lz_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_lz.o:  ${TEST_MODS} 

time_indexing_lz: time_indexing_lz.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

