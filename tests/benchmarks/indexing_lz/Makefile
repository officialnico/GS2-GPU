
sinclude ../Makefile.tests_common
all: time_indexing_lz
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_lz,time_indexing_lz_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_lz,time_indexing_lz_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_lz,time_indexing_lz_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_lz.o:  ${TEST_MODS} 

time_indexing_lz: time_indexing_lz.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

