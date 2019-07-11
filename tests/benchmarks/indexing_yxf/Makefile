
sinclude ../Makefile.tests_common
all: time_indexing_yxf
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_yxf,time_indexing_yxf_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_yxf,time_indexing_yxf_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_indexing_yxf,time_indexing_yxf_yxles.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_yxf.o:  ${TEST_MODS} 

time_indexing_yxf: time_indexing_yxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

