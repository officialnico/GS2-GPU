
sinclude ../Makefile.tests_common
all: time_indexing_yxf
	  echo "NOTE: indexing_yxf benchmark is independent of processor count and always runs on one processor." ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_yxf,time_indexing_yxf_lexys.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_yxf,time_indexing_yxf_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,1,time_indexing_yxf,time_indexing_yxf_yxles.in) ; \

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_indexing_yxf.o:  ${TEST_MODS} 

time_indexing_yxf: time_indexing_yxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

