
sinclude ../Makefile.tests_common
all: time_fields_local
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_fields_local,time_fields_local.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_fields_local,time_fields_local.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_fields_local,time_fields_local.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_fields_local.o:  ${TEST_MODS} 

time_fields_local: time_fields_local.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

