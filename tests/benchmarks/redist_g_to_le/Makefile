
sinclude ../Makefile.tests_common

INPUTFILE_BASE=time_redist_g_to_le_level_
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=${INPUTFILE_BASE}2
else
INPUTFILE=${INPUTFILE_BASE}1
endif

all: time_redist_g_to_le
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_le,$(INPUTFILE)_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_le,$(INPUTFILE)_yxles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_le,$(INPUTFILE)_lexys.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_redist_g_to_le.o:  ${TEST_MODS}

time_redist_g_to_le: time_redist_g_to_le.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

