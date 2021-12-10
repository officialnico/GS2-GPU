
sinclude ../Makefile.tests_common

INPUTFILE_BASE=time_redist_g_to_gf_level_
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=${INPUTFILE_BASE}2
else
INPUTFILE=${INPUTFILE_BASE}1
endif

all: time_redist_g_to_gf
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_gf,$(INPUTFILE)_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_gf,$(INPUTFILE)_yxles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_g_to_gf,$(INPUTFILE)_lexys.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

time_redist_g_to_gf.o:  ${TEST_MODS} 

time_redist_g_to_gf: time_redist_g_to_gf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

