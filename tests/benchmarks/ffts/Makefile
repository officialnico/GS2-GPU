
sinclude ../Makefile.tests_common

EXENAME=time_ffts
INPUTFILE_BASE=time_ffts_level_
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=${INPUTFILE_BASE}2
else
INPUTFILE=${INPUTFILE_BASE}1
endif

all: $(EXENAME)
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,$(EXENAME),$(INPUTFILE)_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,$(EXENAME),$(INPUTFILE)_yxles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,$(EXENAME),$(INPUTFILE)_lexys.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))

$(EXENAME).o:  ${TEST_MODS} 

$(EXENAME): $(EXENAME).o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

