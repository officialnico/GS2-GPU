
sinclude ../Makefile.tests_common

EXENAME=time_init_advance
INPUTFILE_BASE=time_init_advance_level_

INPUTFILE=${INPUTFILE_BASE}1
ifeq "$(BENCHMARK_LEVEL)" "2"
INPUTFILE=${INPUTFILE_BASE}2
endif
ifeq "$(BENCHMARK_LEVEL)" "3"
INPUTFILE=${INPUTFILE_BASE}3
endif

all: $(EXENAME)
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,$(EXENAME),$(INPUTFILE).in) ; \
	done

$(EXENAME).o:  ${TEST_MODS} 

$(EXENAME): $(EXENAME).o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
