
sinclude ../Makefile.tests_common
all: time_redist_xxf_to_yxf
	for d in $(NTESTPROCS); \
	do \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_xxf_to_yxf,time_redists_xyles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_xxf_to_yxf,time_redists_yxles.in) ; \
	  $(call BENCHMARKCOMMAND,$$d,time_redist_xxf_to_yxf,time_redists_lexys.in) ; \
	done

TIMINGFILES=$(wildcard *.timing.*)

UPLOADCOMMAND = $(call UPLOAD_BENCHMARKS,$(TIMINGFILES))


# For building unit tests - to get rid of too much dependency
#TEST_MODS = ../../unit_tests.o ../../general_f0.o ../../species.o ../../file_utils.o ../../spl.o ../../command_line.o ../../mp.o ../../text_options.o ../../constants.o ../../analytical_falpha.o


time_redist_xxf_to_yxf.o:  ${TEST_MODS} 

time_redist_xxf_to_yxf: time_redist_xxf_to_yxf.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

