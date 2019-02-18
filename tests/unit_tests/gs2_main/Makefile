
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_main

all: test_gs2_main 
	$(call RUNCOMMAND,test_gs2_main,test_gs2_main.in)

test_gs2_main.o:  ${TEST_MODS} 

test_gs2_main: test_gs2_main.o $(TEST_MODS)
	$(call messageOp,"Linking ","\\t:",$^,"\\n\\t\-\-\-\-\-\-\>",$@)
	$(QUIETSYM)$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

