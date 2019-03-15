#
#
#  Makefile for GS2 Gyrokinetic Turbulence code 
#
#  (requires GNU's gmake)
#
GK_PROJECT := gs2
#
#  Makefile written by Bill Dorland and Ryusuke Numata
#  Modified by GS2 contributors
#
#
# * Changelogs
# 26/09/10: added targets doc_depend and doc to make doxygen documentation. EGH.
#	12/17/09: make it explicit that the codes are written in
#		  Fortran 95 Standard
#	10/12/09: share the Makefile with gsp gyrokinetic PIC code
#	04/11/09: * share the Makefile with rmhdper reduced MHD code
#		  * include progject specific target definitions
#		    Makefile.target_$(GK_PROJECT)
#	04/06/09: SYSTEM environment variable is replaced by GK_SYSTEM
#	02/25/09: gs2 replaced by trinity (MAB)
#	01/26/09: USE_C_INDEX is imported to gs2 by TT
#       12/11/08: add support for NAGWare and Lahey compilers
#       10/28/08: some non-standard macros respect environment variables
#       10/27/08: default commands for MPI and HDF5 are defined
#       09/26/08: hack for treating intermediate f90 files properly
#	 	 (this is to keep all intermediate f90 files when 
#		 .PRECIOUS is given)
#	07/01/08: new switches USE_LOCAL_SPFUNC and USE_NAGLIB
#	06/16/08: clean up unused statements related to PLATFORM_LINKS
#
# * Available Compilers (tested on limited hosts)
#   (must be Fortran 95 Standard compliant)
#
# Intel ifort
# GNU's gfortran and g95
# IBM XL Fortran xlf90
# PathScale Compiler Suite pathf90
# The Portland Group pgf90
# NAGWare f95 (v5.1)
# Lahey/Fujitsu Fortran lf95
# 
# * Frequently Tested Hosts, Systems
#
# Standard Linux
# Standard Mac OS X with MacPorts
# Franklin at NERSC and Jaguar at NCCS (Cray XT4 with PGI)
# Bassi at NERSC (IBM Power 5 with IBM XL Fortran)
# Ranger (Sun Constellation Linux Cluster with Intel)
#
# * Switches:
#
# Here is a list of switches with simple explanation.
# In the brackets, values accepted are shown,
# where "undefined" means blank.
# Switches with (bin) are checked if they are defined or not
# What values they have do not matter.
# Be careful that DEBUG=off means DEBUG=on.
#
# Fortran Standard Specification (95, 2003, 2008, undefined)
FORTRAN_SPEC ?=
# turns on debug mode (bin)
DEBUG ?=
# turns on scalasca instrumentation mode (bin)
SCAL ?=
# turns on profile mode (gprof,ipm)
PROF ?=
# optimization (on,aggressive,undefined)
OPT ?= on
# prevents linking with shared libraries (bin)
STATIC ?=
# promotes precisions of real and complex (bin)
DBLE ?= on
# promotes precisions of double and double complex (bin)
QUAD ?=
# turns on distributed memory parallelization using MPI (mpi2, mpi3)
USE_MPI ?= mpi2
# turns on SHMEM parallel communications on SGI (bin)
USE_SHMEM ?=
# which FFT library to use (fftw,fftw3,mkl_fftw,undefined) 
USE_FFT ?= fftw3
# uses netcdf library (bin)
USE_NETCDF ?= on
# uses parallel netcdf library
USE_PARALLEL_NETCDF ?= 
# uses hdf5 library (bin)
USE_HDF5 ?=  
# Use function pointer in layouts_indices.c (bin)
# see also README
USE_C_INDEX ?= 
# Use local random number generator (mt,undefined)
# see also README
USE_LOCAL_RAN ?=
# Use posix for command_line (bin)
USE_POSIX ?=
# Use Fortran 2003, 2008 Intrinsics (bin)
# (Note this overwhelms other special function options)
USE_F200X_INTRINSICS ?=
# Use local special functions (bin)
USE_LOCAL_SPFUNC ?= 
# Use nag libraray (spfunc,undefined)
USE_NAGLIB ?= 
# Make GS2 into a library which can be called by external programs
MAKE_LIB ?=
# Include higher-order terms in GK equation arising from low-flow physics
LOWFLOW ?=
# Compile with PETSC/SLEPC support for determining eigensystem (bin).
# TODO : Rename USE_EIG for consistency
WITH_EIG ?=
# Include LAPACK and GSL libraries for genquad scheme
GENQUAD ?=
# Compile with gcov flags to check for code coverage
USE_GCOV ?=
HAS_GCOV=$(USE_GCOV)

# Compile with the new simplified diagnostics module (bin)
USE_NEW_DIAG ?=on
# Does the compiler version support the iso_c_binding features of Fortran?
HAS_ISO_C_BINDING ?= on
# If on, build with position independent code, needed if GS2 is going to be included
# in a shared library or linked dynamically (e.g. it's needed by GRYFX).
# This *may* have performance implications, but it hasn't been tested yet.
USE_FPIC ?= 
# Defines the size of the constant character arrays that are used to store the run
# name. 2000 has been selected as the default as being longer than most people will 
# need. However, this can be overriden by setting this parameter either in the environment
# or when executing make. Make sure you rebuild constants.fpp if you want to change this.
# e.g. 
#       $ touch utils/constants.fpp
#       $ ./build_gs2 -c 'RUN_NAME_SIZE=<newsize>'
RUN_NAME_SIZE?=2000
# If defined, fpp files will be preprocessed and compiled in one step, as
# is normal for modern FORTRAN compilers (bin)
ONE_STEP_PP?=
# If defined then don't show the full command being run everywhere
QUIET           ?= yes
# If defined then use escape commands to colourise some output
# Note hack to support both spellings
WITH_COLOUR     ?= $(WITH_COLOR)
#
# * Targets:
#
#  depend: generate dependency
#  test_make: print values of variables
#  clean: clean up
#  distclean: does "make clean" + removes platform links & executables
#  tar: pack
#
###############################################################
# A default empty rule to say we can't rebuild this makefile so don't bother trying
Makefile: ;
############################################################### DEFAULT VALUE
#
# These variables can be set in platform-dependent Makefile.
#
SHELL           = bash
MAKE		= make
CPP		= cpp
CPPFLAGS        = 
FC		= f90
MPIFC		?= mpif90
H5FC		?= h5fc
H5FC_par	?= h5pfc
F90FLAGS	=
F90OPTFLAGS	= 
CC		= cc
MPICC		?= mpicc
H5CC		?= h5cc
H5CC_par	?= h5pcc
CFLAGS		= 
COPTFLAGS 	=
LD 		= $(FC)
LDFLAGS 	= $(F90FLAGS)
ARCH 		= ar
ARCHFLAGS 	= cr
RANLIB		= ranlib
AWK 		= awk
PERL		= perl
FORD           ?= ford
MKDIR           = mkdir -p

MPI_INC	?=
MPI_LIB ?=
FFT_INC ?=
FFT_LIB ?=
EIG_INC ?=
EIG_LIB ?=
GENQUAD_INC ?=
GENQUAD_LIB ?=
NETCDF_INC ?=
NETCDF_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
MACRO_DEFS ?=

#Patch file
PATCHFILE:= .patch.make

################################################### SET COMPILE MODE SWITCHES

# It makes no sense to set USE_PARALLEL_NETCDF without
# setting USE_HDF5... so we set USE_HDF5 here. 
# If you want to use parallel-netcdf instead of netcdf/hdf5 you can 
# override this behaviour in the system makefile
ifdef USE_PARALLEL_NETCDF
ifndef USE_HDF5
	USE_HDF5=on
endif
endif

######################################################### PLATFORM DEPENDENCE

# compile mode switches (DEBUG, TEST, PROF, OPT, STATIC, DBLE)
# must be set before loading Makefile.$(GK_SYSTEM) because they may affect
# compiler options.
# However, Makefile.local may override some options set in Makefile.$(GK_SYSTEM),
# thus it is included before and after Makefile.$(GK_SYSTEM)
Makefile.local: ;
sinclude Makefile.local

# include system-dependent make variables
ifndef GK_SYSTEM
	ifdef SYSTEM
$(warning SYSTEM environment variable is obsolete)
$(warning use GK_SYSTEM instead)
	GK_SYSTEM = $(SYSTEM)
	else
$(error GK_SYSTEM is not set)
	endif
endif


include Makefile.$(GK_SYSTEM)

# Files that we know are going to be in the submodule directories.
# This is in case the directories exist but the submodules haven't
# actually been initialised
MAKEFILES_SENTINEL = Makefiles/Makefile.bitbucket
UTILS_SENTINEL = utils/mp.fpp

.PRECIOUS: $(MAKEFILES_SENTINEL) $(UTILS_SENTINEL)

# Make sure the submodules exist; they all currently share the same recipe
$(UTILS_SENTINEL) : $(MAKEFILES_SENTINEL)

$(MAKEFILES_SENTINEL) submodules:
	@echo "Downloading submodules"
	git submodule update --init --recursive


.Makefile.makefile_exists_helper: ;

# Here we have a dependency on $(MAKEFILES_SENTINEL) such that this
# will fetch the submodules first. We then try to build a null target
# in a helper makefile which includes Makefile.$(GK_SYSTEM). This helper
# makefile has a recipe to built Makefile.$(GK_SYSTEM) which just raises
# an error -- this allows us to abort if the makefile can't be found even
# after fetching the submodules
Makefile.$(GK_SYSTEM): $(MAKEFILES_SENTINEL)
	@$(MAKE) -f .Makefile.makefile_exists_helper --no-print-directory

# include Makefile.local if exists
sinclude Makefile.local

#Include utility makefile
Makefile.utils: ;
include Makefile.utils

#############################################################################

#Record the top level path
GK_HEAD_DIR:=$(PWD)

UTILS=utils
GEO=geo

ifndef ONE_STEP_PP
CPPFLAGS	+= -P -traditional
endif

GIT_HASH=$(shell git rev-parse HEAD)
MACRO_DEFS+=-DGIT_HASH='"$(GIT_HASH)"'

#The magic sed esacpes either ( or ) by replacing with \( or \)
GIT_BRANCH=$(shell git branch | grep \* | cut -d ' ' -f2- | sed 's|\([()]\)|\\\1|g')
MACRO_DEFS+=-DGIT_BRANCH='"$(GIT_BRANCH)"'

# Find if there are any modified tracked files (except Makefile.depend)
ifeq ($(shell git status --short -uno -- . | wc -l), 0)
	GIT_STATE="clean"
	MACRO_DEFS+=-DGIT_STATE='$(GIT_STATE)'
$(shell touch $(PATCHFILE))
else
	GIT_STATE="modified"
	MACRO_DEFS+=-DGIT_STATE='$(GIT_STATE)'

# Here we make a patch file, which could be read by the program if we wanted to keep a record
$(shell git diff -w -- .  > $(PATCHFILE))
endif
MACRO_DEFS+=-DGIT_PATCH_FILE='"$(PATCHFILE)"'


# Define RUN_NAME_SIZE
MACRO_DEFS+=-DRUN_NAME_SIZE=$(RUN_NAME_SIZE) 

############################################################################# Deal with user flags

# Define RELEASE based on the contents of the RELEASE
# file, which is empty unless we are in a release.
RELEASE_FILE = $(shell cat RELEASE)
ifneq ($(RELEASE_FILE),) # If the file is not empty
	RELEASE=$(RELEASE_FILE)
	MACRO_DEFS+=-DRELEASE=$(RELEASE_FILE)
	DOC_FOLDER=releases/$(RELEASE_FILE)/
endif

# Define GK_SYSTEM for runtime_tests so that we can get the system name
# at runtime
MACRO_DEFS+=-DGK_SYSTEM='"$(GK_SYSTEM)"'

ifdef USE_HDF5
	ifndef USE_MPI
$(error Currently, USE_HDF5 only works with USE_MPI)
	endif
endif


ifndef USE_NETCDF
ifdef USE_NEW_DIAG
$(error 'USE_NETCDF is off/undefined. The new diagnostics module requires netcdf. Please build without the new diagnostics module (i.e. set USE_NEW_DIAG=) ')
endif
endif

ifndef USE_FFT
$(warning USE_FFT is undefined/off)
$(warning Nonlinear runs are therefore not sensible)
endif

ifdef USE_MPI
	FC = $(MPIFC)
	CC = $(MPICC)
	MACRO_DEFS += -DMPI
endif

ifeq ($(USE_MPI),mpi3)
	MACRO_DEFS += -DMPI3
endif

ifdef HAS_ISO_C_BINDING
	MACRO_DEFS += -DISO_C_BINDING
endif

ifeq ($(USE_NEW_DIAG),on)
	MACRO_DEFS += -DNEW_DIAG
endif

ifeq ($(USE_FPIC),on)
	CFLAGS+=-fPIC
	F90FLAGS += -fPIC
endif

ifdef USE_SHMEM
	MACRO_DEFS += -DSHMEM
endif

ifeq ($(USE_FFT),fftw)
	MACRO_DEFS += -DFFT=_FFTW_
	FFT_LIB ?= -lfftw -lrfftw
endif

ifeq ($(USE_FFT),fftw3)
	MACRO_DEFS += -DFFT=_FFTW3_
        CPPFLAGS += $(FFT_INC)
	FFT_LIB ?= -lfftw -lrfftw
endif

ifdef GENQUAD
	GENQUAD_LIB ?= -llapack 
	MACRO_DEFS += -DGENQUAD
endif

ifeq ($(USE_FFT),mkl_fftw)
	MACRO_DEFS += -DFFT=_FFTW_
endif

ifdef USE_NETCDF
	NETCDF_LIB ?= -lnetcdf
	MACRO_DEFS += -DNETCDF
endif

ifdef USE_HDF5
	FC = $(H5FC_par)
	CC = $(H5CC_par)
	ifdef USE_PARALLEL_NETCDF
		MACRO_DEFS += -DNETCDF_PARALLEL
	endif
	MACRO_DEFS += -DHDF
endif

ifdef USE_C_INDEX
	MACRO_DEFS += -DUSE_C_INDEX
endif

ifeq ($(USE_LOCAL_RAN),mt)
	MACRO_DEFS += -DRANDOM=_RANMT_
endif

ifdef USE_POSIX
	MACRO_DEFS += -DPOSIX
endif

ifdef USE_F200X_INTRINSICS
	MACRO_DEFS += -DF200X_INTRINSICS -DSPFUNC=_SPF200X_
else
	ifdef USE_LOCAL_SPFUNC
		MACRO_DEFS += -DSPFUNC=_SPLOCAL_
	else
		ifeq ($(findstring spfunc,$(USE_NAGLIB)),spfunc)
			MACRO_DEFS += -DSPFUNC=_SPNAG_
		endif
	endif
endif

ifdef USE_NAGLIB
	ifeq ($(NAG_PREC),dble)
		ifndef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		MACRO_DEFS += -DNAG_PREC=_NAGDBLE_
	endif
	ifeq ($(NAG_PREC),sngl)
		ifdef DBLE
$(warning Precision mismatch with NAG libarray)	
		endif
		MACRO_DEFS += -DNAG_PREC=_NAGSNGL_
	endif
endif

ifndef DBLE
	MACRO_DEFS += -DSINGLE_PRECISION
endif

ifdef QUAD
	DBLE = on
	MACRO_DEFS += -DQUAD_PRECISION
endif

ifdef MAKE_LIB
	MACRO_DEFS += -DMAKE_LIB
endif

ifdef LOWFLOW
	MACRO_DEFS += -DLOWFLOW
endif

#Setup the flags for using the eigensolver
ifdef WITH_EIG
	EIG_INC += -I$(PETSC_DIR)/include -I$(SLEPC_DIR)/include
	ifdef PETSC_ARCH
		EIG_LIB += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib
		EIG_INC += -I$(PETSC_DIR)/$(PETSC_ARCH)/include
	else
		EIG_LIB += -L$(PETSC_DIR)/lib
	endif
	ifdef SLEPC_ARCH
		EIG_LIB += -L$(SLEPC_DIR)/$(SLEPC_ARCH)/lib
		EIG_INC += -I$(SLEPC_DIR)/$(SLEPC_ARCH)/include
	else
		EIG_LIB += -L$(SLEPC_DIR)/lib
	endif
	EIG_LIB += -lslepc -lpetsc

	MACRO_DEFS += -DWITH_EIG
	CPPFLAGS += $(EIG_INC)
	CFLAGS += -DWITH_EIG 
endif 

ifdef SCAL
   FC:= scalasca -instrument $(FC)
endif

#Make empty targets if not using the new diagnostics
ifeq ($(USE_NEW_DIAG),on)
Makefile.diagnostics/Makefile.diagnostics: ;
sinclude diagnostics/Makefile.diagnostics
else
distclean_simpledataio:
clean_simpledataio:
diagnostics:
simpledataio:
endif

############################################################################# Done dealing with user flags

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(NETCDF_LIB) $(HDF5_LIB) $(FFT_LIB) \
		   $(NAG_LIB) $(EIG_LIB) $(GENQUAD_LIB)
INCLUDES = $(DEFAULT_INC) $(MPI_INC) $(NETCDF_INC) $(HDF5_INC) $(FFT_INC)\
		   $(SIMPLEDATAIO_INC) $(EIG_INC) $(GENQUAD_INC) $(PFUNIT_INC)
F90FLAGS+= $(F90OPTFLAGS) $(INCLUDES)
CFLAGS += $(COPTFLAGS) \
	   $(DEFAULT_INC) $(MPI_INC) $(NETCDF_INC) $(HDF5_INC) $(FFT_INC)


DATE=$(shell date +%y%m%d)
TARDIR=$(GK_PROJECT)_$(DATE)
TOPDIR=$(CURDIR)
ifeq ($(notdir $(CURDIR)), $(UTILS))
	TOPDIR=$(subst /$(UTILS),,$(CURDIR))
endif
ifeq ($(notdir $(CURDIR)), $(GEO))
	TOPDIR=$(subst /$(GEO),,$(CURDIR))
endif
ifneq ($(TOPDIR),$(CURDIR))
	SUBDIR=true
endif

VPATH = $(GK_HEAD_DIR):$(UTILS):$(GEO):Aux:diagnostics:../$(UTILS):../$(GEO):../diagnostics$(SIMPLEDATAIO_VPATH)
# this just removes non-existing directory from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH = .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")
#
ifdef SUBDIR
	VPATH +=:..
endif

# most common include and library directories
DEFAULT_INC_LIST = . $(GK_HEAD_DIR)/$(UTILS) $(GK_HEAD_DIR)/$(GEO) .. ../$(UTILS) ../$(GEO)
DEFAULT_LIB_LIST =
# This default library path list can simplify the procedure of porting,
# however, I found this (actually -L/usr/lib flag) causes an error
# when linking gs2 at bassi (RS6000 with xl fortran)
DEFAULT_INC=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

# List of fortran files generated by ruby. 
# Note these are not deleted at the end of the compilation:
# in fact, they are stored in the repository 
# because some antique systems don't have ruby
FORTFROMRUBY=$(subst generate_,,$(patsubst %.rb,%.f90,$(wildcard *.rb */*.rb)))
FORTFROMRUBY+= gs2_init.f90 overrides.f90

# list of intermediate f90 files generated by preprocessor
F90FROMFPP = $(patsubst %.fpp,%.f90,$(wildcard *.fpp */*.fpp))
F90FROMFPP += tests/unit_tests/nonlinear_terms/test_nonlinear_terms.f90
F90FROMFPP += tests/unit_tests/gs2_diagnostics_new/test_gs2_diagnostics_new.f90
ifdef USE_NEW_DIAG
F90FROMFPP += diagnostics/simpledataio/src/simpledataiof.f90 diagnostics/simpledataio/src/simpledataio_write.f90 diagnostics/simpledataio/src/simpledataio_read.f90
endif

####################################################################### RULES

.SUFFIXES:

.PHONY: all $(GK_PROJECT)_all submodules
.PHONY: depend tar tar_exec test_make diagnostics simpledataio
.PHONY: help TAGS
.PHONY: clean cleanlib cleanconfig clean_coverage distclean unlink
.PHONY: doc docs/html docs/graphs cleandoc check_ford_install
.PHONY: doxygen_doc clean_doxygendoc sync_doc check_coderunner_install sync_input_doc

# comment this out to keep intermediate .f90 files
ifdef KEEP_PROCESSED
.PRECIOUS: $(F90FROMFPP)
endif

.INTERMEDIATE: $(F90FROMFPP)

####################################################################### Source code targets

# On many systems the two steps of preprocessing and compiling
# can be done simultaneously. This has many advantages, the
# most important of which is correctly reporting line numbers in 
# the .fpp when debugging. 
#
# On some systems (e.g. Mac OS X) using cpp doesn't work at all
# and you HAVE to preprocess and compile in one step.
ifdef ONE_STEP_PP
%.o: %.fpp
	$(call messageOp,"Compiling and preprocessing","\\t:",$<)
	$(QUIETSYM)$(FC) $(MACRO_DEFS) $(F90FLAGS) $(F90FLAGS_SFXJUNK) -c $<
%.o: %.F90
	$(call messageOp,"Compiling and preprocessing","\\t:",$<)
	$(QUIETSYM)$(FC) $(MACRO_DEFS) $(F90FLAGS) $(F90FLAGS_SFXJUNK) -c $<

else
# These are special rules for the suffix problem of absoft
# (not tested)
$(GK_PROJECT)_transforms.o: $(GK_PROJECT)_transforms.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<
$(GK_PROJECT)_io.o: $(GK_PROJECT)_io.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
$(GK_PROJECT)_save.o: $(GK_PROJECT)_save.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFX2) -c $<
mp.o: mp.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFX1) -c $<
fft_work.o: fft_work.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFX0) -c $<
endif

%.o : %.pf
	$(PFUNIT_DIR)/bin/pFUnitParser.py $< $(basename $<).F90
	$(QUIETSYM)$(FC) $(F90FLAGS) -c $(basename $<).F90 $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@

%.f90: %.rb
	$(call messageOp,"Generating from ruby","\\t:",$<)
	$(QUIETSYM)$(call RUBY_GENERATE,$@,$<)

%.o: %.f90
	$(call messageOp,"Compiling","\\t:",$<)
	$(QUIETSYM)$(FC) $(F90FLAGS) $(F90FLAGS_SFXJUNK) -c $<

# The depend command needs the preprocessed source files. 
# This is a little clunky, and maybe could be fixed in the future.
%.f90: %.fpp
	$(call messageOp,"Preprocessing","\\t:",$<)
	$(QUIETSYM)$(CPP) $(CPPFLAGS) $(MACRO_DEFS) $< -o $@
%.f90: %.F90
	$(call messageOp,"Preprocessing","\\t:",$<)
	$(QUIETSYM)$(CPP) $(CPPFLAGS) $(MACRO_DEFS) $< -o $@

%.o: %.c
	$(call messageOp,"Compiling C","\\t:",$<)
	$(QUIETSYM)$(CC) $(CFLAGS) $(subst -traditional,,$(subst -C,,$(subst -P,,$(CPPFLAGS)))) -c $<

# This prevents error messages like m2c: Command not found
# This is from a implicit makefile rule that we're overriding
%.o : %.mod


DEPEND=Makefile.depend
DEPEND_CMD=$(PERL) scripts/fortdep

# Only define recipe to build dependency file if it doesn't exist
ifeq ("$(wildcard $(DEPEND))", "")
.PRECIOUS: $(DEPEND)
$(DEPEND): $(UTILS_SENTINEL)
$(DEPEND): $(UTILS_SENTINEL) depend
#Provide a no-op recipe
	@echo -n ""
endif

# If this is true then we keep the preprocessed files generated in the process
# of producing the dependency file. This is a minor optimisation to make
KEEP_PROCESSED_FILES_FROM_DEPEND?="no"
ifdef ONE_STEP_PP
KEEP_PROCESSED_FILES_FROM_DEPEND="no"
endif

depend: $(FORTFROMRUBY) $(UTILS_SENTINEL) Makefile.$(GK_SYSTEM)
	$(call messageOp,"Generating dependency file","\\t:",$(DEPEND))
# Note we force UTILS to be included in the paths that we pass to fortdep
# We touch the dependency file first as fortdep will call make in order to try
# and preprocess files, so we need this to avoid recursion
	$(QUIETSYM)touch $(DEPEND) && $(DEPEND_CMD) -m "$(MAKE) --no-print-directory" -1 -o -v=1\
		--keep=$(KEEP_PROCESSED_FILES_FROM_DEPEND) $(VPATH):$(UTILS) ||\
		rm $(DEPEND)
	$(call messageOp,"Dependency file created at","\\t:",$(DEPEND))
##################################################################### Main TARGETS
# .DEFAULT_GOAL works for GNU make 3.81 (or higher)
# For 3.80 or less, see all target
.DEFAULT_GOAL := $(GK_PROJECT)_all
ifeq ($(notdir $(CURDIR)),utils)
	.DEFAULT_GOAL := utils_all
endif
ifeq ($(notdir $(CURDIR)),geo)
	.DEFAULT_GOAL := geo_all
endif

all: $(.DEFAULT_GOAL)

#TODO : Add a rule to automatically build the depend file when appropriate.
ifneq ("$(MAKECMDGOALS)", "distclean")
ifneq ("$(MAKECMDGOALS)", "depend")
include $(DEPEND)
endif
endif

#If we're using fft(w) then we want to link the fft_save_wisdom object file
#this currently has to be handled specially as it comes from a c source file
#rather than fortran (so isn't picked up by fortdep).
ifdef USE_FFT
    utils_mod += fft_save_wisdom.o
    gs2_mod += fft_save_wisdom.o
    ingen_mod += fft_save_wisdom.o
    generate_fftw_wisdom_mod += fft_save_wisdom.o
endif

ifdef USE_C_INDEX
    gs2_mod += layouts_indices.o
endif

layouts_indices.o: layouts_type.h
layouts_type.h: layouts_type.f90
	$(AWK) -f makehead.awk $^ > $@

# Make sure template dependencies are specified
gs2_init.f90: templates/gs2_init_template.f90

geometry_c_interface.o: geometry_c_interface.h

gryfx_libs: utils.a geo.a geo/geometry_c_interface.o

############################################################# Get project specific targets
Makefile.target_$(GK_PROJECT): ;
include Makefile.target_$(GK_PROJECT)
############################################################# System config targets
ifdef STANDARD_SYSTEM_CONFIGURATION
system_config:
	@echo "#!/bin/bash " > system_config
	@echo "$(STANDARD_SYSTEM_CONFIGURATION)" >> system_config

else
.PHONY: system_config
system_config:
	$(error "STANDARD_SYSTEM_CONFIGURATION is not defined for this system")
endif

############################################################# Cleanup targets
clean: clean_simpledataio clean_coverage
	$(call messageRm,"Object and module files")
	-$(QUIETSYM)rm -f *.o *.mod *.g90 *.h core */core

CLEANCOMMAND=echo $$$$PWD
CLEANCOMMAND=rm -f *.o *.error *.out *.out.nc *.gridgen.200 *.lpc *.vres *.fields *.g *.fort.?? *.mod .*.scratch *.timing.* *.moments *.cdf *.jext *.parity *.heat *.heat2 *.vres2 *.amoments *.eigenfunc *.nc* *.mom2 *.epar .*.in results_of_test.txt *.stop *.fftw_wisdom *.phase *.kpar *.interp *.dist 

ifdef CLEAN_TEXTFILES
	CLEANCOMMAND+= *~ *.orig
endif

export CLEANCOMMAND

cleanlib:
	$(call messageRm,"Libraries")
	-$(QUIETSYM)rm -f *.a

cleanconfig:
	$(call messageRm,"System config")
	-$(QUIETSYM)rm -f system_config .tmp_output

clean_coverage:
ifeq ($(HAS_GCOV),on)
	$(call messageRm,"Gcov files")	
	-$(QUIETSYM)rm -f *.gcov *.gcda *.gcno
endif

clean_depend:
	$(call messageRm,"Dependency file")	
	-$(QUIETSYM)rm -f $(DEPEND)

distclean: unlink clean cleanlib clean_tests clean_benchmarks distclean_simpledataio cleanconfig clean_doxygendoc cleandoc clean_coverage clean_depend pfunit_clean

unlink:
	$(call messageRm,"Preprocessed files")	
	-$(QUIETSYM)rm -f $(F90FROMFPP) layouts_type.h

############################################################# Tar targets
tar:
	@[ ! -d $(TARDIR) ] || echo "ERROR: directory $(TARDIR) exists. Stop."
	@[ -d $(TARDIR) ] || $(MAKE) tar_exec

$(TARDIR):
	$(QUIETSYM)$(MKDIR) $@

### setting tar_exec local $(TARLIST*) variables
# expand wildcards listed $(TARLIST_wild) in ( $(TARLIST_dir) + . )
# directories and add them into TARLIST
tar_exec: TARLIST = makehead.awk scripts/fortdep AstroGK.in
tar_exec: TARLIST_dir = Makefiles utils geo Aux
tar_exec: TARLIST_wild = *.f90 *.fpp *.inc *.c Makefile Makefile.* README README.*
tar_exec: TARLIST += $(foreach dir,. $(TARLIST_dir),$(wildcard $(addprefix $(dir)/,$(TARLIST_wild))))
tar_exec: | $(TARDIR)

tar_exec:
	@for dir in $(TARLIST_dir) ;\
	  do ( [ ! -d $$dir ] ||  $(MKDIR) $(TARDIR)/$$dir ; ) ;\
	done
	@for name in $(TARLIST) ;\
	  do ( [ -f $$name ] && ln $$name $(TARDIR)/$$name ; ) ;\
	done
	@tar cvf - $(TARDIR) | bzip2 -9 > $(TARDIR).tar.bz2
	@rm -rf $(TARDIR) # This looks dangerous

############################################################# Testing targets
#This is the location of the individual test suites
TEST_DIR:=tests
$(TEST_DIR)/Makefile.tests_and_benchmarks: ;
sinclude $(TEST_DIR)/Makefile.tests_and_benchmarks

#This is the location of the pfunit test suites
PFUNIT_TEST_DIR:=$(TEST_DIR)/pfunit_tests
sinclude $(PFUNIT_TEST_DIR)/Makefile

############################################################# Introspection targets
test_make:
	$(call messageInfo,"Global settings:")
	$(call messageKeyVal, "GK_SYSTEM", $(GK_SYSTEM))
	$(call messageKeyVal, "VPATH", $(VPATH))
	$(call messageKeyVal, "CURDIR", $(CURDIR))
	$(call messageKeyVal, "TOPDIR", $(TOPDIR))
	$(call messageKeyVal, "MAKECMDGOALS", $(MAKECMDGOALS))
	$(call messageKeyVal, ".DEFAULT_GOAL", $(.DEFAULT_GOAL))
	$(call messageKeyVal, "DEFAULT_LIB", $(DEFAULT_LIB))
	$(call messageKeyVal, "DEPEND_CMD", $(DEPEND_CMD))
	$(call messageKeyVal, "DEPEND", $(DEPEND))
	$(call messageKeyVal, "QUIET", $(QUIET))
	$(call messageKeyVal, "WITH_COLOUR", $(WITH_COLOUR))
	$(call messageKeyVal, "GIT_HASH", $(GIT_HASH))
	$(call messageKeyVal, "GIT_STATE", $(GIT_STATE))
	$(call messageKeyVal, "GIT_BRANCH", $(GIT_BRANCH))
	$(call messageKeyVal, "PATCHFILE", $(PATHFILE))
	$(call messageKeyVal, "F90FROMFPP", $(F90FROMFPP))
	@echo
	$(call messageInfo,"Commands:")
	$(call messageKeyVal, "COMPILER", $(COMPILER))
	$(call messageKeyVal, "FC", $(FC))
	$(call messageKeyVal, "CC", $(CC))
	$(call messageKeyVal, "LD", $(LD))
	$(call messageKeyVal, "CPP", $(CPP))
	$(call messageKeyVal, "ARCH", $(ARCH))
	$(call messageKeyVal, "RANLIB", $(RANLIB))
	$(call messageKeyVal, "PERL", $(PERL))
	$(call messageKeyVal, "SHELL", $(SHELL))
	$(call messageKeyVal, "MAKE", $(MAKE))
	$(call messageKeyVal, "FORD", $(FORD))
	$(call messageKeyVal, "MKDIR", $(MKDIR))
	@echo
	$(call messageInfo,"Compilation mode:")
	$(call messageKeyVal, "DEBUG", $(DEBUG))
	$(call messageKeyVal, "SCAL", $(SCAL))
	$(call messageKeyVal, "TEST", $(TEST))
	$(call messageKeyVal, "PROF", $(PROF))
	$(call messageKeyVal, "OPT", $(OPT))
	$(call messageKeyVal, "STATIC", $(STATIC))
	$(call messageKeyVal, "DBLE", $(DBLE))
	@echo
	$(call messageInfo,"Functionality:")
	$(call messageKeyVal, "USE_MPI", $(USE_MPI))
	$(call messageKeyVal, "USE_SHMEM", $(USE_SHMEM))
	$(call messageKeyVal, "USE_FFT", $(USE_FFT))
	$(call messageKeyVal, "USE_NETCDF", $(USE_NETCDF))
	$(call messageKeyVal, "USE_HDF5", $(USE_HDF5))
	$(call messageKeyVal, "WITH_EIG", $(WITH_EIG))
	$(call messageKeyVal, "USE_C_INDEX", $(USE_C_INDEX))
	$(call messageKeyVal, "USE_POSIX", $(USE_POSIX))
	$(call messageKeyVal, "USE_LOCAL_RAN", $(USE_LOCAL_RN))
	$(call messageKeyVal, "USE_SPFUNC", $(USE_SPFUNC))
	$(call messageKeyVal, "USE_NAGLIB", $(USE_NAGLIB))
	$(call messageKeyVal, "USE_NEW_DIAG", $(USE_NEW_DIAG))
	$(call messageKeyVal, "USE_GCOV", $(USE_GCOV))
	@echo
	$(call messageInfo,"Flags:")
	$(call messageKeyVal, "F90FLAGS", $(F90FLAGS))
	$(call messageKeyVal, "F90OPTFLAGS", $(F90OPTFLAGS))
	$(call messageKeyVal, "CFLAGS", $(CFLAGS))
	$(call messageKeyVal, "COPTFLAGS", $(COPTFLAGS))
	$(call messageKeyVal, "LDFLAGS", $(LDFLAGS))
	$(call messageKeyVal, "CPPFLAGS", $(CPPFLAGS))
	$(call messageKeyVal, "ARCHFLAGS", $(ARCHFLAGS))
	$(call messageKeyVal, "MACRO_DEFS", $(MACRO_DEFS))
	$(call messageKeyVal, "LD_LIBRARY_PATH", $(LD_LIBRARY_PATH))
	$(call messageKeyVal, "LIBS", $(LIBS))
	$(call messageKeyVal, "INCS", $(INCS))
ifdef USE_MPI
	$(call messageInfo,"MPI:")
	$(call messageKeyVal, "MPI_MOD_NAME", $(MPI_MOD_NAME))
	$(call messageKeyVal, "MPI_LIB", $(MPI_LIB))
	$(call messageKeyVal, "MPI_INC", $(MPI_INC))
else
	$(call messageInfo,"MPI: DISABLED")
endif
	@echo
ifdef USE_NETCDF
	$(call messageInfo,"Netcdf:")
	$(call messageKeyVal, "NETCDF_LIB", $(NETCDF_LIB))
	$(call messageKeyVal, "NETCDF_INC", $(NETCDF_INC))
else
	$(call messageInfo,"Netcdf: DISABLED")
endif
	@echo
ifdef USE_HDF5
	$(call messageInfo,"HDF5:")
	$(call messageKeyVal, "HDF5_LIB", $(HDF5_LIB))
	$(call messageKeyVal, "HDF5_INC", $(HDF5_INC))
else
	$(call messageInfo,"HDF5: DISABLED")
endif
	@echo
ifdef USE_FFT
	$(call messageInfo,"FFT:")
	$(call messageKeyVal, "FFT_LIB", $(FFT_LIB))
	$(call messageKeyVal, "FFT_INC", $(FFT_INC))
	$(call messageKeyVal, "FFT_DIR", $(FFT_DIR))
else
	$(call messageInfo,"FFT: DISABLED")
endif
	@echo
ifdef USE_NEW_DIAG
	$(call messageInfo,"SIMPLEDATAIO:")
	$(call messageKeyVal, "SIMPLEDATAIO_LIB", $(SIMPLEDATAIO_LIB))
	$(call messageKeyVal, "SIMPLEDATAIO_LIB_ABS", $(SIMPLEDATAIO_LIB_ABS))
	$(call messageKeyVal, "SIMPLEDATAIO_INC" $(SIMPLEDATAIO_INC))
else
	$(call messageInfo,"SIMPLEDATAIO: DISABLED")
endif
	@echo
ifdef WITH_EIG
	$(call messageInfo,"Eigensolver:")
	$(call messageKeyVal, "EIG_LIB", $(EIG_LIB))
	$(call messageKeyVal, "EIG_INC", $(EIG_INC))
else
	$(call messageInfo,"Eigensolver: DISABLED")
endif
	@echo

revision:
	@LANG=C git rev-parse HEAD > Revision

#This only works if a helplocal rule has been defined in the specific system makefile.
#It would be nice to provide a default helplocal target that just says no help available
help: helplocal

TAGS:	*.f90 *.fpp */*.f90 */*.fpp
	etags $^

############################################################# Documentation

ifneq ("$(wildcard $(shell which $(FORD) 2>/dev/null))","")
check_ford_install:
	$(call messageInfo,"Using ford at $(shell which $(FORD))")
else
check_ford_install:
	$(call messageInfo,"Ford command $(FORD) not in path -- is it installed?\\n\\tConsider installing with 'pip install --user ford' and add ${HOME}/.local/bin to PATH") ; which $(FORD)
endif

docs/html:
doc: docs/.docs.config.md check_ford_install
	$(FORD) $(subst -D,--macro ,$(MACRO_DEFS)) $(INCLUDES) -r $(GIT_HASH) $<

cleandoc: docs/html docs/graphs
	$(call messageRm,"FORD docs")	
	-$(QUIETSYM)rm -rf $^

doxygen_doc: $(F90FROMFPP) $(FORTFROMRUBY)
	doxygen doxygen_config 

clean_doxygendoc: #Doxygen documentation
	$(call messageRm,"Doxygen documentation")
	-$(QUIETSYM)rm -rf doc/

#TODO: This will need updating when documentation has migrated to bitbucket
sync_doc:
	$(call messageInfo,"Warning: Target sync_doc out of date."
	rsync -av --delete --exclude=releases doc/html/	${USER},gyrokinetics@web.sourceforge.net:htdocs/gs2_documentation/$(DOC_FOLDER)
	$(MKDIR) doc/html/releases
	echo "Options +Indexes" > doc/html/releases/.htaccess
	rsync -av doc/html/releases/ ${USER},gyrokinetics@web.sourceforge.net:htdocs/gs2_documentation/releases/


ifneq ("$(wildcard $(shell which $(coderunner) 2>/dev/null))","")
check_coderunner_install:
	$(call messageInfo,"Using coderunner at $(shell which coderunner)")
else
check_coderunner_install:
	$(call messageInfo,"Coderunner command coderunner not in path -- is it installed?") ; which coderunner
endif

#This has an external dependency on coderunner
sync_input_doc: |check_coderunner_install
	$(call messageInfo,"Warning: Target sync_input_doc out of date."
	coderunner cc synchronise_variables . -C gs2
	curl 'http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=GS2_Input_Parameters&action=edit' | sed 's/&amp;/\&/g' | sed 's/&quot;/"/g' | sed 's/&gt;/>/g' | sed 's/&lt;/</g'  | sed 's/&nbsp;/ /g' > wiki_in.tmp
	coderunner cc read_mediawiki_documentation wiki_in.tmp -C gs2
	coderunner cc write_mediawiki_documentation > wiki_out.txt -C gs2
