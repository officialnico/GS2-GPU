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
# uses lapack or compatible library (bin)
USE_LAPACK ?=
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
# Name of python exe
GS2_PYTHON_NAME ?= python3

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
# We default WITH_COLOR (and hence WITH_COLOUR) to on if
# MAKE_TERMOUT is set. This is an internal make variable which is set
# provided the output is going to a terminal.  In other words we don't
# default the colour on if we're redirecting output to file.
ifdef MAKE_TERMOUT
WITH_COLOR ?= yes
endif
WITH_COLOUR     ?= $(WITH_COLOR)
# If defined then show the full path in compilation statements etc.
SHOW_ABS_PATH ?= yes

###############################################################
# By default we export all variables defined in this Makefile
# or any included makefiles. This gives sub-make invocations
# (i.e. make calls from within our Makefile) access to the same
# variables. Really we should just default what is needed.
export
###############################################################

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
CPPFLAGS	=
FPPFLAGS	= -P -traditional
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
FORD           ?= $(GK_HEAD_DIR)/scripts/ford_interceptor.py
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
LAPACK_INC ?=
LAPACK_LIB ?=
HDF5_INC ?=
HDF5_LIB ?=
NAG_LIB ?=
NAG_PREC ?= dble
MACRO_DEFS ?=

# Makefile.local may override some options set in Makefile.$(GK_SYSTEM),
# thus it is included before and after Makefile.$(GK_SYSTEM)
Makefile.local: ;
sinclude Makefile.local

#### directories

#Record the top level path. Note we don't just use $(PWD) as this
#resolves to the directory from which make was invoked. The approach
#taken here ensures that GK_HEAD_DIR is the location of this Makefile.
#Note the realpath call removes the trailing slash so later we need
#to add a slash if we want to address subdirectories
GK_THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
GK_HEAD_DIR := $(realpath $(dir $(GK_THIS_MAKEFILE)))

#Directories to store outputs of different types
OBJDIR:= $(GK_HEAD_DIR)/objs
MODDIR:= $(GK_HEAD_DIR)/include
BINDIR:= $(GK_HEAD_DIR)/bin
LIBDIR:= $(GK_HEAD_DIR)/lib
GENERATED_SRC_DIR:= $(GK_HEAD_DIR)/generated_src

#Src files exist here
SRCDIR:= $(GK_HEAD_DIR)/src

#Programs exist here
PROGDIR:= $(SRCDIR)/programs

#Extra subdirectories which include code that might be needed
SUBDIRS= geo
ifdef USE_NEW_DIAG
SUBDIRS += diagnostics diagnostics/simpledataio/src
endif
SUBDIRS:= $(addprefix $(SRCDIR)/, $(SUBDIRS))

#External subdirectories which include code that might be needed
EXTERNALDIR:= $(GK_HEAD_DIR)/externals
EXTERNALSRCDIRS:= utils
EXTERNALSRCDIRS:= $(addprefix $(EXTERNALDIR)/, $(EXTERNALSRCDIRS))

SUBDIRS += $(EXTERNALSRCDIRS)

GEO:=$(SRCDIR)/geo
UTILS:=$(EXTERNALDIR)/utils

#Tests
TEST_DIR:=$(GK_HEAD_DIR)/tests
TEST_RESULT_DIR:= $(GK_HEAD_DIR)/test-results
COVERAGE_DIR:= $(GK_HEAD_DIR)/coverage

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
MAKEFILES_SENTINEL := Makefiles/Makefile.bitbucket
UTILS_SENTINEL := $(UTILS)/mp.fpp

.PRECIOUS: $(MAKEFILES_SENTINEL) $(UTILS_SENTINEL)

# Make sure the submodules exist; they all currently share the same recipe
$(UTILS_SENTINEL) $(MAKEFILES_SENTINEL) submodules:
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

# Here we add dependency on the utils_sentinel for a somewhat arbitrary
# file in utils that we will try to include. If $(UTILS)/define.inc is not
# present this dependency will ensure the submodules recipe is run.
$(UTILS)/define.inc: $(UTILS_SENTINEL)

# Here we try to include this file from utils. We don't want (or get) any
# information from including it, it is merely to trigger the submodules
# recipe when utils hasn't been checked out
include $(UTILS)/define.inc

# include Makefile.local if exists
sinclude Makefile.local

#Include utility makefile
Makefile.utils: ;
include Makefile.utils

#Indicate that we don't want to rebuild the COMPILER makefile
Makefiles/Makefile.$(COMPILER): ;
#############################################################################

GIT_HASH:=$(shell git rev-parse HEAD)
MACRO_DEFS+=-DGIT_HASH='"$(GIT_HASH)"'

GIT_HASH_MAKEFILES:=$(shell git submodule status Makefiles | cut -d ' ' -f2)
MACRO_DEFS+=-DGIT_HASH_MAKEFILES='"$(GIT_HASH_MAKEFILES)"'

GIT_HASH_UTILS:=$(shell git submodule status $(UTILS) | cut -d ' ' -f2)
MACRO_DEFS+=-DGIT_HASH_UTILS='"$(GIT_HASH_UTILS)"'

# Version string, including current commit if not on a release, plus
# '-dirty' if there are uncommitted changes
GIT_VERSION := $(shell git describe --tags --long --dirty)
MACRO_DEFS += -DGIT_VERSION='"$(GIT_VERSION)"'

#The magic sed esacpes either ( or ) by replacing with \( or \)
GIT_BRANCH:=$(shell git branch | grep \* | cut -d ' ' -f2- | sed 's|\([()]\)|\\\1|g')
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

######## PYTHON ###############################################################

PYTHON_EXE:=$(shell which $(GS2_PYTHON_NAME))
GS2_HAS_PYTHON:=$(if $(PYTHON_EXE),yes,no)
ifeq ($(GS2_HAS_PYTHON),yes)
PYTHON_VERSION:=$(shell $(PYTHON_EXE) --version 2>&1)
endif

###############################################################################

# Define RUN_NAME_SIZE
MACRO_DEFS+=-DRUN_NAME_SIZE=$(RUN_NAME_SIZE) 

############################################################################# Deal with user flags

# Define IS_RELEASE based on if current commit exactly matches
# a tagged version. If it does then we define RELEASE to reflect
# this tag. If it doesn't then we define RELEASE from the "closest"
# tag plus a modifier.
IS_RELEASE := $(shell git describe --tags --exact-match 2>/dev/null)
MACRO_DEFS+=-DIS_RELEASE=$(IS_RELEASE)

ifneq ($(IS_RELEASE),) # If we're on a commit corresponding to a tag
	RELEASE:=$(IS_RELEASE)
else
	RELEASE:=$(shell git describe --tags 2>/dev/null)
endif
MACRO_DEFS+=-DRELEASE='"$(RELEASE)"'

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

ifeq ($(findstring fftw3,$(USE_FFT)),fftw3)
	MACRO_DEFS += -DFFT=_FFTW3_
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

ifeq ($(USE_LOCAL_RAN),mt)
	MACRO_DEFS += -DRANDOM=_RANMT_
endif

ifdef USE_POSIX
	MACRO_DEFS += -DPOSIX
endif

ifdef USE_LAPACK
	MACRO_DEFS += -DLAPACK
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
		EIG_LIB += -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib
		EIG_INC += -I$(SLEPC_DIR)/$(PETSC_ARCH)/include
	else
		EIG_LIB += -L$(PETSC_DIR)/lib
	endif
#SLEPC_ARCH is non-standard and is unlikely to be set, it should
#really be PETSC_ARCH (as used above)
	ifdef SLEPC_ARCH
		EIG_LIB += -L$(SLEPC_DIR)/$(SLEPC_ARCH)/lib
		EIG_INC += -I$(SLEPC_DIR)/$(SLEPC_ARCH)/include
	else
		EIG_LIB += -L$(SLEPC_DIR)/lib
	endif
	EIG_LIB += -lslepc -lpetsc

	MACRO_DEFS += -DWITH_EIG
	CFLAGS += -DWITH_EIG 
endif 

ifdef SCAL
   FC:= scalasca -instrument $(FC)
endif

#Make empty targets if not using the new diagnostics
ifeq ($(USE_NEW_DIAG),on)
$(SRCDIR)/diagnostics/Makefile.diagnostics: ;
include $(SRCDIR)/diagnostics/Makefile.diagnostics
else
diagnostics:
endif

############################################################################# Done dealing with user flags

LIBS	+= $(DEFAULT_LIB) $(MPI_LIB) $(NETCDF_LIB) $(HDF5_LIB) $(FFT_LIB) \
		   $(NAG_LIB) $(EIG_LIB) $(GENQUAD_LIB) $(LAPACK_LIB)
INCLUDES = $(DEFAULT_INC) $(MPI_INC) $(NETCDF_INC) $(HDF5_INC) $(FFT_INC)\
		   $(EIG_INC) $(GENQUAD_INC) $(PFUNIT_INC) $(LAPACK_INC)

INCLUDES += -I$(MODDIR)

# Here we attempt to remove non-existent directories from INCLUDES. We assume that
# INCLUDES is currently just a list of entries of the form -I<some directory> so for
# each entry we strip the -I and then echo -I<some directory> if <some directory> is
# found to exist. This can help avoid warnings about non-existent directories
INCLUDES := $(foreach tmpvp,$(INCLUDES),$(shell [ -d $(echo $(tmpvp) | sed 's/^-I//g') ] && echo $(tmpvp)))

F90FLAGS+= $(F90OPTFLAGS) $(INCLUDES)
CFLAGS += $(COPTFLAGS) \
	   $(DEFAULT_INC) $(MPI_INC) $(NETCDF_INC) $(HDF5_INC) $(FFT_INC)

CPPFLAGS += $(MACRO_DEFS)
FPPFLAGS += $(MACRO_DEFS)

DATE:=$(shell date +%y%m%d)
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

VPATH = $(SUBDIRS) $(MODDIR) $(SRCDIR) $(PROGDIR) $(GENERATED_SRC_DIR)
# this just removes non-existing directory from VPATH
VPATH_tmp := $(foreach tmpvp,$(subst :, ,$(VPATH)),$(shell [ -d $(tmpvp) ] && echo $(tmpvp)))
VPATH := .:$(shell echo $(VPATH_tmp) | sed "s/ /:/g")

# most common include and library directories
DEFAULT_INC_LIST := . $(GK_HEAD_DIR)/$(UTILS) $(GK_HEAD_DIR)/$(GEO) .. ../$(UTILS) ../$(GEO)
DEFAULT_LIB_LIST :=
# This default library path list can simplify the procedure of porting,
# however, I found this (actually -L/usr/lib flag) causes an error
# when linking gs2 at bassi (RS6000 with xl fortran)
DEFAULT_INC:=$(foreach tmpinc,$(DEFAULT_INC_LIST),$(shell [ -d $(tmpinc) ] && echo -I$(tmpinc)))
DEFAULT_LIB:=$(foreach tmplib,$(DEFAULT_LIB_LIST),$(shell [ -d $(tmplib) ] && echo -L$(tmplib)))

# Commands for generating repetitive Fortran files from python scripts
ifeq ($(GS2_HAS_PYTHON),yes)

$(SRCDIR)/gs2_init_down.inc $(SRCDIR)/gs2_init_up.inc $(SRCDIR)/gs2_init_subroutines.inc $(SRCDIR)/gs2_init_level_list.inc: $(SRCDIR)/generate_gs2_init.py
	$(call messageOp,"Generating from python","\\t:",$(call processName,$@) ($(call processName,$^)))
	$(QUIETSYM)$(PYTHON_EXE) $< $@

$(GENERATED_SRC_DIR)/gs2_init.f90: $(SRCDIR)/gs2_init.fpp $(SRCDIR)/gs2_init_down.inc $(SRCDIR)/gs2_init_up.inc $(SRCDIR)/gs2_init_subroutines.inc $(SRCDIR)/gs2_init_level_list.inc

$(SRCDIR)/overrides.f90: $(SRCDIR)/generate_overrides.py
	$(call messageOp,"Generating from python","\\t:",$(call processName,$<))
	$(QUIETSYM)$(PYTHON_EXE) $< $@
endif

# list of intermediate f90 files generated by preprocessor
FPPFILES = $(wildcard $(addsuffix /*.fpp, $(SRCDIR) $(PROGDIR) $(SUBDIRS)))
FPPFILES += $(wildcard $(addsuffix /*.F90, $(SRCDIR) $(PROGDIR) $(SUBDIRS)))
F90FROMFPP := $(addprefix $(GENERATED_SRC_DIR)/, $(patsubst %.F90,%.f90,$(patsubst %.fpp,%.f90,$(notdir $(FPPFILES)))))

####################################################################### RULES

.SUFFIXES:

.PHONY: all $(GK_PROJECT)_all submodules
.PHONY: depend tar tar_exec test_make diagnostics simpledataio
.PHONY: help TAGS
.PHONY: clean cleanlib cleanconfig clean_coverage distclean unlink
.PHONY: doc docs/html docs/graphs cleandoc check_ford_install namelist_doc

# comment this out to keep intermediate .f90 files
ifdef KEEP_PROCESSED
.PRECIOUS: $(F90FROMFPP)
endif

####################################################################### Source code targets
$(BINDIR) $(LIBDIR) $(OBJDIR) $(MODDIR) $(GENERATED_SRC_DIR):
	@mkdir -p $@

%.o : %.pf
	$(PFUNIT_PREPROCESSOR) $< $(basename $<).F90
ifdef ONE_STEP_PP
	$(QUIETSYM)$(FC) $(FPPFLAGS) $(F90FLAGS) -c $(basename $<).F90 $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@
else
	$(QUIETSYM)$(CPP) $(FPPFLAGS) $(INCLUDES) $(basename $<).F90 -o $(basename $<).f90
	$(QUIETSYM)$(FC) $(F90FLAGS) -c $(basename $<).f90 $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@
endif

# On many systems the two steps of preprocessing and compiling
# can be done simultaneously. This has many advantages, the
# most important of which is correctly reporting line numbers in 
# the .fpp when debugging. 
#
# On some systems (e.g. Mac OS X) using cpp doesn't work at all
# and you HAVE to preprocess and compile in one step.
ifdef ONE_STEP_PP
$(OBJDIR)/%.o: %.fpp |$(MODDIR) $(OBJDIR)
	$(call messageOp,"Compiling and preprocessing","\\t:",$(call processName,$<))
	$(QUIETSYM)$(FC) $(FPPFLAGS) $(F90FLAGS) $(MODDIR_FLAG)$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o: %.F90 |$(MODDIR) $(OBJDIR)
	$(call messageOp,"Compiling and preprocessing","\\t:",$(call processName,$<))
	$(QUIETSYM)$(FC) $(FPPFLAGS) $(F90FLAGS) $(MODDIR_FLAG)$(MODDIR) -c $< -o $@
endif

$(OBJDIR)/%.o:$(GENERATED_SRC_DIR)/%.f90 |$(MODDIR) $(OBJDIR)
	$(call messageOp,"Compiling","\\t:",$(call processName,$<))
	$(QUIETSYM) $(FC) $(F90FLAGS) $(MODDIR_FLAG)$(MODDIR) -c $< -o $@

$(OBJDIR)/%.o:%.f90 |$(MODDIR) $(OBJDIR)
	$(call messageOp,"Compiling","\\t:",$(call processName,$<))
	$(QUIETSYM) $(FC) $(F90FLAGS) $(MODDIR_FLAG)$(MODDIR) -c $< -o $@

# The depend command needs the preprocessed source files. 
# This is a little clunky, and maybe could be fixed in the future.
$(GENERATED_SRC_DIR)/%.f90: %.fpp |$(GENERATED_SRC_DIR)
	$(call messageOp,"Preprocessing","\\t:",$(call processName,$<))
	$(QUIETSYM)$(CPP) $(FPPFLAGS) $(INCLUDES) $< -o $@

$(GENERATED_SRC_DIR)/%.f90: %.F90 |$(GENERATED_SRC_DIR)
	$(call messageOp,"Preprocessing","\\t:",$(call processName,$<))
	$(QUIETSYM)$(CPP) $(FPPFLAGS) $(INCLUDES) $< -o $@

$(OBJDIR)/%.o: %.c |$(OBJDIR)
	$(call messageOp,"Compiling C","\\t:",$(call processName,$<))
	$(QUIETSYM)$(CC) $(CFLAGS) $(INCLUDES) $(CPPFLAGS) -c $< -o $@

# This prevents error messages like m2c: Command not found
# This is from a implicit makefile rule that we're overriding
%.o : %.mod


DEPEND:=Makefile.depend
DEPEND_CMD:=$(PERL) scripts/fortdep

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

depend: $(F90FROMFPP) $(UTILS_SENTINEL) Makefile.$(GK_SYSTEM) |$(GENERATED_SRC_DIR)
	$(call messageOp,"Generating dependency file","\\t:",$(call processName,$(DEPEND)))
# Note we force UTILS to be included in the paths that we pass to fortdep
# We touch the dependency file first as fortdep will call make in order to try
# and preprocess files, so we need this to avoid recursion
	$(QUIETSYM)touch $(DEPEND) && $(DEPEND_CMD) -m "$(MAKE) --no-print-directory" -1 -o -v=0\
		--keep=$(KEEP_PROCESSED_FILES_FROM_DEPEND) $(VPATH):$(UTILS) \
		--generated_src_dir=$(GENERATED_SRC_DIR) || rm $(DEPEND)
	$(call messageOp,"Dependency file created at","\\t:",$(call processName,$(DEPEND)))

# Dump the compilation flags to a file, so we can check if they change between
# invocations of `make`. The `cmp` bit checks if the file contents
# change. Adding a dependency of a file on `compiler_flags` causes it to be
# rebuilt when the flags change. Taken from
# https://stackoverflow.com/a/3237349/2043465
COMPILER_FLAGS_CONTENTS = "FC = $(FC)\nFPPFLAGS = $(FPPFLAGS)\nF90FLAGS = $(F90FLAGS)"
.PHONY: force
compiler_flags: force
	$(QUIETSYM)echo -e $(COMPILER_FLAGS_CONTENTS) | cmp -s - $@ || echo -e $(COMPILER_FLAGS_CONTENTS) > $@

# Things that should be rebuilt if the compilation flags change (for example, on a new commit)
$(OBJDIR)/git_hash.o: compiler_flags

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

# Here we list all targets for which we don't care if there's a dependency file or not
# this can help us avoid the slow generation of the dependency file when not required
DONT_INCLUDE_DEPEND_TARGETS := clean distclean depend test_make submodules doc namelist_doc
# Note we really need to say do any of our targets not match the above as if we have multiple
# targets and any of them require the dependency file then we have to build the dependency file
# We achieve that using filter-out. Because we also want to be able to build when the user
# doesn't specify the target (i.e. just make -j ) it's not enough to check that MAKECMDGOALS is
# empty after filtering. We need to know if this was the case before hand.
ifneq ($(strip $(MAKECMDGOALS)),)
MUST_INCLUDE_DEPEND := $(filter-out $(DONT_INCLUDE_DEPEND_TARGETS),  $(MAKECMDGOALS))
else
MUST_INCLUDE_DEPEND := "YES"
endif

# If MUST_INCLUDE_DEPEND isn't empty then we better include the dependency file
ifneq ($(strip $(MUST_INCLUDE_DEPEND)),)
include $(DEPEND)
endif

#If we're using fft(w) then we want to link the fft_save_wisdom object file
#this currently has to be handled specially as it comes from a c source file
#rather than fortran (so isn't picked up by fortdep).
ifdef USE_FFT
    utils_mod += $(OBJDIR)/fft_save_wisdom.o
    gs2_mod += $(OBJDIR)/fft_save_wisdom.o
    ingen_mod += $(OBJDIR)/fft_save_wisdom.o
    dump_response_mod += $(OBJDIR)/fft_save_wisdom.o
    generate_fftw_wisdom_mod += $(OBJDIR)/fft_save_wisdom.o
    write_default_input_file_mod += $(OBJDIR)/fft_save_wisdom.o
endif

geometry_c_interface.o: geometry_c_interface.h

gryfx_libs: utils.a geo.a $(OBJDIR)/geometry_c_interface.o

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
clean: cleanobjs cleanpreprocess

cleanobjs:
	$(call messageRm,"Object and module files")
	-$(QUIETSYM)rm -rf $(OBJDIR)/*.o $(MODDIR)/*.mod *.o *.mod

cleandirs:
	$(call messageRm,"Directories")
	-$(QUIETSYM)rm -rf $(MODDIR) $(OBJDIR) $(BINDIR) $(LIBDIR) $(GENERATED_SRC_DIR) $(COVERAGE_DIR) $(TEST_RESULT_DIR)

cleanpreprocess:
	$(call messageRm,"Processed files")
	-$(QUIETSYM)rm -f $(F90FROMFPP)

cleanlib:
	$(call messageRm,"Libraries")
	-$(QUIETSYM)rm -f $(LIBDIR)/*.a

cleanconfig:
	$(call messageRm,"System config")
	-$(QUIETSYM)rm -f system_config .tmp_output

clean_depend:
	$(call messageRm,"Dependency file")	
	-$(QUIETSYM)rm -f $(DEPEND)

distclean: unlink clean cleanlib clean_tests clean_benchmarks cleanconfig cleandoc clean_depend pfunit_clean cleandirs

unlink:
	$(call messageRm,"Preprocessed files")	
	-$(QUIETSYM)rm -f $(F90FROMFPP)

############################################################# Tar targets
tar:
	@[ ! -d $(TARDIR) ] || echo "ERROR: directory $(TARDIR) exists. Stop."
	@[ -d $(TARDIR) ] || $(MAKE) tar_exec

$(TARDIR):
	$(QUIETSYM)$(MKDIR) $@

### setting tar_exec local $(TARLIST*) variables
# expand wildcards listed $(TARLIST_wild) in ( $(TARLIST_dir) + . )
# directories and add them into TARLIST
tar_exec: TARLIST = scripts/fortdep
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
$(TEST_DIR)/Makefile.tests_and_benchmarks: ;
sinclude $(TEST_DIR)/Makefile.tests_and_benchmarks

#This is the location of the pfunit test suites
PFUNIT_TEST_DIR:=$(TEST_DIR)/pfunit_tests
$(PFUNIT_TEST_DIR)/Makefile: ;
sinclude $(PFUNIT_TEST_DIR)/Makefile

############################################################# Introspection targets
test_make:
	$(call messageInfo,"Global settings:")
	$(call messageKeyVal, "GK_HEAD_DIR", $(GK_HEAD_DIR))
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
	$(call messageKeyVal, "GIT_HASH_MAKEFILES", $(GIT_HASH_MAKEFILES))
	$(call messageKeyVal, "GIT_HASH_UTILS", $(GIT_HASH_UTILS))
	$(call messageKeyVal, "GIT_STATE", $(GIT_STATE))
	$(call messageKeyVal, "GIT_BRANCH", $(GIT_BRANCH))
	$(call messageKeyVal, "PATCHFILE", $(PATHFILE))
	$(call messageKeyVal, "IS_RELEASE", $(IS_RELEASE))
	$(call messageKeyVal, "RELEASE", $(RELEASE))
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
	$(call messageKeyVal, "USE_POSIX", $(USE_POSIX))
	$(call messageKeyVal, "USE_LAPACK", $(USE_LAPACK))
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
	$(call messageKeyVal, "FPPFLAGS", $(FPPFLAGS))
	$(call messageKeyVal, "ARCHFLAGS", $(ARCHFLAGS))
	$(call messageKeyVal, "MACRO_DEFS", $(MACRO_DEFS))
	$(call messageKeyVal, "LD_LIBRARY_PATH", $(LD_LIBRARY_PATH))
	$(call messageKeyVal, "LIBS", $(LIBS))
	$(call messageKeyVal, "INCLUDES", $(INCLUDES))
	@echo
	$(call messageInfo,"Build info:")
	$(call messageKeyVal, "SRCDIR ", $(SRCDIR))
	$(call messageKeyVal, "PROGDIR", $(PROGDIR))
	$(call messageKeyVal, "BINDIR ", $(BINDIR))
	$(call messageKeyVal, "LIBDIR ", $(LIBDIR))
	$(call messageKeyVal, "OBJDIR ", $(OBJDIR))
	$(call messageKeyVal, "MODDIR ", $(MODDIR))
	$(call messageKeyVal, "SUBDIRS", $(SUBDIRS))
	$(call messageKeyVal, "GENERATED_SRC_DIR", $(GENERATED_SRC_DIR))
	$(call messageKeyVal, "COVERAGE_DIR", $(COVERAGE_DIR))
	$(call messageKeyVal, "TEST_RESULT_DIR", $(TEST_RESULT_DIR))
	$(call messageKeyVal, "FPPFILES", $(FPPFILES))
	$(call messageKeyVal, "F90FROMFPP", $(F90FROMFPP))
	@echo
	$(call messageInfo,"PYTHON:")
	$(call messageKeyVal, "GS2_HAS_PYTHON", $(GS2_HAS_PYTHON))
	$(call messageKeyVal, "GS2_PYTHON_NAME", $(GS2_PYTHON_NAME))
	$(call messageKeyVal, "PYTHON_EXE", $(PYTHON_EXE))
	$(call messageKeyVal, "PYTHON_VERSION", $(PYTHON_VERSION))
	@echo
	$(call messageInfo,"PFUNIT:")
	$(call messageKeyVal, "PFUNIT_VERSION_MAJOR", $(PFUNIT_VERSION_MAJOR))
	$(call messageKeyVal, "PFUNIT_PREPROCESSOR", $(PFUNIT_PREPROCESSOR))
	$(call messageKeyVal, "PFUNIT_TESTS", $(PFUNIT_TESTS))
	$(call messageKeyVal, "PFUNIT_NPROC", $(PFUNIT_NPROC))
	$(call messageKeyVal, "PFUNIT_MOD", $(PFUNIT_MOD))
	$(call messageKeyVal, "PFUNIT_TEST_DIR", $(PFUNIT_TEST_DIR))
	$(call messageKeyVal, "PFUNIT_TEST_EXE", $(PFUNIT_TEST_EXE))
	$(call messageKeyVal, "PFUNIT_DRIVER", $(PFUNIT_DRIVER))
	$(call messageKeyVal, "PFUNIT_INC", $(PFUNIT_INC))
	$(call messageKeyVal, "PFUNIT_LIB", $(PFUNIT_LIB))
	@echo
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
ifdef USE_LAPACK
	$(call messageInfo,"Lapack:")
	$(call messageKeyVal, "LAPACK_LIB", $(LAPACK_LIB))
	$(call messageKeyVal, "LAPACK_INC", $(LAPACK_INC))
else
	$(call messageInfo,"Lapack: DISABLED")
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

namelist_doc: PATHS:=$(foreach dir, $(SRCDIR) $(SUBDIRS), --path ${dir})
namelist_doc: FILES:=$(shell $(GK_HEAD_DIR)/scripts/get_all_matching_files.py --pattern "*.[fF]??" ${PATHS} --absolute --exclude "$(GK_HEAD_DIR)/src/templates")
# for the `doc` target, `check_ford_install` must be before `namelist_doc` to ensure visibility of the error message
# in the case where Ford is not installed and `make` is run with the `-j` option. The most robust way to do this is
# to have `namelist_doc` to depend explicitly on `check_ford_install` as below
namelist_doc: check_ford_install
	$(call messageInfo,"Generating namelist documentation")
	@$(GK_HEAD_DIR)/scripts/generate_namelist_documentation.py ${FILES}

.PHONY: clean_namelist_doc
clean_namelist_doc:
	$(call messageRm,"Generated namelist documentation")
	-$(QUIETSYM)rm -f docs/page/namelists/auto_generated_* docs/page/namelists/index.md

docs/html:
doc: docs/docs.config.md namelist_doc  check_ford_install
	$(FORD) $(subst -D,--macro ,$(MACRO_DEFS)) $(INCLUDES) -r $(GIT_HASH) $<

cleandoc: clean_namelist_doc
cleandoc: docs/html docs/graphs
	$(call messageRm,"FORD docs")	
	-$(QUIETSYM)rm -rf $^
