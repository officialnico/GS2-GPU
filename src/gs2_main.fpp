!> This module provides the external interface to gs2. It contains functions
!! to initialize gs2, functions to run gs2, functions to finalize gs2 and 
!! functions to override/tweak gs2 parameters. 
!!
!! gs2_main contains both the new and the old interface to gs2. The new interface
!! comprises the gs2_program_state_type object, and functions for manipulating
!! that object, which are documented individually. 
!!
!! The basic flow of a gs2 program can be seen in gs2.fpp, reproduced here:
!!
!!     type(gs2_program_state_type) :: state
!!     call initialize_gs2(state)
!!     call initialize_equations(state)
!!     call initialize_diagnostics(state)
!!     if (state%do_eigsolve) then 
!!       call run_eigensolver(state)
!!     else
!!       call evolve_equations(state, state%nstep)
!!     end if
!!     call finalize_diagnostics(state)
!!     call finalize_equations(state)
!!     call finalize_gs2(state)
!!
!! You can manipulate the gs2_program_state before running gs2. A typical
!! manipulation is to pass in an external mpi communicator (which also 
!! means that MPI_Init is not called within initialize_gs2. 
!!
!!     state%mp_comm_external = .true
!!     state%mp_comm = my_mpi_communicator
!!     call initialize_gs2(state)
!!     ...
!! 
!! You can also override parameters: typically this is done either before
!! calling initialize_equations, or between two calls to evolve equations;
!! here we manipulate the temperature gradient of the first species:
!!
!!     call prepare_profiles_overrides(state)
!!     state%init%prof_ov%override_tprim(1) = .true.
!!     state%init%prof_ov%tprim(1) = 5.5
!!     call initialize_equations(state)
!!     call initialize_diagnostics(state)
!!     call evolve_equations(state, state%nstep/2)
!!     call prepare_profiles_overrides(state)
!!     state%init%prof_ov%tprim(1) = 6.0
!!     call initialize_equations(state)
!!     call evolve_equations(state, state%nstep/2)
!!     ...
!!   
!! It is very important to note that just because this interface is presented
!! in an object-oriented way, it does not, unfortunately, mean that the entire
!! program state, i.e. data arrays etc, is encapsulated in the gs2_program_state 
!! object. In fact most data is stored globally as module public variables. 
!! The gs2_program_state object is provided for convenience, as a way of
!! monitoring the execution state of gs2, and because it is hoped that in the 
!! future gs2 will, in fact, be thread safe and have proper data encapsulation.
!! A particular consequence of this is that only one instance of the
!! gs2_program_state object should exist on each process, i.e. in a given 
!! memory space.
module gs2_main
  use gs2_init, only: init_type, init_level_list
  use optimisation_config, only: optimisation_type
  use constants, only: run_name_size
  use exit_codes, only: exit_code, EXIT_NOT_REQUESTED
  implicit none
  public :: parse_command_line

  public :: gs2_program_state_type

  public :: initialize_gs2
  public :: initialize_equations
  public :: initialize_diagnostics, evolve_equations, run_eigensolver
  public :: finalize_diagnostics, finalize_equations, finalize_gs2
  public :: calculate_outputs
  public :: reset_equations
  public :: write_used_inputs_file

  public :: prepare_miller_geometry_overrides
  public :: prepare_profiles_overrides
  public :: prepare_kt_grids_overrides
  public :: prepare_optimisations_overrides
  public :: prepare_initial_values_overrides
  public :: set_initval_overrides_to_current_vals
  public :: finalize_overrides

  public :: initialize_wall_clock_timer


  !> FIXME : Add documentation
  type gs2_timers_type
    real :: init(2) 
    real :: advance(2) = 0.
    real :: timestep(2) = 0.
    real :: finish(2) = 0.
    real :: total(2) = 0. 
    real :: diagnostics(2)=0.
    real :: eigval(2)=0.
    !real :: interval
    real :: main_loop(2)
  end type gs2_timers_type

  !> A type for storing outputs of gs2
  !! for access externally
  type gs2_outputs_type
    !> The gradient of the flux tube volume wrt
    !! the flux label: related to the surface
    !! area of the flux tube
    real :: dvdrho
    !> The average gradient of the flux tube
    !! label <|grad rho|>
    real :: grho
    !> Particle flux by species
    real, dimension (:), allocatable :: pflux 
    !> Heat flux by species
    real, dimension (:), allocatable :: qflux 
    !> Turbulent heating by species
    real, dimension (:), allocatable :: heat 
    !> Momentum flux
    real :: vflux
  end type gs2_outputs_type

  !> The object which specifies and records the gs2 program state.
  !! Some attributes of the object, like mp_comm_external, are
  !! designed to directly manipulated by the user, and some are
  !! designed to store program state information and be 
  !! manipulated internally, like init\%level.
  type gs2_program_state_type

    ! Flags indicating the current state of the 
    ! program (used for error checking)
    logical :: gs2_initialized = .false.
    logical :: equations_initialized = .false.
    logical :: diagnostics_initialized = .false.

    !> A type for keeping track of the current
    !! initialization level of gs2, as well
   !! as storing all the overrides. See 
    !! gs2_init::init_type for more information.
    type(init_type) :: init


    !> Do not set manually. If fewer than the 
    !! available number of processors are being
    !! used, this is true for the processors
    !! that are active and false for those that
    !! lie idle.
    logical :: included = .true.

    ! Timers
    type(gs2_timers_type) :: timers

    !> The exit flag is set to true by any 
    !! part of the main timestep loop that 
    !! wants to cause the loop to exit
    logical :: exit = .false.

    !> Whether the run has converged to a
    !! stationary state
    logical :: converged = .false.

    !> Whether to run the eigenvalue solver
    !! or not. Set equal to the input value
    !! in intialize_equations. Typically
    !! only important for the standalone
    !! gs2 program
    logical :: do_eigsolve = .false.

    !> This parameter is set equal to run_parameters::nstep in 
    !! initialize_equations and is the maximum
    !! number of timesteps that can be run 
    !! without reinitalising the diagnostics. 
    !! Do not modify!
    integer :: nstep

    !> Gets set to the final value of istep
    !! in evolve equations. Any future calls
    !! to evolve_equations will increment this
    !! further. A call to finalize_diagnostics
    !! will set it back to 1. Note that setting
    !! this manually is NOT advised. 
    integer :: istep_end = 0

    !> Whether to print out debug messages
    !logical :: debug = .false.
    integer :: verb = 3

    !> Parameters pertaining to cases when gs2
    !! is being used as library.
    !! external_job_id is not to confused with the parameter
    !! job in mp, which identifies the subjob if
    !! running in list mode or with nensembles > 1
    !integer :: trin_job = -1
    integer :: external_job_id = -1

    !> is_external_job should be set to true when GS2
    !! is being used as a library 
    logical :: is_external_job = .false.

    !> Set true if using trinity. This does several things: 
    !!  * it forces the calculation of the fluxes
    !!  * it causes the species and theta_grid name lists to use parameters
    !!     provided by trinity
    !! Setting this flag true automatically sets is_external_job 
    !! to true
    logical :: is_trinity_job = .false.

    integer :: trinity_timestep = -1
    integer :: trinity_iteration = -1


    !> If true, don't initialize or evolve the equations, 
    !! just read their values from an existing output file.
    !! This allows additional calculations to be done, and 
    !! also allows GS2 to cheaply reproduce its functionality,
    !! e.g. within Trinity.
    logical :: replay = .false.

    !> If true, the equations are fully initialized during
    !! replays. This is typically needed for linear diffusive
    !! flux estimates, and is set true automatically in that 
    !! case.
    logical :: replay_full_initialize = .false.

    !> If true, 
    !! print full timing breakdown. 
    logical :: print_full_timers = .false.

    !> If true, print run time or
    !! full timing breakdown, depending
    !! on the value of print_full_timers
    logical :: print_times = .true.


    !> Parameters to be used when passing in an external communicator
    logical :: mp_comm_external = .false.
    integer :: mp_comm

    !> This is set in initialize_gs2 to the number
    !! of procs actually used
    integer :: nproc_actual

    logical :: run_name_external = .false.
    character(run_name_size) :: run_name

    logical :: skip_diagnostics = .false.
    logical :: dont_change_timestep = .false.

    !> Whether this is a list mode run
    logical :: list
    !> The number of identical runs happening
    !! simultaneously (used for ensemble averaging).
    !! Cannot be used in conjunction with list mode
    integer :: nensembles = 1

    ! Outputs (e.g. for Trinity)
    type(gs2_outputs_type) :: outputs

    ! Optimisation configuration
    type(optimisation_type) :: optim

    !> Reason for end of simulation
    type(exit_code) :: exit_reason = EXIT_NOT_REQUESTED
   
  end type gs2_program_state_type

  private
 
contains

  !> Parse some basic command line arguments. Currently just 'version' and 'help'.
  !>
  !> This should be called before anything else, but especially before initialising MPI.
  subroutine parse_command_line()
    use git_version_mod, only: get_git_version
    use build_config, only : formatted_build_config
    integer :: arg_count, arg_n
    integer :: arg_length
    character(len=:), allocatable :: argument
    character(len=*), parameter :: nl = new_line('a')
    character(len=*), parameter :: usage = &
         "gs2 [--version|-v] [--help|-h] [--build-config] [input file]" // nl // nl // &
         "GS2: A fast, flexible physicist's toolkit for gyrokinetics" // nl // &
         "For more help, see the documentation at https://gyrokinetics.gitlab.io/gs2/" // nl // &
         "or create an issue https://bitbucket.org/gyrokinetics/gs2/issues?status=open" // nl // &
         nl // &
         "  -h, --help     Print this message" // nl // &
         "  -v, --version  Print the GS2 version" // nl // &
         "  --build-config Print the current build configuration"

    arg_count = command_argument_count()

    do arg_n = 0, arg_count
      call get_command_argument(arg_n, length=arg_length)
      if (allocated(argument)) deallocate(argument)
      allocate(character(len=arg_length)::argument)
      call get_command_argument(arg_n, argument)

      if ((argument == "--help") .or. (argument == "-h")) then
        write(*, '(a)') usage
        stop
      else if ((argument == "--version") .or. (argument == "-v")) then
        write(*, '("GS2 version ", a)') get_git_version()
        stop
      else if (argument == "--build-config") then
        write(*, '(a)') formatted_build_config()
        stop
     else if (argument == "--") then
        exit
      else if (argument(1:1) == "-") then
        write(*, '(a)') "Error: Unrecognised argument '" // argument // "'. Usage is:" // nl // nl // usage
        stop 2
      end if
    end do
  end subroutine parse_command_line

  !> Starts the global wall clock timer
  !! used by check time. This is useful
  !! if you have multiple copies of gs2 
  !! running but you don't want to start 
  !! them all at the same time, but you
  !! want them all to respect avail_cpu_time
  subroutine initialize_wall_clock_timer
    use job_manage, only: init_checktime
    call init_checktime
  end subroutine initialize_wall_clock_timer

  !> Initialise message passing, open the input file, split the 
  !! communicator if running in list mode or if nensembles > 1,
  !! initialize the timers. After calling this function, gs2
  !! reaches init\%level = basic. If it is desired to provide
  !! an external commuicator or set the input file name externally,
  !! these should be set before calling this function. 
  subroutine initialize_gs2(state, header, quiet)
    use file_utils, only: init_file_utils
    use file_utils, only: run_name, run_name_target
    use gs2_init, only: init_gs2_init
    use job_manage, only: checktime, time_message
    use job_manage, only: init_checktime, checktime_initialized
    use job_manage, only: job_fork
    use mp, only: init_mp, broadcast
    use mp, only: iproc, nproc, proc0
    use mp, only: trin_flag, job
    use mp, only: use_nproc, included, mp_comm
    use redistribute, only: using_measure_scatter
    use unit_tests, only: debug_message, set_job_id
    use standard_header, only: standard_header_type
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    !> Header for files with build and run information
    type(standard_header_type), intent(in), optional :: header
    !> If true, don't print start-up messages/header
    logical, intent(in), optional :: quiet

    integer :: oldrank
    ! Actual value for optional `header` input
    type(standard_header_type) :: local_header
    ! Internal value for optional `quiet` input
    logical :: quiet_value

    if (state%init%level .ge. init_level_list%basic) then
      write  (*,*) "ERROR: Called initialize_gs2 twice &
        & without calling finalize_gs2"
      stop 1
    end if

    call debug_message(4, 'gs2_main::initialize_gs2 starting initialization')

    if (state%mp_comm_external) then
       call init_mp (state%mp_comm)
    else
       call init_mp
    end if

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    if (present(quiet)) then
      quiet_value = quiet
    else
      quiet_value = .false.
    end if

    oldrank = iproc
    if (state%init%opt_ov%is_initialised() .and. state%init%opt_ov%override_nproc) then
      state%init%opt_ov%old_comm = mp_comm
      call use_nproc(state%init%opt_ov%nproc)
    end if
    state%included = included
    state%nproc_actual = nproc
    if (.not. state%included) return

    if (state%is_trinity_job) state%is_external_job = .true.
    if (state%is_external_job) then 
      call broadcast(state%external_job_id)
      call set_job_id(state%external_job_id)
    end if

    trin_flag = .false. 

    call debug_message(state%verb, 'gs2_main::initialize_gs2 initialized mp')

    call reset_timers(state%timers)
    if (.not. checktime_initialized) call init_checktime ! <doc> Initialize timer </doc>
    call debug_message(state%verb, &
      'gs2_main::initialize_gs2 called init_checktime')
  
    if (proc0) then
      ! Report number of processors being used
      if (.not. quiet_value) then
        write(*, '(a)') local_header%to_string(comment_character=" ")

        if (state%external_job_id >= 0) then
          write(*, '(a, i0, a)', advance="no") "External job ", state%external_job_id, " "
        end if
        write(*, '(a, i0, a)', advance="no") "Running on ", nproc, " processor"
        if (nproc > 1) then
          ! Pluralise processor
          write(*, '(a)') "s" // new_line('a')
        else
          ! Make sure we still get a blank line
          write(*,*) new_line('a')
        end if
      end if

      ! Call init_file_utils, ie. initialize the run_name, checking
      !  whether we are doing a Trinity run or a list of runs.
      ! Until the logic of init_file_utils is fixed we set trin_run
      ! when an external file name has been provided... this prevents
      ! it overriding the name from the command line.
      if (state%run_name_external) then
        call init_file_utils (state%list, trin_run=state%run_name_external, &
             name=state%run_name, n_ensembles=state%nensembles)
      else
        call init_file_utils (state%list, name="gs")
      end if
      call debug_message(state%verb, &
           'gs2_main::initialize_gs2 initialized file_utils')
    end if
    
    call broadcast (state%list)
    call debug_message(state%verb, 'gs2_main::initialize_gs2 broadcasted list')

    if (state%list) then
       call job_fork
       call set_job_id(job)
       call debug_message(state%verb, 'gs2_main::initialize_gs2 called job fork')
    else if (state%nensembles > 1) then 
       call job_fork (n_ensembles=state%nensembles)
       call debug_message(state%verb, &
         'gs2_main::initialize_gs2 called job fork (nensembles>1)')
    end if

    call time_message(.false.,state%timers%total,' Total')


    ! Pass run name to other procs
    if (proc0) then
       run_name_target = trim(run_name)
    end if
    call broadcast (run_name_target)
    if (.not. proc0) run_name => run_name_target
    call debug_message(state%verb, 'gs2_main::initialize_gs2 run_name = '//run_name)

    !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
    using_measure_scatter=.false.

    ! Initialize the gs2 initialization system
    call init_gs2_init(state%init)

    !state%gs2_initialized = .true.
    state%init%level = init_level_list%basic

    call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, 'gs2_main::initialize_gs2 finished')

  end subroutine initialize_gs2

  !> Initialize all the modules which are used to evolve the 
  !! equations. After calling this function, gs2 reaches 
  !! init\%level = full. 
  subroutine initialize_equations(state)
    use fields, only: init_fields
    use theta_grid, only: surfarea, dvdrhon
    use gs2_time, only: init_tstart
    use gs2_init, only: init
    use init_g, only: tstart
    use job_manage, only: time_message
    use mp, only: proc0, broadcast
    use parameter_scan, only: init_parameter_scan
    use run_parameters, only: nstep, do_eigsolve
    use unit_tests, only: debug_message
    use gs2_reinit, only: init_gs2_reinit
    use fields_implicit, only: skip_initialisation
    use fields, only: fields_allocate_arrays=>allocate_arrays
    use fields_local, only: fieldmat
    use run_parameters, only: trinity_linear_fluxes
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    call debug_message(state%verb, 'gs2_main::initialize_equations starting')
    call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, &
      'gs2_main::initialize_equations calling init_parameter_scan')
    call init_parameter_scan

    call time_message(.false., state%timers%init,' Initialization')
    !Set using_measure_scatter to indicate we want to use in "gather/scatter" timings
    call debug_message(state%verb, 'gs2_main::initialize_equations calling init_fields')

    !> If we are doing a replay then we are not evolving the equations
    !! so we disable calculating the implicit response
    !! matrix
    if (state%replay) then
      skip_initialisation = .true.
      fieldmat%no_prepare = .true.
      fieldmat%no_populate = .true.
      call init(state%init, init_level_list%le_grids)

      if (trinity_linear_fluxes) then
        state%replay_full_initialize = .true.
        call init(state%init, init_level_list%full)
      end if
        
    else
      ! This triggers initializing of all the grids, all the physics parameters
      ! and all the modules which solve the equations
      call init(state%init, init_level_list%full)
    end if

    ! Set the initial simulation time (must be after init_fields
    ! because initial time may be read from a restart file)
    call init_tstart(tstart)


    ! Here we copy some geometric information required by 
    ! Trinity to state%outputs
    if (proc0) then
       state%outputs%dvdrho = dvdrhon
       ! Note surfarea isn't defined/set for s-alpha geometry so
       ! the following is not well defined.
       ! Note this isn't actually "grho" as defined by the rest of the
       ! code but rather "grhoavg", i.e. an average value. It may be
       ! appropriate to replace this with something like `mean(grho)`
       state%outputs%grho = surfarea / dvdrhon
    end if
    call broadcast (state%outputs%dvdrho)
    call broadcast (state%outputs%grho)
    
    call time_message(.false.,state%timers%init,' Initialization')


    ! Set defaults. These are typically only important
    ! for the standalone gs2 program
    state%nstep = nstep
    state%do_eigsolve = do_eigsolve

    call debug_message(state%verb, 'gs2_main::initialize_equations finished')

    ! Since the number of species may have changed
    ! since the last call to initialize_equations 
    ! we reallocate outputs every time.
    if (allocated(state%outputs%pflux)) then
      call deallocate_outputs(state)
    end if
    call allocate_outputs(state)

    call init_gs2_reinit

    call time_message(.false.,state%timers%total,' Total')

  end subroutine initialize_equations

  !> Initialize the diagostics modules. This function can only be
  !! called when the init level is 'full', i.e., after calling 
  !! initialize_equations.  However, gs2 can be partially
  !! uninitialized without finalizing the diagnostics, for example to
  !! change parameters such as the timestep.
  subroutine initialize_diagnostics(state, header)
    use gs2_diagnostics, only: init_gs2_diagnostics
    use gs2_diagnostics, only: nwrite, write_nl_flux
    use gs2_diagnostics, only: loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: diagnostics_init_options_type
    use gs2_diagnostics_new, only: init_gs2_diagnostics_new
    use gs2_diagnostics_new, only: run_diagnostics
#endif
    use job_manage, only: time_message
    use mp, only: proc0, mp_abort, broadcast
    use parameter_scan, only: allocate_target_arrays
    use run_parameters, only: nstep, use_old_diagnostics
    use unit_tests, only: debug_message
    use standard_header, only: standard_header_type
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    !> Header for files with build and run information
    type(standard_header_type), intent(in), optional :: header
    ! Actual value for optional `header` input
    type(standard_header_type) :: local_header

#ifdef NEW_DIAG
    ! Configuration for the new diagnostics module
    type(diagnostics_init_options_type) :: diagnostics_init_options
    real :: precision_test
#endif
    if (.not. state%included) return

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    state%exit = .false. 

    call time_message(.false.,state%timers%total,' Total')

    if (state%init%diagnostics_initialized) &
      call mp_abort('Calling initialize_diagnostics twice', .true.)

#ifndef NEW_DIAG
    if (.not. use_old_diagnostics)then
       if(proc0)then
          write(*,'("  WARNING: Compiled without new diagnostics and runnning with use_old_diagnostics=.false. will lead to no output --> Forcing use_old_diagnostics=.true.")')
          write(*,'("           Please consider recompiling with the new diagnostics enabled.")')
          write(*,*)
       endif
       use_old_diagnostics=.true.
    endif
#endif

    if (.not. use_old_diagnostics) then
#ifdef NEW_DIAG
      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling init_gs2_diagnostics_new')

#ifdef NETCDF_PARALLEL
      diagnostics_init_options%parallel_io_capable = .true.
      diagnostics_init_options%parallel_io_capable = .false.
#else
      diagnostics_init_options%parallel_io_capable = .false.
#endif

      ! Here we check if reals have been promoted to doubles
      diagnostics_init_options%default_double =  (precision(precision_test).gt.10)
      ! Check whether this is a Trinity run... enforces calculation of the
      ! fluxes
      diagnostics_init_options%is_trinity_run = state%is_trinity_job
      ! Specifiy whether we are doing a replay.
      diagnostics_init_options%replay = state%replay

      call init_gs2_diagnostics_new(diagnostics_init_options, local_header)

      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling run_diagnostics')
      ! Create variables and write constants
      call run_diagnostics(-1,state%exit, state%replay)
      !call broadcast(state%exit)
      ! Write initial values
      call run_diagnostics(0,state%exit, state%replay)
#else
      call mp_abort("use_old_diagnostics is .false. but you have &
        & not built gs2 with new diagnostics enabled", .true.)
#endif
    else ! if use_old_diagnostics

      call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics calling init_gs2_diagnostics')
      call init_gs2_diagnostics (state%list, nstep, header=local_header)
      call allocate_target_arrays(nwrite,write_nl_flux) ! must be after init_gs2_diagnostics
      call loop_diagnostics(0,state%exit)
    end if

    state%init%diagnostics_initialized = .true.

    call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, &
        'gs2_main::initialize_diagnostics finished')

  end subroutine initialize_diagnostics



  !> Run the initial value solver. nstep_run must
  !! be less than or equal to state\%nstep, which is 
  !! set from the input file. The cumulative number of
  !! steps that can be run cannot exceed state\%nstep, 
  !! without a call to reset_equations.
  !! Examples:
  !!
  !!      call evolve_equations(state, state%nstep) ! Basic run
  !!
  !! Note that these two calls: 
  !!  
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!
  !! will in general have the identical effect to the single call above.
  !! There are circumstances when this is not true,
  !! for example, if the first of the two calls to evolve_equations  
  !! exits without running nstep/2 steps (perhaps because the
  !! growth rate has converged). 
  !!
  !! This example will cause an error because the total number
  !! of steps exceeds state\%nstep:
  !!
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!      call evolve_equations(state, state%nstep/2)
  !!
  !! This is OK:
  !!
  !!      call evolve_equations(state, state%nstep/2) 
  !!      call evolve_equations(state, state%nstep/2)
  !!      call reset_equations(state)
  !!      call evolve_equations(state, state%nstep/2)
  subroutine evolve_equations(state, nstep_run, header)
    use collisions, only: vnmult
    use dist_fn_arrays, only: gnew
    use fields, only: advance
    use gs2_diagnostics, only: nsave, loop_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: run_diagnostics
    use gs2_diagnostics_new, only: gnostics
#endif
    use gs2_reinit, only: reset_time_step
    use gs2_reinit, only: check_time_step
    use gs2_save, only: gs2_save_for_restart
    use gs2_time, only: user_time, user_dt, update_time, write_dt
    use job_manage, only: time_message, checkstop, checktime
    use mp, only: proc0
    use mp, only: scope, subprocs, broadcast
    use parameter_scan, only: update_scan_parameter_value
    use run_parameters, only: reset, fphi, fapar, fbpar, nstep, max_sim_time
    use run_parameters, only: avail_cpu_time, margin_cpu_time
    use run_parameters, only: use_old_diagnostics
    use unit_tests, only: ilast_step, debug_message
    use gs2_init, only: init
    use nonlinear_terms, only: gryfx_zonal
    use standard_header, only: standard_header_type
    use exit_codes, only: EXIT_MAX_SIM_TIME
    type(gs2_program_state_type), intent(inout) :: state
    type(standard_header_type), optional, intent(in) :: header
    integer :: istep, istatus
    integer, intent(in) :: nstep_run
    logical :: temp_initval_override_store
    integer :: istep_loop_max, istep_offset
    ! Actual value for optional `header` input
    type(standard_header_type) :: local_header

    if (.not. state%included) return

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    call debug_message(state%verb, &
        'gs2_main::evolve_equations starting')
    
    temp_initval_override_store = state%init%initval_ov%override
    
    call time_message(.false.,state%timers%total,' Total')
    
    if (state%nensembles > 1) &
          call scope (subprocs)

    if (state%is_trinity_job) call write_trinity_parameters(state)

    call time_message(.false.,state%timers%main_loop,' Main Loop')

    ! Make sure exit is false before entering
    ! timestep loop
    state%exit = .false.

    call debug_message(state%verb, &
        'gs2_main::evolve_equations starting loop')
    ! We run for nstep_run iterations, starting from whatever istep we got
    ! to in previous calls to this function. Note that calling
    ! finalize_diagnostics resets state%istep_end
    istep_loop_max = state%istep_end + nstep_run
    do istep = state%istep_end+1, istep_loop_max

       if (istep .gt. nstep) then
         if (proc0) write (*,*) 'Reached maximum number of steps allowed &
           & (set by nstep) without restarting diagnostics.'
         exit
       end if

       call time_message(.false.,state%timers%advance,' Advance time step')

       ! If we are doing a replay we don't actually evolve 
       ! the equations
       if (.not. state%replay) then 

         call debug_message(state%verb+1, &
            'gs2_main::evolve_equations calling advance')

         call time_message(.false.,state%timers%timestep,' Timestep')
         ! Try to take a time step
         call advance (istep)
         call time_message(.false.,state%timers%timestep,' Timestep')
         call debug_message(state%verb+1, &
              'gs2_main::evolve_equations called advance')

         if(state%dont_change_timestep) reset = .false.

         !If the advance call has triggered a reset then change the timestep.
         !Currently this is triggered if immediate_reset is true and the nonlinear
         !cfl limit is exceeded.
         if (reset) then
            call prepare_initial_values_overrides(state)
            call debug_message(state%verb+1, &
                 'gs2_main::evolve_equations resetting timestep')
            call set_initval_overrides_to_current_vals(state%init%initval_ov)
            state%init%initval_ov%override = .true.
            if (state%is_external_job) then
               call reset_time_step (state%init, istep, state%exit, state%external_job_id)
            else
               call reset_time_step (state%init, istep, state%exit)
            end if

            ! Now we've reset set the flag back to .false.
            reset = .false.

            ! If we've changed the timestep and have not decided to exit then
            ! we need to call advance to actually complete this timestep.
            if (.not. state%exit) then
               call advance(istep)
            end if
         end if
         if(state%exit) exit

         call debug_message(state%verb+1, &
            'gs2_main::evolve_equations calling gs2_save_for_restart')
         
         if (use_old_diagnostics) then
           if (nsave > 0 .and. mod(istep, nsave) == 0) &
                call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar, header=local_header)
         else
#ifdef NEW_DIAG
           if (gnostics%nsave > 0 .and. mod(istep, gnostics%nsave) == 0) &
                call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar, header=local_header)
#endif
         end if
         call update_time
         call time_message(.false.,state%timers%diagnostics,' Diagnostics')

       end if ! if (.not. state%replay)
       call debug_message(state%verb+1, &
          'gs2_main::evolve_equations calling diagnostics')
       if (.not. state%skip_diagnostics) then 
         if (use_old_diagnostics) then 
           call loop_diagnostics (istep, state%exit)
#ifdef NEW_DIAG
         else 
           call run_diagnostics (istep, state%exit, state%replay)
#endif
         end if
       end if

       if(state%exit) state%converged = .true.
       if (state%exit) call debug_message(state%verb+1, &
         'gs2_main::evolve_equations exit true after diagnostics')

       call time_message(.false.,state%timers%diagnostics,' Diagnostics')
       call time_message(.false.,state%timers%advance,' Advance time step')

       if(.not.state%exit)then
          if (.not.state%replay) then
            call debug_message(state%verb+1, &
              'gs2_main::evolve_equations checking time step')

            !Note this should only trigger a reset for timesteps too small
            !as timesteps too large are already handled when immediate_reset is
            !.true. If immediate_reset is .false. then this is where we handle
            !all checking of the time step.
            call check_time_step(reset,state%exit)
            !call update_scan_parameter_value(istep, reset, state%exit)
            call debug_message(state%verb+1, &
              'gs2_main::evolve_equations checked time step')

            !If something has triggered a reset then reset here
            if(state%dont_change_timestep) reset = .false.
            if (reset) then
               call prepare_initial_values_overrides(state)
               call set_initval_overrides_to_current_vals(state%init%initval_ov)
               state%init%initval_ov%override = .true.
               ! if called within trinity, do not dump info to screen
               if (state%is_external_job) then
                  call reset_time_step (state%init, istep, state%exit, state%external_job_id)
               else       
                  call reset_time_step (state%init, istep, state%exit)
               end if
            end if
          end if ! if (.not. state%replay)

          if ((mod(istep,5) == 0).and.(.not.state%exit)) call checkstop(state%exit,state%exit_reason)
          if (state%exit) call debug_message(state%verb+1, &
            'gs2_main::evolve_equations exit true after checkstop')
          if (.not.state%exit) call checktime(avail_cpu_time,state%exit,margin_cpu_time)
          if (state%exit) call debug_message(state%verb+1, &
            'gs2_main::evolve_equations exit true after checktime')
       endif
       call debug_message(state%verb+1, &
        'gs2_main::evolve_equations after reset and checks')

       state%istep_end = istep

       ! Check if the next step would exceed the maximum requested simulation time,
       ! and if so flag that we wish to exit. Only check this if we're not already exiting.
       if ((.not. state%exit) .and. user_time + user_dt > max_sim_time) then
          state%exit_reason = EXIT_MAX_SIM_TIME
          state%exit = .true.
       end if

       if (state%exit) then
           call debug_message(state%verb, &
                'gs2_main::evolve_equations exiting loop')
          exit
       end if

    end do

    ! Try to ensure that we write the time dependent diagnostics at the end of the simulation
    ! even in cases where nstep does not have nwrite as a factor.
    if (.not. state%skip_diagnostics) then
       ! Note when we have completed all our iterations we pass istep - 1 as we assume istep has been 
       ! incremented once from the last step actually taken. In cases where we explicitly exit the loop
       ! we need to pass istep instead.
       istep_offset = 1
       if(state%exit) istep_offset = 0
       if (use_old_diagnostics) then
          call loop_diagnostics (istep - istep_offset, state%exit, force = .true.)
#ifdef NEW_DIAG
       else
          call run_diagnostics (istep - istep_offset, state%exit, state%replay, force = .true.)
#endif
       end if
    end if

  ! NRM: When running inside GryfX, evolve_equations is called every timestep to
  ! advance a single timestep.
  ! So anything that is supposed to happen after GS2 completes the main nstep
  ! loop should not be done while running inside GryfX.
  if(.not. gryfx_zonal%on) then
    call time_message(.false.,state%timers%main_loop,' Main Loop')

    if (proc0 .and. .not. state%is_external_job) call write_dt
    
    call time_message(.false.,state%timers%total,' Total')

    if (state%print_times) call print_times(state, state%timers)

    ilast_step = state%istep_end

    ! We also restore the value of the override
    ! switch to what it was when this function was
    ! called.
    state%init%initval_ov%override = temp_initval_override_store

    call debug_message(state%verb,'gs2_main::evolve_equations finished')
  end if


  end subroutine evolve_equations

  !> FIXME : Add documentation
  subroutine run_eigensolver(state)
#ifdef WITH_EIG
    use eigval, only: init_eigval, finish_eigval
    use eigval, only: BasicSolve
#endif 
    use job_manage, only: time_message
    use mp, only: mp_abort, proc0
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
#ifdef WITH_EIG
   call time_message(.false.,state%timers%total,' Total')

   !Start timer
   call time_message(.false.,state%timers%eigval,' Eigensolver')

   !Initialise slepc and the default/input eigensolver parameters
   call init_eigval

   !Create a solver based on input paramters, use it to solve and
   !then destroy it.
   call BasicSolve

   !Tidy up
   call finish_eigval

   !Stop timer
   call time_message(.false.,state%timers%eigval,' Eigensolver')
#else
   call mp_abort("Require slepc/petsc")
#endif
   call time_message(.false.,state%timers%total,' Total')

  end subroutine run_eigensolver

  !> FIXME : Add documentation
  subroutine reset_equations (state)
    use gs2_diagnostics, only: gd_reset => reset_init
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: reset_averages_and_counters
#endif
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes, use_old_diagnostics
    use gs2_time, only: code_dt, save_dt
    use mp, only: scope, subprocs, allprocs
    use mp, only: mp_abort
    use unit_tests, only: debug_message

    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return


    if (state%nensembles > 1) call scope (subprocs)

    call debug_message(state%verb, 'gs2_main::reset_equations starting')

    call save_dt (code_dt)
    ! This call to gs2_diagnostics::reset_init sets some time averages
    ! and counters to zero... used mainly for trinity convergence checks.
    if (use_old_diagnostics) then
      call gd_reset
#ifdef NEW_DIAG
    else 
      call reset_averages_and_counters
#endif
    end if

    if (trinity_linear_fluxes .and. &
       nonlinear_mode_switch .eq. nonlinear_mode_none) &
        call reset_linear_magnitude

    state%istep_end = 0

    call debug_message(state%verb, 'gs2_main::reset_equations finished')

    if (state%nensembles > 1) call scope (allprocs)

  end subroutine reset_equations

  !> FIXME : Add documentation  
  subroutine finalize_diagnostics(state)
    use gs2_diagnostics, only: finish_gs2_diagnostics
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: finish_gs2_diagnostics_new 
#endif
    use job_manage, only: time_message
    use mp, only: proc0, mp_abort
    use parameter_scan, only: deallocate_target_arrays
    use run_parameters, only: use_old_diagnostics
    use unit_tests, only: debug_message
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return

    call time_message(.false.,state%timers%total,' Total')
    call time_message(.false.,state%timers%finish,' Finished run')

    call debug_message(state%verb, 'gs2_main::finalize_diagnostics starting')

    if (.not. state%init%diagnostics_initialized) &
      call mp_abort('Calling finalize_diagnostics when &
      & diagnostics are not initialized', .true.)

    if (.not. use_old_diagnostics) then
#ifdef NEW_DIAG
      call finish_gs2_diagnostics_new 
#endif
    else
      call finish_gs2_diagnostics (state%istep_end)
    end if

    call deallocate_target_arrays

    state%istep_end = 0
    state%init%diagnostics_initialized = .false.

    call time_message(.false.,state%timers%finish,' Finished run')
    call time_message(.false.,state%timers%total,' Total')
  end subroutine finalize_diagnostics

  !> FIXME : Add documentation  
  subroutine finalize_equations(state)
    use gs2_init, only: init
    use job_manage, only: time_message
    use mp, only: proc0
    use parameter_scan, only: finish_parameter_scan
    use unit_tests, only: debug_message
    use gs2_reinit, only: finish_gs2_reinit
    !use fields_arrays, only: phi, apar, bpar, phinew, aparnew, bparnew
    implicit none
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return

    call time_message(.false.,state%timers%finish,' Finished run')
    call time_message(.false.,state%timers%total,' Total')

    call debug_message(state%verb, 'gs2_main::finalize_equations starting')

    call finish_gs2_reinit

    call deallocate_outputs(state)
    call finish_parameter_scan
    !if (state%replay .and. allocated(phi)) deallocate (phi, apar, bpar, phinew, aparnew, bparnew)

    call init(state%init, init_level_list%basic)
    call time_message(.false.,state%timers%finish,' Finished run')
    call time_message(.false.,state%timers%total,' Total')
  end subroutine finalize_equations

  !> FIXME : Add documentation  
  subroutine finalize_gs2(state, quiet)
    use file_utils, only: finish_file_utils
    use gs2_init, only: finish_gs2_init
    use job_manage, only: time_message
    use mp, only: finish_mp, proc0, barrier, unsplit_all, included
    use unit_tests, only: debug_message
    use standard_header, only: date_iso8601
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    !> If true, don't print start-up messages/header
    logical, intent(in), optional :: quiet

    ! Internal value for optional `quiet` input
    logical :: quiet_value

    if (present(quiet)) then
      quiet_value = quiet
    else
      quiet_value = .false.
    end if

    if (state%included) then

      if (state%init%level .ne. init_level_list%basic) then
        write  (*,*) "ERROR: Called finalize_gs2 at the &
          & wrong init_level (perhaps you have called finalize_gs2 &
          & without calling initialize_gs2, or without calling &
          & finalize_equations"
        stop 1
      end if

      !if ((.not. state%gs2_initialized) .or. &
          !state%equations_initialized .or. &
          !state%diagnostics_initialized) then
          !write (*,*) 'ERROR: initialize_gs2 can only be called when &
          !& gs2_initialized is true, and equations_initialized &
          !& and diagnostics_initialized, &
          !& are all false. '
         !stop 1
       !end if

      call time_message(.false.,state%timers%total,' Total')

      call finish_gs2_init(state%init)

      if (proc0) call finish_file_utils

      call time_message(.false.,state%timers%finish,' Finished run')
      call time_message(.false.,state%timers%total,' Total')
      
      call debug_message(state%verb, 'gs2_main::finalize_gs2 calling print_times')

      if (state%print_times .and. .not. quiet_value) call print_times(state, state%timers)

    end if

    if (state%init%opt_ov%is_initialised() .and. state%init%opt_ov%override_nproc) &
      call unsplit_all(state%init%opt_ov%old_comm)
    call barrier
    state%included = included

    call debug_message(state%verb, 'gs2_main::finalize_gs2 calling finish_mp')

    if (.not. state%mp_comm_external) call finish_mp

    !state%gs2_initialized = .false.

    state%init%level = 0

    if (proc0 .and. (.not. quiet_value) ) write(*, '("Run finished at ",A)') date_iso8601()

    !if (.not. present(mpi_comm) .and. .not. nofin) call finish_mp
  end subroutine finalize_gs2

  !> FIXME : Add documentation  
  subroutine prepare_optimisations_overrides(state)
    use overrides, only: init_optimisations_overrides
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init_optimisations_overrides(state%init%opt_ov)
  end subroutine prepare_optimisations_overrides

  !> FIXME : Add documentation  
  subroutine prepare_miller_geometry_overrides(state)
    use overrides, only: init_miller_geometry_overrides
    use gs2_init, only: init, init_level_list
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_miller_geometry-1)
    call init_miller_geometry_overrides(state%init%mgeo_ov)
  end subroutine prepare_miller_geometry_overrides

  !> FIXME : Add documentation  
  subroutine prepare_profiles_overrides(state)
    use overrides, only: init_profiles_overrides
    use gs2_init, only: init, init_level_list
    use species, only: nspec
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_profiles-1)
    call init_profiles_overrides(state%init%prof_ov, nspec)
  end subroutine prepare_profiles_overrides

  !> FIXME : Add documentation  
  subroutine prepare_kt_grids_overrides(state)
    use overrides, only: init_kt_grids_overrides
    use gs2_init, only: init, init_level_list
    type(gs2_program_state_type), intent(inout) :: state

    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_kt_grids-1)
    call init_kt_grids_overrides(state%init%kt_ov)
  end subroutine prepare_kt_grids_overrides

  !> FIXME : Add documentation  
  subroutine prepare_initial_values_overrides(state)
    use overrides, only: init_initial_values_overrides
    use fields, only: force_maxwell_reinit
    use gs2_init, only: in_memory
    use gs2_init, only: init, init_level_list
    use gs2_layouts, only: g_lo
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gexp_3
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. state%included) return
    ! Initialize to the level below so that overrides are triggered
    call init(state%init, init_level_list%override_initial_values-1)
    call init_initial_values_overrides(state%init%initval_ov,&
      ntgrid, ntheta0, naky, g_lo%llim_proc, g_lo%ulim_alloc, &
      force_maxwell_reinit=force_maxwell_reinit, &
      in_memory=in_memory, has_explicit_terms = allocated(gexp_3))
  end subroutine prepare_initial_values_overrides

  !> FIXME : Add documentation  
  subroutine set_initval_overrides_to_current_vals(initval_ov)
    use dist_fn_arrays, only: gnew, gexp_1, gexp_2, gexp_3
    use gs2_save, only: gs2_save_for_restart
    use mp, only: proc0, broadcast, mp_abort
    use collisions, only: vnmult
    use file_utils, only: error_unit
    use run_parameters, only: fphi, fapar, fbpar
    use fields, only:  force_maxwell_reinit
    use fields_arrays, only: phinew, aparnew, bparnew
    use gs2_time, only: user_time, user_dt
    use overrides, only: initial_values_overrides_type
    use antenna, only: dump_ant_amp
    implicit none
    type(initial_values_overrides_type), intent(inout) :: initval_ov
    integer :: istatus
    !if (.not. state%included) return

    if (.not.initval_ov%is_initialised()) &
      call mp_abort("Trying to set initial value overrides &
      & before they are initialized... have you called &
      & prepare_initial_values_overrides ? ", .true.)

    if(initval_ov%in_memory)then
      initval_ov%g=gnew
      initval_ov%vnmult = vnmult
      if (allocated(initval_ov%gexp_1)) initval_ov%gexp_1 = gexp_1
      if (allocated(initval_ov%gexp_2)) initval_ov%gexp_2 = gexp_2
      if (allocated(initval_ov%gexp_3)) initval_ov%gexp_3 = gexp_3

      ! Do not use the value from [old_iface_]state%init%initval_ov%force_maxwell_reinit = initval_ov%force_maxwell_reinit as
      ! follows:
      !     if (.not. initval_ov%force_maxwell_reinit) then
      ! as was done previously. This is because this causes problems in gs2_init so we use the value from
      ! fields::force_maxwell_reinit directly (see comments in gs2_init for details). To use here the value from
      ! initval_ov%force_maxwell_reinit could result in inconsistent values being used in different parts of the code. Therefore, we
      ! now use the value from fields::force_maxwell_reinit directly here too.  Now that the reference to
      ! initval_ov%force_maxwell_reinit here has also been removed, this component is no longer referenced anywhere in the codebase.
      ! Therefore, the force_maxwell_reinit component has been completely removed from the initial_values_overrides_type derived
      ! type to avoid confusion in the future.
      if (.not. force_maxwell_reinit) then
        if(fphi.gt.0) then
          initval_ov%phi=phinew
        end if
        if(fapar.gt.0) then 
          initval_ov%apar=aparnew
        end if
        if(fbpar.gt.0) then 
          initval_ov%bpar=bparnew
        end if
      end if
    else ! if(.not.in_memory)then
      !Should really do this with in_memory=.true. as well but
      !not sure that we really need to as we never read in the dumped data.
      if (proc0) call dump_ant_amp

      call gs2_save_for_restart (gnew, user_time, user_dt, vnmult, istatus, fphi, fapar, fbpar)
    endif

  end subroutine set_initval_overrides_to_current_vals

  !> FIXME : Add documentation  
  subroutine finalize_overrides(state)
    use overrides, only: finish_profiles_overrides
    use overrides, only: finish_miller_geometry_overrides
    use overrides, only: finish_initial_values_overrides
    type(gs2_program_state_type), intent(inout) :: state

    if (state%init%mgeo_ov%is_initialised()) call state%init%mgeo_ov%finish()
    if (state%init%prof_ov%is_initialised()) call state%init%prof_ov%finish()
    if (state%init%initval_ov%is_initialised()) call state%init%initval_ov%finish()
  end subroutine finalize_overrides

  !> FIXME : Add documentation
  subroutine calculate_outputs(state)
    use gs2_diagnostics, only: start_time
    use gs2_diagnostics, only: pflux_avg, qflux_avg, heat_avg, vflux_avg !, start_time, nwrite, write_nl_flux
    use gs2_diagnostics, only: diffusivity
    use gs2_diagnostics, only: ensemble_average
#ifdef NEW_DIAG
    use gs2_diagnostics_new, only: gnostics
#endif 
    use gs2_time, only: user_time
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use run_parameters, only: trinity_linear_fluxes, trinity_ql_fluxes
    use species, only: nspec, spec
    use run_parameters, only: use_old_diagnostics
    use mp, only: mp_abort
    implicit none
    type(gs2_program_state_type), intent(inout) :: state
    real :: time_interval
    real :: diff
    real, dimension(nspec) :: qf, pf, ht, vf
    integer :: is

    if (.not. state%included) return

    !write (*,*) 'GETTING FLUXES', 'user_time', user_time, start_time, 'DIFF', time_interval

    if (state%nensembles > 1) &
      call ensemble_average (state%nensembles, time_interval)

      
   ! Use simple gamma / k^2 estimates for transport in the 
   ! linear case. This is for testing Trinity
   if (trinity_linear_fluxes .and. .not. trinity_ql_fluxes .and. &
         nonlinear_mode_switch .eq. nonlinear_mode_none) then
     if (use_old_diagnostics) then
       diff = diffusivity()
     else
#ifdef NEW_DIAG
       diff = gnostics%current_results%diffusivity
#endif 
     endif
     !write(*,*) 'difff is', diff
     do is = 1,nspec
       ! Q = n chi grad T = n (gamma / k^2) dT / dr
       ! = dens  n_r (gamma_N v_thr / k_N**2 rho_r a) dT / drho drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) T a / L_T drho/dr
       ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) temp T_r tprim drho/dr_N/a
       ! Q / (n_r  v_r T_r rho_r**2/a**2) 
       ! = dens (gamma_N / k_N**2) temp tprim grho
       !   
       ! grho factored to diffusivity in diagnostics
       qf(is) =  diff * spec(is)%dens * spec(is)%temp * spec(is)%tprim 
       pf(is) =  diff * spec(is)%dens**2.0 * spec(is)%fprim 
       ht = 0.0
       vf = 0.0
     end do
       !time_interval = user_time-start_time
     time_interval = 1.0
   else 
     if (use_old_diagnostics) then
       time_interval = user_time-start_time
       qf = qflux_avg
       pf = pflux_avg
       ht = heat_avg
       vf = vflux_avg
     else
#ifdef NEW_DIAG
       time_interval = user_time - gnostics%start_time
       qf = gnostics%current_results%species_heat_flux_avg
       pf = gnostics%current_results%species_particle_flux_avg
       ht = gnostics%current_results%species_heating_avg
       vf = gnostics%current_results%species_momentum_flux_avg
#endif
     endif
   end if

   state%outputs%pflux = pf/time_interval
   state%outputs%qflux = qf/time_interval
   state%outputs%heat = ht/time_interval
   state%outputs%vflux = vf(1)/time_interval
   
   if (trinity_ql_fluxes) then
#ifdef NEW_DIAG
     qf = gnostics%current_results%species_heat_flux
     pf = gnostics%current_results%species_particle_flux
     ht = gnostics%current_results%species_heating
     vf = gnostics%current_results%species_momentum_flux
     state%outputs%pflux = pf/gnostics%current_results%phi2
     state%outputs%qflux = qf/gnostics%current_results%phi2
     state%outputs%heat = ht/gnostics%current_results%phi2
     state%outputs%vflux = vf(1)/gnostics%current_results%phi2
#else
     call mp_abort('ERROR: trinity_ql_fluxes only works with new diagnostics', .true.)
#endif 
   end if
  end subroutine calculate_outputs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> FIXME : Add documentation
  subroutine reset_timers(timers)
    type(gs2_timers_type), intent(inout) :: timers
    timers%init = 0.
    timers%advance = 0.
    timers%timestep = 0.
    timers%finish = 0.
    timers%total = 0. 
    timers%diagnostics=0.
    !timers%interval
    timers%main_loop = 0.
  end subroutine reset_timers

  !> Print timing information to standard output.
  !>
  !> Prints a full breakdown of timing information if `state%print_full_timers`
  !> is true, otherwise just prints the total time.
  subroutine print_times(state, timers)
    use mp, only: proc0, min_reduce, max_reduce, sum_reduce, job
    use mp, only: get_mp_times
    use redistribute, only: get_redist_times
    use fft_work, only: time_fft
    use fields_arrays, only: time_field, time_field_mpi
    use gs2_reinit, only: time_reinit
    use le_derivatives, only: time_vspace_derivatives, time_vspace_derivatives_mpi
    use nonlinear_terms, only: time_add_explicit_terms, time_add_explicit_terms_mpi
#ifdef GS2_EXPERIMENTAL_ENABLE_MEMORY_USAGE
    use memory_usage, only: print_memory_usage
#endif
    implicit none
    !> Program state information, needed to print
    type(gs2_program_state_type), intent(in) :: state
    !> Timing information
    type(gs2_timers_type), intent(in) :: timers
    real :: mp_total, mp_over, mp_coll, mp_ptp, mp_sync
    real :: redist_total, redist_mpi, redist_copy
    if (.not. proc0) return

    if (state%print_full_timers) then
      call get_mp_times(mp_total, mp_over, mp_coll, mp_ptp, mp_sync)
      call get_redist_times(redist_total, redist_mpi, redist_copy)

      !Blank line to force separation from other outputs
      write(*, '("")')
      write(*, '(" ", a)') formatted_time("Initialization", timers%init(1))
      write(*, '(" ", a)') formatted_time("Advance steps", timers%advance(1))
#ifdef WITH_EIG
      write(*, '(" ", a)') formatted_time("Eigensolver", timers%eigval(1))
#endif
      write(*, '("(", a, ")")') formatted_time("field solve", time_field(1))
      write(*, '("(", a, ")")') formatted_time("field solve mpi", time_field_mpi)
      write(*, '("(", a, ")")') formatted_time("collisions", time_vspace_derivatives(1))
      write(*, '("(", a, ")")') formatted_time("collisions mpi", time_vspace_derivatives_mpi)
      write(*, '("(", a, ")")') formatted_time("explicit/nl", time_add_explicit_terms(1))
      write(*, '("(", a, ")")') formatted_time("explicit/nl mpi", time_add_explicit_terms_mpi)
      write(*, '("(", a, ")")') formatted_time("explicit/nl fft", time_fft(1))
      write(*, '("(", a, ")")') formatted_time("diagnostics", timers%diagnostics(1))
      write(*, '(" ", a)') formatted_time("Re-initialize", time_reinit(1))
      write(*, '(" ", a)') formatted_time("Redistribute", redist_total)
      write(*, '("(", a, ")")') formatted_time("redistribute mpi", redist_mpi)
      write(*, '("(", a, ")")') formatted_time("redistribute copy", redist_copy)
      write(*, '(" ", a)') formatted_time("Finishing", timers%finish(1))
      write(*, '(" ", a)') formatted_time("MPI" ,mp_total)
      write(*, '("(", a, ")")') formatted_time("Overheads", mp_over)
      write(*, '("(", a, ")")') formatted_time("Collectives", mp_coll)
      write(*, '("(", a, ")")') formatted_time("PTP", mp_ptp)
      write(*, '("(", a, ")")') formatted_time("Sync", mp_sync)
      write(*, '(" total from timer is:", 0pf9.2, " min", /)') timers%total(1) / 60.
    else
      if (state%external_job_id >= 0) then
        write(*, '(a, i0, ", ")', advance="no") " External job ID: ", state%external_job_id
      end if
      if (state%list) then
        write(*, '(a, i0, ", ")', advance="no") " Job ID: ", job
      end if

      write(*, '("Total from timer is:", 0pf9.2," min")') state%timers%total(1)/60.
    endif

#ifdef GS2_EXPERIMENTAL_ENABLE_MEMORY_USAGE
    call print_memory_usage
#endif
  contains
    !> Print a single formatted timer line
    character(len=46) function formatted_time(name, timer)
      !> Name of the timer
      character(len=*), intent(in) :: name
      !> Time in seconds
      real, intent(in) :: timer

      write(formatted_time, "(a, T25, 0pf9.3, ' min', T40, 2pf5.1, ' %')") name, timer / 60., timer/timers%total(1)
    end function formatted_time
  end subroutine print_times

  !> FIXME : Add documentation  
  subroutine allocate_outputs(state)
    use species, only: nspec
    type(gs2_program_state_type), intent(inout) :: state
    if (.not. allocated(state%outputs%pflux)) then
      allocate(state%outputs%pflux(nspec))
      allocate(state%outputs%qflux(nspec))
      allocate(state%outputs%heat(nspec))
    end if
  end subroutine allocate_outputs

  !> FIXME : Add documentation  
  subroutine deallocate_outputs(state)
    type(gs2_program_state_type), intent(inout) :: state
    if (allocated(state%outputs%pflux)) then
      deallocate(state%outputs%pflux)
      deallocate(state%outputs%qflux)
      deallocate(state%outputs%heat)
    end if
  end subroutine deallocate_outputs

  !> FIXME : Add documentation  
  subroutine reset_linear_magnitude
    use dist_fn_arrays, only: g, gnew
    use fields, only: set_init_fields
    use fields_arrays, only: phi, apar, bpar, phinew, aparnew,  bparnew
    use run_parameters, only: fphi, fapar, fbpar
    use init_g, only: ginit
    implicit none
    logical :: dummy
    
    if (fphi > 0.) then 
       phi = 0.0 
       phinew = 0.0
    end if
    
    if (fapar > 0.) then 
       apar = 0.0
       aparnew = 0.0
    end if
    if (fbpar > 0.) then 
       bpar = 0.0
       bparnew = 0.0
    end if
    g = 0.0; gnew = 0.0
    
    call ginit(dummy)
    call set_init_fields
  end subroutine reset_linear_magnitude

  !> FIXME : Add documentation  
  subroutine write_trinity_parameters(state)
      use file_utils, only: open_output_file, close_output_file
      use theta_grid_params, only: write_theta_grid => write_trinity_parameters
      use species, only: write_species => write_trinity_parameters
      use run_parameters, only: write_run_parameters => write_trinity_parameters
      use mp, only: proc0
      implicit none
      type(gs2_program_state_type), intent(in) :: state
      integer :: trinpars_unit
      character(len=1000) :: ext

      if (state%trinity_timestep .gt. 0) then
        write(ext,'(".",I0,".",I0,".trinpars")') state%trinity_timestep, &
          state%trinity_iteration
      else
        ext = '.trinpars'
      end if
      
      if (proc0) then
        call open_output_file(trinpars_unit, trim(ext))
        call write_theta_grid(trinpars_unit)
        call write_species(trinpars_unit)
        call write_run_parameters(trinpars_unit)
        call close_output_file(trinpars_unit)
      end if
  end subroutine write_trinity_parameters

  !> Write an input file containing the current state of all input
  !> parameters
  subroutine write_used_inputs_file(header)
    use mp, only: proc0
    use file_utils, only: open_output_file, close_output_file, run_name
    use standard_header, only: standard_header_type
    use config_collection, only : gs2_config_type

    !> Header for files with build and run information
    type(standard_header_type), intent(in), optional :: header

    ! Actual value for optional `header` input
    type(standard_header_type) :: local_header
    ! File unit to write to
    integer :: used_inputs_unit
    type(gs2_config_type) :: input_config

    ! Only write this file on rank 0
    if (.not. proc0) return

    if (present(header)) then
      local_header = header
    else
      local_header = standard_header_type()
    end if

    call input_config%get_configs()
    call open_output_file(used_inputs_unit, ".used_inputs.in")
    write(unit=used_inputs_unit, fmt='(a)') &
         local_header%to_string( &
         comment_character="! ", &
         file_description="Input parameters as used by " // trim(run_name) &
         )

    call input_config%write_to_unit(used_inputs_unit)
    call close_output_file(used_inputs_unit)

  end subroutine write_used_inputs_file

end module gs2_main
