!> This module is basically a store for the input parameters that are
!> specified in the namelists knobs and parameters.
!> In general, the names of the public variables in this module are
!> the same as the name of the input parameter they correspond to.
module run_parameters
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: init_run_parameters, finish_run_parameters, check_run_parameters
  public :: wnml_run_parameters
  public :: beta, zeff, tite, reset, immediate_reset, delt, ieqzip, k0, user_comments
  public :: fphi, fapar, fbpar, has_phi, has_apar, has_bpar, progress_frequency
  public :: nstep, max_sim_time, wstar_units, rhostar, neo_test, ncheck_stop
  public :: delt_option_switch, delt_option_auto, avail_cpu_time, margin_cpu_time
  public :: do_eigsolve, save_timer_statistics, save_init_times
  !> If true use old diagnostics. Included for testing
  !> only and will eventually be removed. Please use
  !> the new diagnostics module!
  public :: use_old_diagnostics
  public :: set_overrides, get_knobs_config, knobs_config_type, set_run_parameters_config
  
  real :: beta, zeff, tite, fphi, fapar, fbpar, k0, rhostar
  logical :: has_phi, has_apar, has_bpar, neo_test
  real :: delt, avail_cpu_time, margin_cpu_time, max_sim_time
  integer :: ncheck_stop, nstep, seed, progress_frequency
  logical :: reset = .false., initialized = .false.
  logical :: immediate_reset, wstar_units, save_timer_statistics, save_init_times
  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: knexist, do_eigsolve, use_old_diagnostics
  character(len=100000) :: user_comments

  !> Used to indicate which wavenumbers we do not want to evolve.
  !> This choice is controlled by the [[eqzip_option]]
  !> input. Primarily intended for secondary mode analysis.
  logical, allocatable :: ieqzip(:,:)
  integer :: eqzip_option_switch
  integer, parameter :: &
       eqzip_option_none = 1, &
       eqzip_option_secondary = 2, &
       eqzip_option_tertiary = 3, &
       eqzip_option_equilibrium = 4

  !> Used to represent the input configuration of run_parameters through "knobs"
  type, extends(abstract_config_type) :: knobs_config_type
     ! namelist : knobs
     ! indexed : false     
     !> The maximum wall clock time available to GS2
     real :: avail_cpu_time = 1.0e10
     !> \(\beta\) is the ratio of the reference pressure to the reference
     !> magnetic energy density, \(\beta=2\mu_0n_{ref}T_{ref}/B^2_{ref}\).
     real :: beta = 0.0
     !> Initial timestep
     real :: delt = 0.1
     !> Determines how initial timestep is set. Possible options are:
     !>
     !>  - "default" or "set_by_hand": use timestep from input file
     !>  - "check_restart": read timestep(s) from restart file
     !>
     !> If "check_restart" is used but the restart files aren't
     !> present for any reason (for instance, the current run is not a
     !> restart), GS2 will error when trying to read the
     !> timesteps. You can run the [[ingen]] tool on your input file
     !> to check this issue before running your job
     character(len = 20) :: delt_option = 'default'
     !> If true then use eigensolver instead of initial value solver.
     logical :: do_eigsolve = .false.
     !> Advanced option to freeze evolution of certain modes. Possible values are:
     !>
     !> - "secondary": don't evolve `ik == 2, it == 1`;
     !> - "tertiary": dont' evolve `ik == 1, it == 2 or ntheta0`;
     !> - "harris": don't evolve `ik == 1`.
     !>
     !> Only used in [[dist_fn]], should be moved to [[dist_fn_knobs]].
     character(len = 20) :: eqzip_option = 'none'
     !> Multiplies \(A_\|\) throughout.
     real :: fapar = 0.0
     !> Multiplies \(B_\|\) throughout.
     real :: fbpar = -1.0
     !> Multiplies \(\phi\) throughout.
     real :: fphi = 1.0
     !> Determines the behaviour when the CFL condition is broken in the nonlinear term:
     logical :: immediate_reset = .true.
     !> Used in [[init_g_knobs:ginit_option]] "harris" and "convect". Should be
     !> moved to init_g namelist.
     real :: k0 = 1.0
     !> How close to avail_cpu_time can we go before trying to stop the run cleanly
     real :: margin_cpu_time = 300.0
     !> The simulation time after which the run should stop.
     real :: max_sim_time = 1.0e6
     !> Sets the time step interval at which check for the existence
     !> of the <runname>.stop file.
     integer :: ncheck_stop = 5
     !> Lowflow only. If true, GS2 exits after writing out neoclassical
     !> quantities of interest
     logical :: neo_test = .false.
     !> Maximum number of steps to take
     integer :: nstep = 100
     !> How frequently, in integer steps, do we display a progress message.
     !> If <= 0 then no progress messages are displayed
     integer :: progress_frequency = 0
     !> Normalised gyro-radius, only used for low flow builds
     real :: rhostar = 3.0e-3
     !> If true then save init timer information to
     !> <run_name>.init_times.
     logical :: save_init_times = .false.
     !> If true then save some extra timer information to
     !> <run_name>.timing_stats giving min/max/mean values.
     logical :: save_timer_statistics = .true.
     !> Used to set the random seed used in the random number
     !> generators provided by [[ran]]. If this input is 0 (the
     !> default) then we don't set the seed and use the system default
     !> instead.
     integer :: seed = 0
     !> Only used when there is an adiabatic species. Sets `n q ^2 / T` of the adiabatic
     !> species, normalised to the reference values.
     real :: tite = 1.0
     !> If true use original diagnostics rather than new form
     logical :: use_old_diagnostics = .false.
     !> Custom description of run to be added to netcdf output (new diagnostics only)
     character(len = 65000) :: user_comments = ''
     !> If true makes timestep proportional to ky*rho. Only sensible for linear runs.
     logical :: wstar_units = .false.
     !> Effective ionic charge appearing in electron collision frequency
     real :: zeff = 1.0
   contains
     procedure, public :: read => read_knobs_config
     procedure, public :: write => write_knobs_config
     procedure, public :: reset => reset_knobs_config
     procedure, public :: broadcast => broadcast_knobs_config
     procedure, public, nopass :: get_default_name => get_default_name_knobs_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_knobs_config
  end type knobs_config_type
  
  type(knobs_config_type) :: knobs_config

contains

  !> FIXME : Add documentation  
  subroutine check_run_parameters(report_unit)
    use species, only: has_hybrid_electron_species, spec
    implicit none
    integer, intent(in) :: report_unit
    if (fphi /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fphi in the knobs namelist = ',e11.4)") fphi
       write (report_unit, fmt="('fphi is a scale factor of all instances of Phi (the electrostatic potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fapar == 0.) then
       write (report_unit, fmt="('A_parallel will not be included in the calculation.')")
    end if
    if (fapar == 1.) then
       write (report_unit, fmt="('A_parallel will be included in the calculation.')")
    end if
    if (fapar /= 0. .and. fapar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fapar in the knobs namelist = ',e11.4)") fapar
       write (report_unit, fmt="('fapar is a scale factor of all instances of A_parallel (the parallel vector potential).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (fbpar == 0.) then
       write (report_unit, fmt="('B_parallel will not be included in the calculation.')")
    end if
    if (fbpar == 1.) then
       write (report_unit, fmt="('B_parallel will be included in the calculation.')")
    end if
    if (fbpar /= 0. .and. fbpar /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('fbpar in the knobs namelist = ',e11.4)") fbpar
       write (report_unit, fmt="('fbpar is a scale factor of all instances of B_parallel &
           & (the perturbed parallel magnetic field).')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    write (report_unit, *) 
    if(immediate_reset)then
       write (report_unit, fmt="('The time step will be reset immediately after cfl violation detected.')") 
    else
       write (report_unit, fmt="('The time step will be reset just before the next time step after cfl violation detected.')") 
    endif

    write (report_unit, *)
    if(seed .ne. 0)then
       write (report_unit, fmt="('The random number generator will use a seed derived from use input, ',I0,'.')") seed
    else
       write (report_unit, fmt="('The random number generator will use the default seed.')")
    endif
    write (report_unit, *)

    if ( has_hybrid_electron_species(spec) .and. beta /= 0.0 &
         .and. (fapar /= 0 .or. fbpar /=0) ) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You are using a hybrid electron species in an electromagnetic simulation.')")
       write (report_unit, fmt="('This is probably not physical.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    end if

  end subroutine check_run_parameters

  !> FIXME : Add documentation    
  subroutine wnml_run_parameters(unit,electrons,collisions)
    implicit none
    integer, intent(in) :: unit
    logical, intent(in) :: electrons, collisions
    if (knexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "knobs"
       write (unit, fmt="(' fphi   = ',f6.3)") fphi
       write (unit, fmt="(' fapar  = ',f6.3)") fapar
       write (unit, fmt="(' fbpar = ',f6.3)") fbpar
       write (unit, fmt="(' delt = ',e17.10)") delt
       write (unit, fmt="(' max_sim_time = ',e17.10)") max_sim_time
       write (unit, fmt="(' nstep = ',i8)") nstep
       write (unit, fmt="(' wstar_units = ',L1)") wstar_units
       select case (delt_option_switch)
       case (delt_option_auto)
          write (unit, fmt="(' delt_option = ',a)") '"check_restart"'
       case (delt_option_hand)
          ! nothing
       end select
       write (unit, fmt="(' immediate_reset = ',L1)") immediate_reset
       write (unit, fmt="(' beta = ',e17.10)") beta       ! if zero, fapar, fbpar should be zero
       if (collisions) write (unit, fmt="(' zeff = ',e17.10)") zeff
       if (.not. electrons)  write (unit, fmt="(' tite = ',e17.10)") tite
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_run_parameters

  !> FIXME : Add documentation  
  subroutine init_run_parameters(knobs_config_in)
    use kt_grids, only: init_kt_grids, naky, nakx => ntheta0
    implicit none
    type(knobs_config_type), intent(in), optional :: knobs_config_in

    if (initialized) return
    initialized = .true.

    call read_parameters(knobs_config_in)

    call init_kt_grids

    if(.not.allocated(ieqzip)) allocate(ieqzip(nakx,naky))
    ieqzip = .false.
    select case (eqzip_option_switch)
    case (eqzip_option_secondary)
       ! suppress evolution of secondary mode
       ieqzip(1,2) = .true.
    case (eqzip_option_tertiary)
       ! suppress evolution of tertiary mode
       ieqzip(2,1) = .true.
       ieqzip(nakx,1) = .true.
    case (eqzip_option_equilibrium)
       ! suppress evolution of 1D equilibrium (x dependent)
       ieqzip(1:nakx,1) = .true.
    end select
  end subroutine init_run_parameters

  !> FIXME : Add documentation    
  subroutine read_parameters(knobs_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast, mp_abort
    use text_options, only: text_option, get_option_value
    use ran, only: set_seed_from_single_integer
    implicit none
    type(knobs_config_type), intent(in), optional :: knobs_config_in
    
    type (text_option), dimension (4), parameter :: eqzipopts = &
         (/ text_option('none', eqzip_option_none), &
            text_option('secondary', eqzip_option_secondary), &
            text_option('tertiary', eqzip_option_tertiary), &
            text_option('equilibrium', eqzip_option_equilibrium) /)
    character (len=20) :: eqzip_option
    type (text_option), dimension (3), parameter :: deltopts = &
         (/ text_option('default', delt_option_hand), &
            text_option('set_by_hand', delt_option_hand), &
            text_option('check_restart', delt_option_auto) /)
    character(20) :: delt_option
    integer :: ierr, in_file
    logical :: rpexist
    !> Identify if we need to warn about removal of parameters
    rpexist = .false.
    if (proc0) in_file = input_unit_exist('parameters', rpexist)
    call broadcast(rpexist)
    if (rpexist) call mp_abort('Input file contains deprecated "parameters" namelist -- Please move inputs to "knobs".', .true.)
    
    if (present(knobs_config_in)) knobs_config = knobs_config_in
    call knobs_config%init(name = 'knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    avail_cpu_time = knobs_config%avail_cpu_time
    beta = knobs_config%beta
    delt = knobs_config%delt
    delt_option = knobs_config%delt_option
    do_eigsolve = knobs_config%do_eigsolve
    eqzip_option = knobs_config%eqzip_option
    fapar = knobs_config%fapar
    fbpar = knobs_config%fbpar
    fphi = knobs_config%fphi
    immediate_reset = knobs_config%immediate_reset
    k0 = knobs_config%k0
    margin_cpu_time = knobs_config%margin_cpu_time
    max_sim_time = knobs_config%max_sim_time
    ncheck_stop = knobs_config%ncheck_stop
    neo_test = knobs_config%neo_test
    nstep = knobs_config%nstep
    progress_frequency = knobs_config%progress_frequency
    rhostar = knobs_config%rhostar
    save_init_times = knobs_config%save_init_times
    save_timer_statistics = knobs_config%save_timer_statistics
    seed = knobs_config%seed
    tite = knobs_config%tite
    use_old_diagnostics = knobs_config%use_old_diagnostics
    user_comments = knobs_config%user_comments
    wstar_units = knobs_config%wstar_units
    zeff = knobs_config%zeff
  
    knexist = knobs_config%exist

    ! Now handle inputs to calculate derived quantities

    ! Override fapar, fbpar and beta if set inconsistently:
    
    !> If we have zero beta then disable apar and bpar fields
    !! as these contribute nothing to result but slow down calculation.
    if(beta.eq.0) then
       if(((fapar.ne.0) .or. (fbpar.ne.0)).and. proc0) then
          ierr = error_unit()
          write(ierr,'("Warning: Disabling apar and bpar as beta = 0.")')
       endif
       fapar = 0.
       fbpar = 0.
    endif

    has_phi = fphi > epsilon(0.0)
    has_apar = fapar > epsilon(0.0)
    has_bpar = fbpar > epsilon(0.0)

    !> If we have non-zero beta then alert the user that they have some
    !! of the perturbed magnetic fields disabled. This may be intended
    !! behaviour, so we throw a warning. However, if both are disabled,
    !! we set beta=0. This will make the simulation electrostatic, as 
    !! probably intended if both apar and bpar are set to zero.
    if( beta /= 0.0 .and. .not. (has_apar .or. has_bpar)) then
       if(proc0) then
          ierr = error_unit()
          write(ierr,'("Warning: Both fapar and fbpar are zero: setting beta = 0.")')
       endif
       beta = 0.0
    endif
    if((beta.ne.0).and.proc0) then
       ierr = error_unit()
       if(.not. has_apar) write(ierr,'("Warning: Running with finite beta but fapar=0.")')
       if(.not. has_bpar) write(ierr,'("Warning: Running with finite beta but fbpar=0.")')
    endif

    ierr = error_unit()
    call get_option_value &
         (delt_option, deltopts, delt_option_switch, ierr, &
         "delt_option in knobs",.true.)
    
    call get_option_value ( &
         eqzip_option, eqzipopts, eqzip_option_switch, error_unit(), &
         "eqzip_option in knobs",.true.)

    if (seed .ne. 0 ) call set_seed_from_single_integer(seed)
  end subroutine read_parameters

  !> FIXME : Add documentation    
  subroutine set_overrides(tstep_ov)
    use overrides, only: timestep_overrides_type
    type(timestep_overrides_type), intent(in) :: tstep_ov
    if (tstep_ov%override_immediate_reset) immediate_reset=tstep_ov%immediate_reset
  end subroutine set_overrides

  !> FIXME : Add documentation    
  subroutine finish_run_parameters
    implicit none

    if (allocated(ieqzip)) deallocate (ieqzip)

    initialized = .false.
    call knobs_config%reset()
  end subroutine finish_run_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_run_parameters_config(knobs_config_in)
    use mp, only: mp_abort
    type(knobs_config_type), intent(in), optional :: knobs_config_in
    
    if (initialized) then
       call mp_abort("Trying to set run_parameters config when already initialized.", to_screen = .true.)
    end if
    if (present(knobs_config_in)) then
       knobs_config = knobs_config_in
    end if
  end subroutine set_run_parameters_config

#include "knobs_auto_gen.inc"  
end module run_parameters
