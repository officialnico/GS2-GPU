!> This module is basically a store for the input parameters that are
!> specified in the namelists knobs and parameters.
!> In general, the names of the public variables in this module are
!> the same as the name of the input parameter they correspond to.
module run_parameters
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: init_run_parameters, finish_run_parameters
  public :: check_run_parameters, wnml_run_parameters
  public :: write_trinity_parameters
  public :: beta, zeff, tite, reset, immediate_reset, delt
  public :: fphi, fapar, fbpar, has_phi, has_apar, has_bpar
  public :: nstep, max_sim_time, wstar_units, margin
  public :: ieqzip, k0
  public :: delt_option_switch, delt_option_hand
  public :: delt_option_auto
  public :: avail_cpu_time, margin_cpu_time
  public :: rhostar, neo_test
  public :: do_eigsolve, save_timer_statistics, save_init_times
  !> If true and nonlinear_mode is "off", return
  !! simple diffusive estimates of fluxes to trinity
  public :: trinity_linear_fluxes
  !> If true and nonlinear_mode is "off", 
  !! and trinity_linear_fluxes is true, return
  !! quasilinear estimates of fluxes to trinity
  public :: trinity_ql_fluxes
  public :: user_comments
  
  !> If true use old diagnostics. Included for testing
  !! only and will eventually be removed. Please use 
  !! the new diagnostics module!
  public :: use_old_diagnostics

  public :: set_overrides

  public :: knobs_config_type
  public :: parameters_config_type
  public :: set_run_parameters_config
  public :: get_knobs_config
  public :: get_parameters_config
  
  real :: beta, zeff, tite
  real :: fphi, fapar, fbpar, faperp
  logical :: has_phi, has_apar, has_bpar
  real :: margin, delt
  real :: avail_cpu_time, margin_cpu_time
  real :: max_sim_time
  integer :: nstep, seed
  logical :: reset=.false.
  logical :: immediate_reset
  logical :: wstar_units
  logical :: save_timer_statistics, save_init_times
  real :: k0
  integer :: delt_option_switch
  integer, parameter :: delt_option_hand = 1, delt_option_auto = 2
  logical :: initialized = .false.
  logical :: rpexist, knexist
  real :: rhostar
  logical :: neo_test
  logical :: trinity_linear_fluxes, do_eigsolve
  logical :: trinity_ql_fluxes
  logical :: use_old_diagnostics

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
     !> magnetic energy density, \(\beta=2\mu_0n_{ref}T_{ref}/B^2_{ref}\).  This
     !> gets a smart default from parameters
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
     !> Deprecated, synonym for `fapar`
     real :: faperp = 0.0
     !> Multiplies \(B_\|\) throughout.
     real :: fbpar = -1.0
     !> Multiplies \(\phi\) throughout.
     real :: fphi = 1.0
     !> Determines the behaviour when the CFL condition is broken in the nonlinear term:
     logical :: immediate_reset = .true.
     !> Used in [[init_g_knobs:ginit_option]] "harris" and "convect". Should be
     !> moved to init_g namelist. This gets a smart default from [[parameters]].
     real :: k0 = 1.0
     !> Not used anywhere -- to be removed.
     real :: margin = 0.05
     !> How close to avail_cpu_time can we go before trying to stop the run cleanly
     real :: margin_cpu_time = 300.0
     !> The simulation time after which the run should stop.
     real :: max_sim_time = 1.0e6
     !> Lowflow only. If true, GS2 exits after writing out neoclassical
     !> quantities of interest
     logical :: neo_test = .false.
     !> Maximum number of steps to take
     integer :: nstep = 100
     !> Normalised gyro-radius, only used for low flow builds
     !> This gets a smart default from parameters
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
     !> Deprecated
     !> This gets a smart default from parameters
     real :: teti = -100.0
     !> Only used when there is an adiabatic species. Sets `n q ^2 / T` of the adiabatic
     !> species, normalised to the reference values.
     !> This gets a default from the parameters namelist.
     real :: tite = 1.0
     !> If true and running linearly, return linear diffusive flux estimates to Trinity.
     logical :: trinity_linear_fluxes = .false.
     !> Unknown
     logical :: trinity_ql_fluxes = .false.
     !> If true use original diagnostics rather than new form
     logical :: use_old_diagnostics = .false.
     !> Custom description of run to be added to netcdf output (new diagnostics only)
     !> This gets a smart default from parameters
     character(len = 100000) :: user_comments = ''
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

  !> Used to represent the input configuration of run_parameters through "parameters"
  !!
  !! @todo This should probably be merged with knobs. This is actually
  !! implemented here by adding the contents of "parameters" to "knobs".
  !! This should be documented clearly and eventually we should remove parameters
  type, extends(abstract_config_type) :: parameters_config_type
     ! namelist : parameters
     ! indexed : false     
     !> \(\beta\) is the ratio of the reference pressure to the reference magnetic energy density, \(\beta=2\mu_0n_{ref}T_{ref}/B^2_{ref}\).
     real :: beta = 0.0
     !> Unknown, only uses in [[init_g]] -- should be moved to init_g namelist.
     real :: k0 = 1.0
     !> Normalised gyro-radius, only used for low flow builds
     real :: rhostar = 3.0e-3
     !> Deprecated
     real :: teti = -100.0
     !> See [[knobs:tite]].
     real :: tite = 1.0
     !> Custom description of run to be added to netcdf output (new diagnostics only)
     character(len = 100000) :: user_comments = ''
     !> Effective ionic charge appearing in electron collision frequency
     real :: zeff = 1.0
   contains
     procedure, public :: read => read_parameters_config
     procedure, public :: write => write_parameters_config
     procedure, public :: reset => reset_parameters_config
     procedure, public :: broadcast => broadcast_parameters_config
     procedure, public, nopass :: get_default_name => get_default_name_parameters_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_parameters_config
  end type parameters_config_type
  
  type(parameters_config_type) :: parameters_config

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
    if (rpexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "parameters"
       write (unit, fmt="(' beta = ',e17.10)") beta       ! if zero, fapar, fbpar should be zero
       if (collisions) write (unit, fmt="(' zeff = ',e17.10)") zeff
       if (.not. electrons)  write (unit, fmt="(' tite = ',e17.10)") tite
       write (unit, fmt="(' /')")
    endif
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
       write (unit, fmt="(' margin = ',e17.10)") margin
       select case (delt_option_switch)
       case (delt_option_auto)
          write (unit, fmt="(' delt_option = ',a)") '"check_restart"'
       case (delt_option_hand)
          ! nothing
       end select
       write (unit, fmt="(' immediate_reset = ',L1)") immediate_reset
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_run_parameters

  !> FIXME : Add documentation  
  subroutine init_run_parameters(knobs_config_in, parameters_config_in)
    use kt_grids, only: init_kt_grids, naky, nakx => ntheta0
    implicit none
    type(knobs_config_type), intent(in), optional :: knobs_config_in
    type(parameters_config_type), intent(in), optional :: parameters_config_in

    if (initialized) return
    initialized = .true.

    call read_parameters(knobs_config_in, parameters_config_in)

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
  subroutine read_parameters(knobs_config_in, parameters_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use text_options, only: text_option, get_option_value
    use kt_grids, only: gryfx
    use ran, only: set_seed_from_single_integer
    implicit none
    type(knobs_config_type), intent(in), optional :: knobs_config_in
    type(parameters_config_type), intent(in), optional :: parameters_config_in
    
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
    integer :: ierr
    real :: teti  ! for back-compatibility

    if (present(parameters_config_in)) parameters_config = parameters_config_in

    call parameters_config%init(name = 'parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    beta = parameters_config%beta
    k0 = parameters_config%k0
    rhostar = parameters_config%rhostar
    teti = parameters_config%teti
    tite = parameters_config%tite
    user_comments = parameters_config%user_comments
    zeff = parameters_config%zeff

    rpexist = parameters_config%exist

    if (present(knobs_config_in)) knobs_config = knobs_config_in

    ! Set smart defaults
    if (.not.knobs_config%is_initialised()) then    
       knobs_config%beta = beta
       knobs_config%k0 = k0
       knobs_config%rhostar = rhostar
       knobs_config%teti = teti
       knobs_config%tite = tite
       knobs_config%user_comments = user_comments
       knobs_config%zeff = zeff       
    end if
    
    call knobs_config%init(name = 'knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    avail_cpu_time = knobs_config%avail_cpu_time
    beta = knobs_config%beta
    delt = knobs_config%delt
    delt_option = knobs_config%delt_option
    do_eigsolve = knobs_config%do_eigsolve
    eqzip_option = knobs_config%eqzip_option
    fapar = knobs_config%fapar
    faperp = knobs_config%faperp
    fbpar = knobs_config%fbpar
    fphi = knobs_config%fphi
    immediate_reset = knobs_config%immediate_reset
    k0 = knobs_config%k0
    margin = knobs_config%margin
    margin_cpu_time = knobs_config%margin_cpu_time
    max_sim_time = knobs_config%max_sim_time
    neo_test = knobs_config%neo_test
    nstep = knobs_config%nstep
    rhostar = knobs_config%rhostar
    save_init_times = knobs_config%save_init_times
    save_timer_statistics = knobs_config%save_timer_statistics
    seed = knobs_config%seed
    teti = knobs_config%teti
    tite = knobs_config%tite
    trinity_linear_fluxes = knobs_config%trinity_linear_fluxes
    trinity_ql_fluxes = knobs_config%trinity_ql_fluxes
    use_old_diagnostics = knobs_config%use_old_diagnostics
    user_comments = knobs_config%user_comments
    wstar_units = knobs_config%wstar_units
    zeff = knobs_config%zeff
  
    knexist = knobs_config%exist

    ! Now handle inputs to calculate derived quantities
    
    if (teti /= -100.0) tite = teti

    ! Allow faperp-style initialization for backwards compatibility.
    ! Only fbpar is used outside of this subroutine.
    if (fbpar == -1.) then
       fbpar = faperp
    end if
    
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

    ! FOR GRYFX
    if (gryfx()) then
       delt = delt/sqrt(2.0)
       max_sim_time = max_sim_time / sqrt(2.0)
       nstep = nstep*2
    end if

    if (seed .ne. 0 ) call set_seed_from_single_integer(seed)
  end subroutine read_parameters

  !> FIXME : Add documentation    
  subroutine write_trinity_parameters(trinpars_unit)
    integer, intent(in) :: trinpars_unit
    write(trinpars_unit, "(A15)") '&run_parameters'
    write (trinpars_unit, *) ' beta = ', beta
    write (trinpars_unit, "(A1)") '/'

  end subroutine write_trinity_parameters

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
    call parameters_config%reset()    
  end subroutine finish_run_parameters


  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_run_parameters_config(knobs_config_in, parameters_config_in)
    use mp, only: mp_abort
    type(knobs_config_type), intent(in), optional :: knobs_config_in
    type(parameters_config_type), intent(in), optional :: parameters_config_in    
    
    if (initialized) then
       call mp_abort("Trying to set run_parameters config when already initialized.", to_screen = .true.)
    end if
    if (present(knobs_config_in)) then
       knobs_config = knobs_config_in
    end if
    if (present(parameters_config_in)) then
       parameters_config = parameters_config_in
    end if
  end subroutine set_run_parameters_config

  !---------------------------------------
  ! Following is for the knobs config_type
  !---------------------------------------

  !> Reads in the knobs namelist and populates the member variables
  subroutine read_knobs_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(knobs_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    character(len = 100000) :: user_comments
    character(len = 20) :: delt_option, eqzip_option
    integer :: nstep, seed
    logical :: do_eigsolve, immediate_reset, neo_test, save_init_times, save_timer_statistics, &
         trinity_linear_fluxes, trinity_ql_fluxes, use_old_diagnostics, wstar_units
    real :: avail_cpu_time, beta, delt, fapar, faperp, fbpar, fphi, k0, margin, margin_cpu_time, max_sim_time, rhostar, teti, tite, zeff

    namelist /knobs/ avail_cpu_time, beta, delt, delt_option, do_eigsolve, eqzip_option, fapar, faperp, fbpar, fphi, immediate_reset, k0, margin, margin_cpu_time, max_sim_time, neo_test, nstep, rhostar, &
         save_init_times, save_timer_statistics, seed, teti, tite, trinity_linear_fluxes, trinity_ql_fluxes, use_old_diagnostics, user_comments, wstar_units, zeff

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    avail_cpu_time = self%avail_cpu_time
    beta = self%beta
    delt = self%delt
    delt_option = self%delt_option
    do_eigsolve = self%do_eigsolve
    eqzip_option = self%eqzip_option
    fapar = self%fapar
    faperp = self%faperp
    fbpar = self%fbpar
    fphi = self%fphi
    immediate_reset = self%immediate_reset
    k0 = self%k0
    margin = self%margin
    margin_cpu_time = self%margin_cpu_time
    max_sim_time = self%max_sim_time
    neo_test = self%neo_test
    nstep = self%nstep
    rhostar = self%rhostar
    save_init_times = self%save_init_times
    save_timer_statistics = self%save_timer_statistics
    seed = self%seed
    teti = self%teti
    tite = self%tite
    trinity_linear_fluxes = self%trinity_linear_fluxes
    trinity_ql_fluxes = self%trinity_ql_fluxes
    use_old_diagnostics = self%use_old_diagnostics
    user_comments = self%user_comments
    wstar_units = self%wstar_units
    zeff = self%zeff

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = knobs)

    ! Now copy from local variables into type members
    self%avail_cpu_time = avail_cpu_time
    self%beta = beta
    self%delt = delt
    self%delt_option = delt_option
    self%do_eigsolve = do_eigsolve
    self%eqzip_option = eqzip_option
    self%fapar = fapar
    self%faperp = faperp
    self%fbpar = fbpar
    self%fphi = fphi
    self%immediate_reset = immediate_reset
    self%k0 = k0
    self%margin = margin
    self%margin_cpu_time = margin_cpu_time
    self%max_sim_time = max_sim_time
    self%neo_test = neo_test
    self%nstep = nstep
    self%rhostar = rhostar
    self%save_init_times = save_init_times
    self%save_timer_statistics = save_timer_statistics
    self%seed = seed
    self%teti = teti
    self%tite = tite
    self%trinity_linear_fluxes = trinity_linear_fluxes
    self%trinity_ql_fluxes = trinity_ql_fluxes
    self%use_old_diagnostics = use_old_diagnostics
    self%user_comments = user_comments
    self%wstar_units = wstar_units
    self%zeff = zeff

    self%exist = exist
  end subroutine read_knobs_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_knobs_config(self, unit)
    implicit none
    class(knobs_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("avail_cpu_time", self%avail_cpu_time, unit_internal)
    call self%write_key_val("beta", self%beta, unit_internal)
    call self%write_key_val("delt", self%delt, unit_internal)
    call self%write_key_val("delt_option", self%delt_option, unit_internal)
    call self%write_key_val("do_eigsolve", self%do_eigsolve, unit_internal)
    call self%write_key_val("eqzip_option", self%eqzip_option, unit_internal)
    call self%write_key_val("fapar", self%fapar, unit_internal)
    call self%write_key_val("faperp", self%faperp, unit_internal)
    call self%write_key_val("fbpar", self%fbpar, unit_internal)
    call self%write_key_val("fphi", self%fphi, unit_internal)
    call self%write_key_val("immediate_reset", self%immediate_reset, unit_internal)
    call self%write_key_val("k0", self%k0, unit_internal)
    call self%write_key_val("margin", self%margin, unit_internal)
    call self%write_key_val("margin_cpu_time", self%margin_cpu_time, unit_internal)
    call self%write_key_val("max_sim_time", self%max_sim_time, unit_internal)
    call self%write_key_val("neo_test", self%neo_test, unit_internal)
    call self%write_key_val("nstep", self%nstep, unit_internal)
    call self%write_key_val("rhostar", self%rhostar, unit_internal)
    call self%write_key_val("save_init_times", self%save_init_times, unit_internal)
    call self%write_key_val("save_timer_statistics", self%save_timer_statistics, unit_internal)
    call self%write_key_val("seed", self%seed, unit_internal)
    call self%write_key_val("teti", self%teti, unit_internal)
    call self%write_key_val("tite", self%tite, unit_internal)
    call self%write_key_val("trinity_linear_fluxes", self%trinity_linear_fluxes, unit_internal)
    call self%write_key_val("trinity_ql_fluxes", self%trinity_ql_fluxes, unit_internal)
    call self%write_key_val("use_old_diagnostics", self%use_old_diagnostics, unit_internal)
    call self%write_key_val("user_comments", self%user_comments, unit_internal)
    call self%write_key_val("wstar_units", self%wstar_units, unit_internal)
    call self%write_key_val("zeff", self%zeff, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_knobs_config

  !> Resets the config object to the initial empty state
  subroutine reset_knobs_config(self)
    class(knobs_config_type), intent(in out) :: self
    type(knobs_config_type) :: empty
    select type (self)
    type is (knobs_config_type)
       self = empty
    end select
  end subroutine reset_knobs_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_knobs_config(self)
    use mp, only: broadcast
    implicit none
    class(knobs_config_type), intent(in out) :: self
    call broadcast(self%avail_cpu_time)
    call broadcast(self%beta)
    call broadcast(self%delt)
    call broadcast(self%delt_option)
    call broadcast(self%do_eigsolve)
    call broadcast(self%eqzip_option)
    call broadcast(self%fapar)
    call broadcast(self%faperp)
    call broadcast(self%fbpar)
    call broadcast(self%fphi)
    call broadcast(self%immediate_reset)
    call broadcast(self%k0)
    call broadcast(self%margin)
    call broadcast(self%margin_cpu_time)
    call broadcast(self%max_sim_time)
    call broadcast(self%neo_test)
    call broadcast(self%nstep)
    call broadcast(self%rhostar)
    call broadcast(self%save_init_times)
    call broadcast(self%save_timer_statistics)
    call broadcast(self%seed)
    call broadcast(self%teti)
    call broadcast(self%tite)
    call broadcast(self%trinity_linear_fluxes)
    call broadcast(self%trinity_ql_fluxes)
    call broadcast(self%use_old_diagnostics)
    call broadcast(self%user_comments)
    call broadcast(self%wstar_units)
    call broadcast(self%zeff)

    call broadcast(self%exist)
  end subroutine broadcast_knobs_config

  !> Gets the default name for this namelist
  function get_default_name_knobs_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_knobs_config
    get_default_name_knobs_config = "knobs"
  end function get_default_name_knobs_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_knobs_config()
    implicit none
    logical :: get_default_requires_index_knobs_config
    get_default_requires_index_knobs_config = .false.
  end function get_default_requires_index_knobs_config

  !> Get the module level config instance
  pure function get_knobs_config()
    type(knobs_config_type) :: get_knobs_config
    get_knobs_config = knobs_config
  end function get_knobs_config
  
  !---------------------------------------
  ! Following is for the parameters config_type
  !---------------------------------------
  
  !> Reads in the parameters namelist and populates the member variables
  subroutine read_parameters_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(parameters_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 100000) :: user_comments
    real :: beta, k0, rhostar, teti, tite, zeff

    namelist /parameters/ beta, k0, rhostar, teti, tite, user_comments, zeff

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    beta = self%beta
    k0 = self%k0
    rhostar = self%rhostar
    teti = self%teti
    tite = self%tite
    user_comments = self%user_comments
    zeff = self%zeff

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = parameters)

    ! Now copy from local variables into type members
    self%beta = beta
    self%k0 = k0
    self%rhostar = rhostar
    self%teti = teti
    self%tite = tite
    self%user_comments = user_comments
    self%zeff = zeff

    self%exist = exist
  end subroutine read_parameters_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_parameters_config(self, unit)
    implicit none
    class(parameters_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("beta", self%beta, unit_internal)
    call self%write_key_val("k0", self%k0, unit_internal)
    call self%write_key_val("rhostar", self%rhostar, unit_internal)
    call self%write_key_val("teti", self%teti, unit_internal)
    call self%write_key_val("tite", self%tite, unit_internal)
    call self%write_key_val("user_comments", self%user_comments, unit_internal)
    call self%write_key_val("zeff", self%zeff, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_parameters_config

  !> Resets the config object to the initial empty state
  subroutine reset_parameters_config(self)
    class(parameters_config_type), intent(in out) :: self
    type(parameters_config_type) :: empty
    select type (self)
    type is (parameters_config_type)
       self = empty
    end select
  end subroutine reset_parameters_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_parameters_config(self)
    use mp, only: broadcast
    implicit none
    class(parameters_config_type), intent(in out) :: self
    call broadcast(self%beta)
    call broadcast(self%k0)
    call broadcast(self%rhostar)
    call broadcast(self%teti)
    call broadcast(self%tite)
    call broadcast(self%user_comments)
    call broadcast(self%zeff)

    call broadcast(self%exist)
  end subroutine broadcast_parameters_config
  
  !> Gets the default name for this namelist
  function get_default_name_parameters_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_parameters_config
    get_default_name_parameters_config = "parameters"
  end function get_default_name_parameters_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_parameters_config()
    implicit none
    logical :: get_default_requires_index_parameters_config
    get_default_requires_index_parameters_config = .false.
  end function get_default_requires_index_parameters_config

  !> Get the module level config instance
  function get_parameters_config()
    type(parameters_config_type) :: get_parameters_config
    get_parameters_config = parameters_config
  end function get_parameters_config

end module run_parameters
