!> Types with bound routines for representing/holding the "input configuration"
!! set from the contents of the provided input file. This allows us to override
!! what's in the input file by changing values in instances of these types.
!!
!! @todo Allow a specific input file to be specified rather than relying on file_utils choice
!!
!! @todo Testing
module config_types
  use abstract_config, only: abstract_config_type
  implicit none
  private

  public :: run_parameters_config_type

  !> Used to represent the input configuration of run_parameters
  !! This actually involves two namelists, "parameters" and "knobs".
  !!
  !! @todo These two namelists should probably be merged. This is actually
  !! implemented here by adding the contents of "parameters" to "knobs".
  !! This should be documented clearly. Eventually we should remove parameters
  type, extends(abstract_config_type) :: run_parameters_config_type
     ! Parameters
     !> \(\beta\) is the ratio of the reference pressure to the reference magnetic energy density, \(\beta=2\mu_0n_{ref}T_{ref}/B^2_{ref}\).
     real :: beta
     !> Effective ionic charge appearing in electron collision frequency
     real :: zeff
     !> Ratio of ion to electron temperature, only used when there is an adiabatic species
     real :: tite
     !> Deprecated
     real :: teti
     !> Unknown, only uses in [[init_g]] -- should be moved to init_g namelist.
     real :: k0
     !> Normalised gyro-radius, only used for low flow builds
     real :: rhostar
     !> Custom description of run to be added to netcdf output (new diagnostics only)
     character(len = 100000) :: user_comments
     ! Knobs
     !> Multiplies \(\phi\) throughout.
     real :: fphi
     !> Multiplies \(A_\|\) throughout.
     real :: fapar
     !> Multiplies \(B_\|\) throughout.
     real :: fbpar
     !> Initial timestep
     real :: delt
     !> Maximum number of steps to take
     integer :: nstep
     !>If true makes timestep proportional to ky*rho. Only sensible for linear runs.
     logical :: wstar_units
     !> Advanced option to freeze evolution of certain modes
     logical :: eqzip
     !> Determines how initial timestep is set
     character(len = 20) :: delt_option
     !> Not used anywhere -- to be removed.
     real :: margin
     !> Advanced
     logical :: secondary
     !> Advanced
     logical :: tertiary
     !> Deprecated
     real :: faperp
     !> Advanced
     logical :: harris
     !> The maximum wall clock time available to GS2
     real :: avail_cpu_time
     !> How close to avail_cpu_time can we go before trying to stop the run cleanly
     real :: margin_cpu_time
     !> What types of modes are frozen with eqzip
     character(len = 20) :: eqzip_option
     !> Lowflow only.
     logical :: neo_test
     !> If true and running linearly, return linear diffusive flux estimates to Trinity.
     logical :: trinity_linear_fluxes
     !> If true then use eigensolver instead of initial value solver.
     logical :: do_eigsolve
     !> Determines the behaviour when the CFL condition is broken in the nonlinear term:
     logical :: immediate_reset
     !> If true use original diagnostics rather than new form
     logical :: use_old_diagnostics
     !> Unknown
     logical :: trinity_ql_fluxes
   contains
     procedure, public :: init => init_run_parameters
     procedure, public :: read => read_run_parameters
     procedure, public :: write => write_run_paramters
     procedure, public :: broadcast => broadcast_run_parameters
     procedure, public :: set_to_defaults => set_to_defaults_run_parameters
  end type run_parameters_config_type
contains

  !> Almost entirely generic currently
  subroutine init_run_parameters(self, index)
    implicit none
    class (run_parameters_config_type), intent(in out) :: self
    integer, intent(in), optional :: index

    if (self%is_initialised()) return
    call self%setup(name = 'knobs', requires_index = .false., index = index)
    call self%set_to_defaults()
    call self%read()
    call self%broadcast()
  end subroutine init_run_parameters

  subroutine read_run_parameters(self)
    use file_utils, only: input_unit_exist
    use mp, only: proc0
    implicit none
    class (run_parameters_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Parameters
    real :: beta, zeff, tite, teti, k0, rhostar
    character(len = 100000) :: user_comments
    ! Knobs
    real :: fphi, fapar, fbpar, delt, faperp, avail_cpu_time, margin_cpu_time, margin
    integer :: nstep
    logical :: wstar_units, eqzip, secondary, tertiary, harris, neo_test
    logical :: trinity_linear_fluxes, do_eigsolve, immediate_reset
    logical :: use_old_diagnostics, trinity_ql_fluxes
    character(len = 20) :: delt_option, eqzip_option

    namelist /parameters/ beta, zeff, tite, teti, k0, rhostar, user_comments
    namelist /knobs/ fphi, fapar, fbpar, delt, nstep, wstar_units, eqzip, &
         delt_option, margin, secondary, tertiary, faperp, harris, &
         avail_cpu_time, margin_cpu_time, eqzip_option, neo_test, &
         trinity_linear_fluxes, do_eigsolve, immediate_reset, &
         use_old_diagnostics,  trinity_ql_fluxes, beta, zeff, tite, &
         teti, k0, rhostar, user_comments

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    beta  = self%beta
    zeff  = self%zeff
    tite  = self%tite
    teti  = self%teti
    k0    = self%k0
    rhostar  = self%rhostar
    user_comments  = self%user_comments
    fphi  = self%fphi
    fapar  = self%fapar
    fbpar  = self%fbpar
    delt  = self%delt
    nstep  = self%nstep
    wstar_units  = self%wstar_units
    eqzip  = self%eqzip
    delt_option  = self%delt_option
    margin = self%margin
    secondary  = self%secondary
    tertiary  = self%tertiary
    faperp  = self%faperp
    harris  = self%harris
    avail_cpu_time  = self%avail_cpu_time
    margin_cpu_time  = self%margin_cpu_time
    eqzip_option  = self%eqzip_option
    neo_test  = self%neo_test
    trinity_linear_fluxes  = self%trinity_linear_fluxes
    do_eigsolve  = self%do_eigsolve
    immediate_reset  = self%immediate_reset
    use_old_diagnostics  = self%use_old_diagnostics
    trinity_ql_fluxes  = self%trinity_ql_fluxes

    ! Now deal with old parameters namelist
    in_file = input_unit_exist('parameters', exist)
    if (exist) then
       read(in_file, nml = parameters)
       !write(error_unit(), '("Warning: The parameters namelist is deprecated, the contents has moved to the knobs namelist.")')
    end if

    ! Now read in the main namelist
    in_file = input_unit_exist(self%name, exist)
    if (exist) read(in_file, nml = knobs)

    ! Now copy from local variables into type members
    self%beta  = beta
    self%zeff  = zeff
    self%tite  = tite
    self%teti  = teti
    self%k0    = k0
    self%rhostar  = rhostar
    self%user_comments  = user_comments
    self%fphi  = fphi
    self%fapar  = fapar
    self%fbpar  = fbpar
    self%delt  = delt
    self%nstep  = nstep
    self%wstar_units  = wstar_units
    self%eqzip  = eqzip
    self%delt_option  = delt_option
    self%margin = margin
    self%secondary  = secondary
    self%tertiary  = tertiary
    self%faperp  = faperp
    self%harris  = harris
    self%avail_cpu_time  = avail_cpu_time
    self%margin_cpu_time  = margin_cpu_time
    self%eqzip_option  = eqzip_option
    self%neo_test  = neo_test
    self%trinity_linear_fluxes  = trinity_linear_fluxes
    self%do_eigsolve  = do_eigsolve
    self%immediate_reset  = immediate_reset
    self%use_old_diagnostics  = use_old_diagnostics
    self%trinity_ql_fluxes  = trinity_ql_fluxes

    self%exist = exist
  end subroutine read_run_parameters

  subroutine write_run_paramters(self, unit)
    implicit none
    class(run_parameters_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    if (.not.self%is_initialised()) return

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("beta", self%beta, unit_internal)
    call self%write_key_val("zeff", self%zeff, unit_internal)
    call self%write_key_val("tite", self%tite, unit_internal)
    call self%write_key_val("teti", self%teti, unit_internal)
    call self%write_key_val("k0", self%k0, unit_internal)
    call self%write_key_val("rhostar", self%rhostar, unit_internal)
    call self%write_key_val("user_comments", self%user_comments, unit_internal)
    call self%write_key_val("fphi", self%fphi, unit_internal)
    call self%write_key_val("fapar", self%fapar, unit_internal)
    call self%write_key_val("fbpar", self%fbpar, unit_internal)
    call self%write_key_val("delt", self%delt, unit_internal)
    call self%write_key_val("nstep", self%nstep, unit_internal)
    call self%write_key_val("wstar_units", self%wstar_units, unit_internal)
    call self%write_key_val("eqzip", self%eqzip, unit_internal)
    call self%write_key_val("delt_option", self%delt_option, unit_internal)
    call self%write_key_val("margin", self%margin, unit_internal)
    call self%write_key_val("secondary", self%secondary, unit_internal)
    call self%write_key_val("tertiary", self%tertiary, unit_internal)
    call self%write_key_val("faperp", self%faperp, unit_internal)
    call self%write_key_val("harris", self%harris, unit_internal)
    call self%write_key_val("avail_cpu_time", self%avail_cpu_time, unit_internal)
    call self%write_key_val("margin_cpu_time", self%margin_cpu_time, unit_internal)
    call self%write_key_val("eqzip_option", self%eqzip_option, unit_internal)
    call self%write_key_val("neo_test", self%neo_test, unit_internal)
    call self%write_key_val("trinity_linear_fluxes", self%trinity_linear_fluxes, unit_internal)
    call self%write_key_val("do_eigsolve", self%do_eigsolve, unit_internal)
    call self%write_key_val("immediate_reset", self%immediate_reset, unit_internal)
    call self%write_key_val("use_old_diagnostics", self%use_old_diagnostics, unit_internal)
    call self%write_key_val("trinity_ql_fluxes", self%trinity_ql_fluxes, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_run_paramters

  subroutine set_to_defaults_run_parameters(self)
    implicit none
    class (run_parameters_config_type), intent(in out) :: self
    self%beta  = 0.0
    self%zeff  = 1.0
    self%tite  = 1.0
    self%teti  = -100.0
    self%k0    = 1.0
    self%rhostar  = 3.0e-3
    self%user_comments  = ''
    self%fphi  = 1.0
    self%fapar  = 0.0
    self%fbpar  = -1.0
    self%delt  = 0.1
    self%nstep  = 100
    self%wstar_units = .false.
    self%eqzip  = .false.
    self%delt_option  = 'default'
    self%margin = 0.05
    self%secondary  = .true.
    self%tertiary  = .false.
    self%faperp  = 0.0
    self%harris  = .false.
    self%avail_cpu_time  = 1.0e10
    self%margin_cpu_time  = 300.0
    self%eqzip_option  = 'none'
    self%neo_test  = .false.
    self%trinity_linear_fluxes  = .false.
    self%do_eigsolve  = .false.
    self%immediate_reset  = .true.
    self%use_old_diagnostics  = .false.
    self%trinity_ql_fluxes  = .false.
  end subroutine set_to_defaults_run_parameters

  subroutine broadcast_run_parameters(self)
    use mp, only: broadcast
    implicit none
    class (run_parameters_config_type), intent(in out) :: self
    call broadcast(self%beta)
    call broadcast(self%zeff)
    call broadcast(self%tite)
    call broadcast(self%teti)
    call broadcast(self%k0)
    call broadcast(self%rhostar)
    call broadcast(self%user_comments)
    call broadcast(self%fphi)
    call broadcast(self%fapar)
    call broadcast(self%fbpar)
    call broadcast(self%delt)
    call broadcast(self%nstep)
    call broadcast(self%wstar_units)
    call broadcast(self%eqzip)
    call broadcast(self%delt_option)
    call broadcast(self%margin)
    call broadcast(self%secondary)
    call broadcast(self%tertiary)
    call broadcast(self%faperp)
    call broadcast(self%harris)
    call broadcast(self%avail_cpu_time)
    call broadcast(self%margin_cpu_time)
    call broadcast(self%eqzip_option)
    call broadcast(self%neo_test)
    call broadcast(self%trinity_linear_fluxes)
    call broadcast(self%do_eigsolve)
    call broadcast(self%immediate_reset)
    call broadcast(self%use_old_diagnostics)
    call broadcast(self%trinity_ql_fluxes)

    call broadcast(self%exist)
  end subroutine broadcast_run_parameters

end module config_types
