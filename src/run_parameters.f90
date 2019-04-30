!> This module is basically a store for the input parameters that are specified in the namelists \a knobs and \a parameters.
!! In general, the names of the public variables in this module are the same as the name of the input parameter they correspond to.
module run_parameters
  use abstract_config, only: abstract_config_type
  
  implicit none

  private

  public :: init_run_parameters, finish_run_parameters
  public :: check_run_parameters, wnml_run_parameters
  public :: write_trinity_parameters
  public :: beta, zeff, tite, reset, immediate_reset
  public :: fphi, fapar, fbpar
  public :: code_delt_max, wunits, woutunits, tunits
  public :: nstep, wstar_units, eqzip, margin
  public :: secondary, tertiary, harris
  public :: ieqzip, k0
  public :: delt_option_switch, delt_option_hand
  public :: delt_option_auto, vnm_init
  public :: avail_cpu_time, margin_cpu_time
  public :: rhostar, neo_test
  public :: do_eigsolve
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


  !> Used to represent the input configuration of run_parameters
  !! This actually involves two namelists, "parameters" and "knobs".
  !!
  !! @todo These two namelists should probably be merged. This is actually
  !! implemented here by adding the contents of "parameters" to "knobs".
  !! This should be documented clearly. Eventually we should remove parameters
  type, extends(abstract_config_type) :: run_parameters_config_type
     ! Parameters
     !> \(\beta\) is the ratio of the reference pressure to the reference magnetic energy density, \(\beta=2\mu_0n_{ref}T_{ref}/B^2_{ref}\).
     real :: beta = 0.0
     !> Effective ionic charge appearing in electron collision frequency
     real :: zeff = 1.0
     !> Ratio of ion to electron temperature, only used when there is an adiabatic species
     real :: tite = 1.0
     !> Deprecated
     real :: teti = -100.0
     !> Unknown, only uses in [[init_g]] -- should be moved to init_g namelist.
     real :: k0 = 1.0
     !> Normalised gyro-radius, only used for low flow builds
     real :: rhostar = 3.0e-3
     !> Custom description of run to be added to netcdf output (new diagnostics only)
     character(len = 100000) :: user_comments = ''
     ! Knobs
     !> Multiplies \(\phi\) throughout.
     real :: fphi = 1.0
     !> Multiplies \(A_\|\) throughout.
     real :: fapar = 0.0
     !> Multiplies \(B_\|\) throughout.
     real :: fbpar = -1.0
     !> Initial timestep
     real :: delt = 0.1
     !> Maximum number of steps to take
     integer :: nstep = 100
     !>If true makes timestep proportional to ky*rho. Only sensible for linear runs.
     logical :: wstar_units = .false.
     !> Advanced option to freeze evolution of certain modes
     logical :: eqzip = .false.
     !> Determines how initial timestep is set
     character(len = 20) :: delt_option = 'default'
     !> Not used anywhere -- to be removed.
     real :: margin = 0.05
     !> Advanced
     logical :: secondary = .true.
     !> Advanced
     logical :: tertiary = .false.
     !> Deprecated
     real :: faperp = 0.0
     !> Advanced
     logical :: harris = .false.
     !> The maximum wall clock time available to GS2
     real :: avail_cpu_time = 1.0e10
     !> How close to avail_cpu_time can we go before trying to stop the run cleanly
     real :: margin_cpu_time = 300.0
     !> What types of modes are frozen with eqzip
     character(len = 20) :: eqzip_option = 'none'
     !> Lowflow only.
     logical :: neo_test = .false.
     !> If true and running linearly, return linear diffusive flux estimates to Trinity.
     logical :: trinity_linear_fluxes = .false.
     !> If true then use eigensolver instead of initial value solver.
     logical :: do_eigsolve = .false.
     !> Determines the behaviour when the CFL condition is broken in the nonlinear term:
     logical :: immediate_reset = .true.
     !> If true use original diagnostics rather than new form
     logical :: use_old_diagnostics = .false.
     !> Unknown
     logical :: trinity_ql_fluxes = .false.
   contains
     procedure, public :: init => init_run_parameters_config
     procedure, public :: read => read_run_parameters_config
     procedure, public :: write => write_run_paramters_config
     procedure, public :: reset => reset_run_parameters_config
     procedure, public :: broadcast => broadcast_run_parameters_config
  end type run_parameters_config_type
  
  type(run_parameters_config_type) :: config_internal

  real :: beta, zeff, tite
  real :: fphi, fapar, fbpar, faperp
  real :: delt, delt_prev1, delt_prev2
  real :: code_delt_max, user_delt_max, margin
  real, dimension (:), allocatable :: wunits, woutunits, tunits
  real, dimension (2) :: vnm_init
  real :: avail_cpu_time, margin_cpu_time
  integer :: nstep
  logical :: reset=.false.
  logical :: immediate_reset
  logical :: wstar_units, eqzip
  logical :: secondary, tertiary, harris
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

  integer, allocatable :: ieqzip(:,:)
  integer :: eqzip_option_switch
  integer, parameter :: &
       eqzip_option_none = 1, &
       eqzip_option_secondary = 2, &
       eqzip_option_tertiary = 3, &
       eqzip_option_equilibrium = 4

contains

  !> FIXME : Add documentation  
  subroutine check_run_parameters(report_unit)
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

    if (eqzip) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('eqzip = T in the knobs namelist.')")
       write (report_unit, fmt="('This freezes some modes in time for a secondary stability analysis.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
       if (secondary) write (report_unit, fmt="('Mode with kx = 0, ky = ky_min fixed in time')")
       if (tertiary)  write (report_unit, fmt="('Mode with ky = 0, kx = kx_min fixed in time')")
    end if

    write (report_unit, *) 
    if(immediate_reset)then
       write (report_unit, fmt="('The time step will be reset immediately after cfl violation detected.')") 
    else
       write (report_unit, fmt="('The time step will be reset just before the next time step after cfl violation detected.')") 
    endif
    write (report_unit, *) 
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
       write (unit, fmt="(' nstep = ',i8)") nstep
       write (unit, fmt="(' wstar_units = ',L1)") wstar_units
       if (eqzip) then
          write (unit, fmt="(' eqzip = ',L1)") eqzip
          write (unit, fmt="(' secondary = ',L1)") secondary
          write (unit, fmt="(' tertiary = ',L1)") tertiary
       end if
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
  subroutine init_run_parameters
    use kt_grids, only: init_kt_grids, naky, nakx => ntheta0
    use gs2_time, only: init_delt, user2code
    implicit none

    if (initialized) return
    initialized = .true.

    call read_parameters

    call init_kt_grids
    call init_delt (delt, delt_prev1, delt_prev2)
    call user2code (user_delt_max, code_delt_max)

    if(.not. allocated(wunits)) allocate (wunits(naky))
    if(.not. allocated(woutunits)) allocate (woutunits(naky))
    if(.not. allocated(tunits)) allocate (tunits(naky))

! omega_* normalization of time: 
    call adjust_time_norm

    if(.not.allocated(ieqzip)) allocate(ieqzip(nakx,naky))
    ieqzip(1:nakx,1:naky)=1
    select case (eqzip_option_switch)
    case (eqzip_option_secondary)
       ! suppress evolution of secondary mode
       ieqzip(1,2) = 0
    case (eqzip_option_tertiary)
       ! suppress evolution of tertiary mode
       ieqzip(2,1) = 0
       ieqzip(nakx,1) = 0
    case (eqzip_option_equilibrium)
       ! suppress evolution of 1D equilibrium (x dependent)
       ieqzip(1:nakx,1) = 0
    end select
  end subroutine init_run_parameters

  !> FIXME : Add documentation    
  subroutine read_parameters(config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use mp, only: proc0, broadcast
    use gs2_time, only: dt_not_set
    use gs2_save, only: init_dt, init_vnm
    use text_options, only: text_option, get_option_value
    use kt_grids, only: gryfx
    implicit none
    type(run_parameters_config_type), intent(in), optional :: config_in
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
    integer :: ierr, istatus, in_file
    real :: delt0, delt1, delt2
    real, dimension (2) :: vnm_saved

    real :: teti  ! for back-compatibility

    if (present(config_in)) config_internal = config_in

    call config_internal%init()

    ! Copy out internal values into module level parameters
    beta  = config_internal%beta 
    zeff  = config_internal%zeff 
    tite  = config_internal%tite 
    teti  = config_internal%teti 
    k0    = config_internal%k0   
    rhostar  = config_internal%rhostar 
    user_comments  = config_internal%user_comments 
    fphi  = config_internal%fphi 
    fapar  = config_internal%fapar 
    fbpar  = config_internal%fbpar 
    delt  = config_internal%delt 
    nstep  = config_internal%nstep 
    wstar_units  = config_internal%wstar_units 
    eqzip  = config_internal%eqzip 
    delt_option  = config_internal%delt_option 
    margin = config_internal%margin
    secondary  = config_internal%secondary 
    tertiary  = config_internal%tertiary 
    faperp  = config_internal%faperp 
    harris  = config_internal%harris 
    avail_cpu_time  = config_internal%avail_cpu_time 
    margin_cpu_time  = config_internal%margin_cpu_time 
    eqzip_option  = config_internal%eqzip_option 
    neo_test  = config_internal%neo_test 
    trinity_linear_fluxes  = config_internal%trinity_linear_fluxes 
    do_eigsolve  = config_internal%do_eigsolve 
    immediate_reset  = config_internal%immediate_reset 
    use_old_diagnostics  = config_internal%use_old_diagnostics 
    trinity_ql_fluxes  = config_internal%trinity_ql_fluxes 
  
    if(proc0) in_file = input_unit_exist("parameters", rpexist)
    knexist = config_internal%exist

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

    !> If we have non-zero beta then alert the user that they have some
    !! of the perturbed magnetic fields disabled. This may be intended
    !! behaviour, so we throw a warning. However, if both are disabled,
    !! we set beta=0. This will make the simulation electrostatic, as 
    !! probably intended if both apar and bpar are set to zero.
    if(fapar.eq.0 .and. fbpar.eq.0 .and. beta.ne.0) then
       if(proc0) then
          ierr = error_unit()
          write(ierr,'("Warning: Both fapar and fbpar are zero: setting beta = 0.")')
       endif
       beta = 0.0
    endif
    if((beta.ne.0).and.proc0) then
       ierr = error_unit()
       if(fapar.eq.0) write(ierr,'("Warning: Running with finite beta but fapar=0.")')
       if(fbpar.eq.0) write(ierr,'("Warning: Running with finite beta but fbpar=0.")')
    endif
    
    if (eqzip) then
       if (secondary .and. tertiary) then
          ierr = error_unit()
          write (ierr, *) 'Forcing secondary = FALSE'
          write (ierr, *) 'because you have chosen tertiary = TRUE'
          secondary = .false.
       end if
       if (secondary .and. harris) then
          ierr = error_unit()
          write (ierr, *) 'Forcing secondary = FALSE'
          write (ierr, *) 'because you have chosen harris = TRUE'
          secondary = .false.
       end if
       if (tertiary .and. harris) then
          ierr = error_unit()
          write (ierr, *) 'Forcing tertiary = FALSE'
          write (ierr, *) 'because you have chosen harris = TRUE'
          tertiary = .false.
       end if
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
       delt = sqrt(2.0)*delt/2.0
       nstep = nstep*2
    end if

    user_delt_max = delt

    !! @todo: Add input variables to allow these to be set directly
    delt_prev1 = dt_not_set
    delt_prev2 = dt_not_set
    
    if (delt_option_switch == delt_option_auto) then
       vnm_init = 1.0
       call init_vnm (vnm_saved, istatus)
       if (istatus == 0) vnm_init = vnm_saved
       call init_dt (delt0, delt1, delt2, istatus)
       if (istatus == 0) then
          delt  = delt0
          delt_prev1 = delt1
          delt_prev2 = delt2
       endif
    endif
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
  subroutine adjust_time_norm
    use file_utils, only: error_unit
    use mp, only: proc0
    use kt_grids, only: aky
    implicit none
!CMR: Sep 2010
! Attempt to understand time normalisation variables, which are arrays(naky)
!    TUNITS: DT(KY)=TUNITS(KY).CODE_DT
!            This is a generally very useful variable to store ky dependent 
!            timestep in the code time normalisation.
!            Used to multiply ky independent source terms on RHS of GKE.
!    WUNITS: WUNITS(KY)=AKY(KY)*TUNITS(KY)/2
!            Auxiliary variable.  Used to save compute operations when 
!            evaluating source terms on RHS of GKE that are proportional to ky.
!            !! The Mysterious factor 1/2 Explained !!
!            The factor of 1/2 arises because those source terms were first
!            specified using the normalisation Tref=mref vtref^2 
! [R Numata et al, "AstroGK: Astrophysical gyrokinetics code", JCP, 2010].
!CMRend
    if (wstar_units) then
       wunits = 1.0
       where (aky /= 0.0)
          tunits = 2.0/aky
       elsewhere
          tunits = 0.0
       end where
       if (any(tunits == 0.0) .and. proc0) then
          write (error_unit(), *) &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
          print *, &
               "WARNING: wstar_units=.true. and aky=0.0: garbage results"
       end if
!CMR: Sep 2010
!  Changes to allow wstar_units to be used
!CMRend
    else
       tunits = 1.0
       wunits = aky/2.0
    end if
    woutunits = 1.0/tunits
  end subroutine adjust_time_norm

  !> FIXME : Add documentation    
  subroutine finish_run_parameters
    implicit none

    if (allocated(wunits)) deallocate (wunits, woutunits, tunits)
    if (allocated(ieqzip)) deallocate (ieqzip)

    initialized = .false.
  end subroutine finish_run_parameters

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------

  !> Fully initialises the config object.
  subroutine init_run_parameters_config(self, index)
    implicit none
    class(run_parameters_config_type), intent(in out) :: self
    integer, intent(in), optional :: index

    if (self%is_initialised()) return
    call self%setup(name = 'knobs', requires_index = .false., index = index)
    call self%read()
    call self%broadcast()
  end subroutine init_run_parameters_config
  
  !> Reads in the parameters and knobs namelists and populates the member variables
  subroutine read_run_parameters_config(self)
    use file_utils, only: input_unit_exist
    use mp, only: proc0
    implicit none
    class(run_parameters_config_type), intent(in out) :: self
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
  end subroutine read_run_parameters_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_run_paramters_config(self, unit)
    implicit none
    class(run_parameters_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

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
  end subroutine write_run_paramters_config

  !> Resets the config object to the initial empty state
  subroutine reset_run_parameters_config(self)
    class(run_parameters_config_type), intent(in out) :: self
    type(run_parameters_config_type) :: empty
    select type (self)
    type is (run_parameters_config_type)
       self = empty
    end select
  end subroutine reset_run_parameters_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_run_parameters_config(self)
    use mp, only: broadcast
    implicit none
    class(run_parameters_config_type), intent(in out) :: self
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
  end subroutine broadcast_run_parameters_config  
end module run_parameters
