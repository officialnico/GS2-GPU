!> Routines for implementing the model collision operator
!! defined by Barnes, Abel et al. 2009. The collision operator causes
!! physically motivated smoothing of structure in velocity space which
!! is necessary to prevent buildup of structure at fine scales in
!! velocity space, while conserving energy and momentum.
module collisions
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use redistribute, only: redist_type

  implicit none

  private 

  public :: init_collisions, finish_collisions
  public :: read_parameters, wnml_collisions, check_collisions
  public :: solfp1
  public :: reset_init
  public :: dtot, fdf, fdb, vnmult, vnfac
  public :: ncheck, vnslow, vary_vnew
  public :: etol, ewindow, etola, ewindowa
  public :: timesteps_between_collisions
  public :: init_lorentz, init_ediffuse
  public :: init_lorentz_conserve, init_diffuse_conserve
  public :: init_lorentz_error, collision_model_switch
  public :: colls, hyper_colls, heating, adjust
  public :: use_le_layout, collision_model_ediffuse
  public :: collision_model_lorentz, collision_model_none
  public :: collision_model_lorentz_test, collision_model_full
  public :: split_collisions

  public :: collisions_config_type
  public :: set_collisions_config
  public :: get_collisions_config
  
  interface solfp1
     module procedure solfp1_le_layout
     module procedure solfp1_standard_layout
  end interface

  interface solfp_lorentz
     module procedure solfp_lorentz_le_layout
     module procedure solfp_lorentz_standard_layout
  end interface

  interface conserve_lorentz
     module procedure conserve_lorentz_le_layout
     module procedure conserve_lorentz_standard_layout
  end interface

  interface conserve_diffuse
     module procedure conserve_diffuse_le_layout
     module procedure conserve_diffuse_standard_layout
  end interface 

  ! knobs
  logical :: use_le_layout = .true.
  logical :: const_v, conserve_moments
  logical :: conservative, resistivity
  integer :: collision_model_switch
  integer :: lorentz_switch, ediff_switch
  logical :: adjust
  logical :: heating
  logical :: hyper_colls
  logical :: ei_coll_only
  logical :: test
  logical :: special_wfb_lorentz
  logical :: vpar_zero_mean
  logical :: conserve_forbid_zero
  logical :: le_nxilim_fix
  
  integer, parameter :: collision_model_lorentz = 1      ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_none = 3
  integer, parameter :: collision_model_lorentz_test = 5 ! if this changes, check gs2_diagnostics
  integer, parameter :: collision_model_full = 6
  integer, parameter :: collision_model_ediffuse = 7

  integer, parameter :: lorentz_scheme_default = 1
  integer, parameter :: lorentz_scheme_old = 2

  integer, parameter :: ediff_scheme_default = 1
  integer, parameter :: ediff_scheme_old = 2

  real, dimension (2) :: vnmult = -1.0
  integer :: ncheck
  logical :: vary_vnew
  real :: vnfac, vnslow
  real :: etol, ewindow, etola, ewindowa
  integer :: timesteps_between_collisions
  logical :: force_collisions, has_lorentz, has_diffuse
  real, dimension (2) :: vnm_init = 1.0

  real, dimension (:,:,:), allocatable :: dtot
  ! (-ntgrid:ntgrid,nlambda,max(ng2,nlambda-ng2)) lagrange coefficients for derivative error estimate

  real, dimension (:,:), allocatable :: fdf, fdb
  ! (-ntgrid:ntgrid,nlambda) finite difference coefficients for derivative error estimate

  real, dimension (:,:,:), allocatable :: vnew, vnew_s, vnew_D, vnew_E, delvnew
  ! (naky,negrid,nspec) replicated

  real, dimension (:,:,:), allocatable :: vpdiff
  ! (-ntgrid:ntgrid,2,nlambda) replicated

  ! only for hyper-diffusive collisions
  real, dimension (:,:,:,:), allocatable :: vnewh
  ! (-ntgrid:ntgrid,ntheta0,naky,nspec) replicated

  ! only for momentum conservation due to Lorentz operator (8.06)
  real, dimension(:,:,:), allocatable :: s0, w0, z0
  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  real, dimension (:,:,:), allocatable :: s0le, w0le, z0le
  ! Copies of aj0, aj1 and vpa from dist_fn_arrays stored in
  ! le layout shaped arrays.
  real, dimension (:,:,:), allocatable :: aj0le, vperp_aj1le, vpa_aj0_le
  ! finish

  ! needed for momentum and energy conservation due to energy diffusion (3.08)
  real, dimension(:,:,:), allocatable :: bs0, bw0, bz0

  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  real, dimension (:,:,:), allocatable :: bs0le, bw0le, bz0le
  ! finish

  real :: cfac

  ! The following (between the start and finish comments) are only used for LE layout
  ! start
  ! only for lorentz
  real, dimension (:,:,:), allocatable :: c1le, betaale, qle, d1le, h1le
  ! only for energy diffusion
  real, dimension (:,:,:), allocatable :: ec1le, ebetaale, eqle
  ! finish
  ! The following (between the start and finish comments) are only used for none LE layout
  ! start
  ! only for lorentz
  real, dimension (:,:), allocatable :: c1, betaa, ql, d1, h1
  ! only for energy diffusion
  real, dimension (:,:), allocatable :: ec1, ebetaa, eql
  ! finish

  real, dimension (:, :), allocatable :: pitch_weights

  logical :: drag = .false.
  logical :: colls = .true.
  logical :: split_collisions

  logical :: hypermult
  logical :: initialized = .false.
  logical :: ediffinit = .false.
  logical :: lzinit = .false., leinit = .false.
  logical :: accelerated_x = .false.
  logical :: accelerated_v = .false.
  logical :: exist

  !> Used to represent the input configuration of collisions
  type, extends(abstract_config_type) :: collisions_config_type
     ! namelist : collisions_knobs
     ! indexed : false
     !> If true (default) then transform from the gyro-averaged
     !> distribution function \(g\) evolved by GS2 to the
     !> non-Boltzmann part of \(\delta f\), (\(h\)), when applying the
     !> collision operator. This is the physically appropriate choice,
     !> this parameter is primarily for numerical testing.
     logical :: adjust = .true.
     !> Factor multipyling the finite Larmor radius terms in the
     !> collision operator. This term is essentially just classical
     !> diffusion. Set `cfac` to 0 to turn off this diffusion.
     !>
     !> @note Default changed to 1.0 in order to include classical
     !> diffusion April 18 2006
     real :: cfac = 1.0
     !> Selects the collision model used in the simulation. Can be one
     !> of
     !>
     !> - `'default'` : Include both pitch angle scattering and energy
     !> diffusion.
     !> - `'collisionless'` : Disable the collision operator.
     !> - `'none'` : Equivalent to `'collisionless'`.
     !> - `'lorentz'` : Only pitch angle scattering.
     !> - `'lorentz-test'` : Only pitch angle scattering. For testing,
     !> disables some components of the operator.
     !> - `'ediffuse'` : Only energy diffusion.
     !> - `'krook'` : "Home made" Krook operator (this is not recommended).
     !>
     !> If no species have a non-zero collision frequency, `vnewk`, then
     !> the collision operator is also automatically disabled.
     character(len = 20) :: collision_model = 'default'
     !> If true (default) then guarantee exact conservation
     !> properties.
     logical :: conservative = .true.
     !> If true (default) then forces conservation corrections to zero
     !> in the forbidden region to avoid introducing unphysical
     !> non-zero values for the distribution function in the forbidden
     !> region of trapped particles.
     !>
     !> @todo Confirm above documentation is accurate.
     !>
     !> @note The conserving terms calculated as part of the field
     !> particle collision operator should respect the forbidden of
     !> the distribution function for trapped particles. This is a
     !> cosmetic change, but has the result that plots of the
     !> distribution function for trapped particles makes sense. Note
     !> that terms involving vpa do not need to be modified Because
     !> vpa = 0 in the forbidden region because of this explicit
     !> forbid statements have not been added to the drag term
     !> involving apar.
     logical :: conserve_forbid_zero = .true.
     !> If true (default) then guarantee collision operator conserves
     !> momentum and energy.
     !>
     !> @todo Clarify the difference with [[collisions_knobs:conservative]].
     !>
     !> @note Default changed to reflect improved momentum and energy
     !> conversation 07/08
     logical :: conserve_moments = .true.  
     !> If true (not the default) then use the thermal velocity to
     !> evaluate the collision frequencies to be used. This results in
     !> an energy independent collision frequency being used for all
     !> species.
     logical :: const_v = .false.
     !> Controls how the coefficients in the matrix representing the
     !> energy diffusion operator are obtained. Can be one of
     !>
     !> - `'default'` : Use a conservative scheme.
     !> - `'old'` : Use the original non-conservative scheme.
     !>
     !> @todo Consider deprecating/removing the `'old'` option.
     character(len = 20) :: ediff_scheme = 'default'
     !> If true (not the default) then force the collision frequency
     !> used for all non-electron species to zero and force all
     !> electron-electron terms to zero.
     logical :: ei_coll_only = .false.
     !> Only used in [[get_verr]] as a part of the adaptive
     !> collisionality algorithm. Sets the maximum relative error
     !> allowed, above which the collision frequency must be
     !> increased.
     !>
     !> @todo Confirm this is really to set the relative error limit.
     real :: etol = 2.e-2
     !> Only used in [[get_verr]] as a part of the adaptive
     !> collisionality algorithm. Sets the maximum absolute error
     !> allowed, above which the collision frequency must be
     !> increased.
     !>
     !> @todo Confirm this is really to set the absolute error limit.     
     real :: etola = 2.e-2
     !> Only used in [[get_verr]] as a part of the adaptive
     !> collisionality algorithm. Sets an offset to apply to the
     !> relative error limit set by [[collisions_knobs:etol]]. This is used to provide
     !> hysteresis is the adaptive collisionality algorithm so to
     !> avoid adjusting the collision frequency up and down every step
     !> (similar to [[reinit_knobs:delt_cushion]]).
     !>
     !> @todo Confirm the above description.
     real :: ewindow = 1.e-2
     !> Only used in [[get_verr]] as a part of the adaptive
     !> collisionality algorithm. Sets an offset to apply to the
     !> absolute error limit set by [[collisions_knobs:etola]]. This is used to provide
     !> hysteresis is the adaptive collisionality algorithm so to
     !> avoid adjusting the collision frequency up and down every step
     !> (similar to the [[reinit_knobs:delt_cushion]]).
     !>
     !> @todo Confirm the above description.     
     real :: ewindowa = 1.e-2
     !> Currently we skip the application of the selected collision
     !> operator for species with zero collision frequency and disable
     !> collisions entirely if no species have non-zero collision
     !> frequency. This is generally sensible, but it can be useful to
     !> force the use of the collisional code path in some situations
     !> such as in code testing. Setting this flag to `.true.` forces
     !> the selected collision model operator to be applied even if
     !> the collision frequency is zero.
     logical :: force_collisions = .false.
     !> If true (not the default) then calculate collisional heating
     !> when applying the collion operator. This is purely a
     !> diagnostic calculation. It should not change the evolution.
     !>
     !> @todo : Verify this does not influence the evolution.
     logical :: heating = .false.
     !> If true (not the default) then multiply the hyper collision
     !> frequency by the species' collision frequency
     !> [[species_parameters:nu_h]]. This only impacts the pitch
     !> angle scattering operator.
     !>
     !> @note The hyper collision frequency is only non-zero if any
     !> species have a non-zero `nu_h` value set in the input. If any
     !> are set then the hyper collision frequency is simply `nu_h *
     !> kperp2 * kperp2` (where `kperp2` here is normalised to the
     !> maximum `kperp2`).
     logical :: hypermult = .false.
     !> This flag controls the limit of ixi loops in the le_layout
     !> decomposition.  When true (default) the limit is set to `nxi +
     !> 1`, which has been observed to remove numerical instabilities.
     !>
     !> @note The "old" behavior was a limit of `nxi`, which is
     !> incorrect. With this limit and trapped particles, there are
     !> numerical instabilities caused by missing terms for totally
     !> trapped particles. By changing the limit to `nxi + 1` these
     !> numerical instabilities are removed and the le_layout gives
     !> the same result as the standard layout. The limit `nxi + 1`
     !> should used as we then include the `ixi` index for the
     !> duplicate `vpar=0` point which must be included For velocity
     !> integration to give the correct value.
     !>
     !> @todo Consider removing this flag and just use the correct
     !> limit.  The motivation for retaining this flag is to allow old
     !> results to be reproduced.
     logical :: le_nxilim_fix = .true.
     !> Controls how the coefficients in the matrix representing the
     !> pitch angle scattering operator are obtained. Can be one of
     !>
     !> - `'default'` : Use a conservative scheme.
     !> - `'old'` : Use the original non-conservative scheme.
     !>
     !> @todo Consider deprecating/removing the `'old'` option.
     character(len = 20) :: lorentz_scheme = 'default'
     !> Used as a part of the adaptive collisionality algorithm. When
     !> active we check the velocity space error with [[get_verr]]
     !> every `ncheck` time steps. This check can be relatively
     !> expensive so it is recommended to avoid small values of
     !> `ncheck`.
     !>
     !> @warning The new diagnostics module currently ignores this
     !> value and instead uses its own input variable named `ncheck`
     !> (which has a different default). See [this
     !> bug](https://bitbucket.org/gyrokinetics/gs2/issues/88).
     integer :: ncheck = 100
     !> If true (default) then potentially include the drag term in
     !> the pitch angle scattering operator. This is a necessary but
     !> not sufficient criteria. For the drag term to be included we
     !> also require \(\beta\neq 0\), more than one simulated species
     !> and finite \(A_\|\) perturbations included in the simulation
     !> (i.e. `fapar /= 0`).
     logical :: resistivity = .true.
     !> If true (not the default) then use special handling for the
     !> wfb particle in the pitch angle scattering operator.
     !>
     !> @note MRH changed default 16/08/2018. Previous default of true
     !> seemed to cause a numerical issue in flux tube simulations for
     !> the zonal modes at large \(k_x\).
     !>
     !> @todo Improve this documentation
     logical :: special_wfb_lorentz = .false.
     !> If true (not the default) then remove the collision operator
     !> from the usual time advance algorithm. Instead the collision
     !> operator is only applied every
     !> [[collisions_knobs:timesteps_between_collisions]]
     !> timesteps. This can potentially substantially speed up
     !> collisional simulations, both in the initialisation and
     !> advance phases.
     !>
     !> @warning Currently the input
     !> [[collisions_knobs:timesteps_between_collisions]] is ignored
     !> so collisions are applied every time step. The primary result
     !> of `split_collision = .true.` currently is that collisions are
     !> not applied in the first linear solve used as a part of a
     !> single time step. Hence the cost of collisions in advance are
     !> roughly halved. The full saving in the initialisation phase is
     !> still realised.
     logical :: split_collisions = .false.
     !> If true (not the default) then performs some additional checks
     !> of the data redistribution routines associated with
     !> transforming being the standard and collisional data
     !> decompositions.
     logical :: test = .false.
     !> Should set the number of timesteps between application of the
     !> collision operator if [[collisions_knobs:split_collisions]] is
     !> true. Currently this is ignored.
     !>
     !> @warning This value is currently ignored.
     integer :: timesteps_between_collisions = 1
     !> If true (default) then use a data decomposition for collisions
     !> that brings both pitch angle and energy local to a
     !> processor. This is typically the most efficient option,
     !> however for collisional simulations that only use part of the
     !> collision operator (either energy diffusion or pitch angle
     !> scattering) then it may be beneficial to set this flag to
     !> false such that we use a decomposition that only brings either
     !> energy or pitch angle local.
     logical :: use_le_layout = .true.
     !> Set to true (not the default) to activate the adaptive
     !> collisionality algorithm.
     !>
     !> @todo Provide more documentation on the adaptive
     !> collisionality algorithm.
     logical :: vary_vnew = .false.
     !> If the collisionality is to be increased as a part of the
     !> adaptive collisionality algorithm then increase it by this
     !> factor.
     real :: vnfac = 1.1
     !> If the collisionality is to be decreased as a part of the
     !> adaptive collisionality algorithm then decrease it by this
     !> factor.
     real :: vnslow = 0.9
     !> Controls how the duplicate `vpar = 0` point is handled.  When
     !> `vpar_zero_mean = .true.` (the default) the average of `g(vpar
     !> = 0)` for both signs of the parallel velcoity (`isgn`) is used
     !> in the collision operator instead of just `g(vpar = 0)` at
     !> `isgn=2`.  This is seen to suppress a numerical instability
     !> when `special_wfb_lorentz =.false.`. With these defaults pitch
     !> angle scattering at \(\theta = \pm \pi \) is now being handled
     !> physically i.e. `vpar = 0` at this theta location is no longer
     !> being skipped.
     !>
     !> @todo Consider removing this option.
     logical :: vpar_zero_mean = .true.
   contains
     procedure, public :: read => read_collisions_config
     procedure, public :: write => write_collisions_config
     procedure, public :: reset => reset_collisions_config
     procedure, public :: broadcast => broadcast_collisions_config
     procedure, public, nopass :: get_default_name => get_default_name_collisions_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_collisions_config
  end type collisions_config_type

  type(collisions_config_type) :: collisions_config
contains
  
  !> FIXME : Add documentation
  subroutine check_collisions(report_unit)
    implicit none
    integer, intent(in) :: report_unit
    select case (collision_model_switch)
    case (collision_model_lorentz,collision_model_lorentz_test)
       write (report_unit, fmt="('A Lorentz collision operator has been selected.')")
       if (cfac > 0) write (report_unit, fmt="('This has both terms of the Lorentz collision operator: cfac=',e12.4)") cfac
       if (cfac == 0) write (report_unit, fmt="('This is only a partial Lorentz collision operator (cfac=0.0)')")
       if (const_v) write (report_unit, fmt="('This is an energy independent Lorentz collision operator (const_v=true)')")  
!          if (hypercoll) call init_hyper_lorentz
    case (collision_model_full)
       write (report_unit, fmt="('Full GS2 collision operator has been selected.')")
    end select
  end subroutine check_collisions

  !> FIXME : Add documentation
  subroutine wnml_collisions(unit)
    implicit none
    integer, intent(in) :: unit
    if (.not.exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "collisions_knobs"
    select case (collision_model_switch)
    case (collision_model_lorentz)
       write (unit, fmt="(' collision_model = ',a)") '"lorentz"'
       if (hypermult) write (unit, fmt="(' hypermult = ',L1)") hypermult
    case (collision_model_lorentz_test)
       write (unit, fmt="(' collision_model = ',a)") '"lorentz-test"'
    case (collision_model_none)
       write (unit, fmt="(' collision_model = ',a)") '"collisionless"'
    end select
    write (unit, fmt="(' cfac = ',f6.3)") cfac
    write (unit, fmt="(' heating = ',L1)") heating
    write (unit, fmt="(' /')")
  end subroutine wnml_collisions

  !> FIXME : Add documentation
  subroutine init_collisions(collisions_config_in)
    use species, only: init_species, nspec, spec
    use theta_grid, only: init_theta_grid, ntgrid
    use kt_grids, only: init_kt_grids, naky, ntheta0
    use le_grids, only: init_le_grids, nlambda, negrid
    use run_parameters, only: init_run_parameters
    use gs2_layouts, only: init_dist_fn_layouts, init_gs2_layouts
    use mp, only: nproc, iproc
    implicit none
    type(collisions_config_type), intent(in), optional :: collisions_config_in

    if (initialized) return
    initialized = .true.
    call init_gs2_layouts
    call init_species

    hyper_colls = .false.
    if (any(spec%nu_h > epsilon(0.0))) hyper_colls = .true.

    call init_theta_grid
    call init_kt_grids
    call init_le_grids (accelerated_x, accelerated_v)
    call init_run_parameters
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nproc, iproc)
    call read_parameters(collisions_config_in)

    call init_arrays

  end subroutine init_collisions

  !> FIXME : Add documentation
  subroutine read_parameters(collisions_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use run_parameters, only: beta, fapar
    use species, only: nspec
    implicit none
    type(collisions_config_type), intent(in), optional :: collisions_config_in
    
    type (text_option), dimension (6), parameter :: modelopts = &
         (/ text_option('default', collision_model_full), &
            text_option('lorentz', collision_model_lorentz), &
            text_option('ediffuse', collision_model_ediffuse), &
            text_option('lorentz-test', collision_model_lorentz_test), &
            text_option('none', collision_model_none), &
            text_option('collisionless', collision_model_none) /)
    type (text_option), dimension (2), parameter :: schemeopts = &
         (/ text_option('default', lorentz_scheme_default), &
            text_option('old', lorentz_scheme_old) /)
    type (text_option), dimension (2), parameter :: eschemeopts = &
         (/ text_option('default', ediff_scheme_default), &
            text_option('old', ediff_scheme_old) /)
    character(20) :: collision_model, lorentz_scheme, ediff_scheme
    integer :: ierr

    if (present(collisions_config_in)) collisions_config = collisions_config_in

    call collisions_config%init(name = 'collisions_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    adjust = collisions_config%adjust
    cfac = collisions_config%cfac
    collision_model = collisions_config%collision_model
    conservative = collisions_config%conservative
    conserve_forbid_zero = collisions_config%conserve_forbid_zero
    conserve_moments = collisions_config%conserve_moments
    const_v = collisions_config%const_v
    ediff_scheme = collisions_config%ediff_scheme
    ei_coll_only = collisions_config%ei_coll_only
    etol = collisions_config%etol
    etola = collisions_config%etola
    ewindow = collisions_config%ewindow
    ewindowa = collisions_config%ewindowa
    force_collisions = collisions_config%force_collisions
    heating = collisions_config%heating
    hypermult = collisions_config%hypermult
    le_nxilim_fix = collisions_config%le_nxilim_fix
    lorentz_scheme = collisions_config%lorentz_scheme
    ncheck = collisions_config%ncheck
    resistivity = collisions_config%resistivity
    special_wfb_lorentz = collisions_config%special_wfb_lorentz
    split_collisions = collisions_config%split_collisions
    test = collisions_config%test
    timesteps_between_collisions = collisions_config%timesteps_between_collisions
    use_le_layout = collisions_config%use_le_layout
    vary_vnew = collisions_config%vary_vnew
    vnfac = collisions_config%vnfac
    vnslow = collisions_config%vnslow
    vpar_zero_mean = collisions_config%vpar_zero_mean

    exist = collisions_config%exist


    ierr = error_unit()
    call get_option_value &
         (collision_model, modelopts, collision_model_switch, &
         ierr, "collision_model in collisions_knobs",.true.)

    call get_option_value &
         (lorentz_scheme, schemeopts, lorentz_switch, &
         ierr, "lorentz_scheme in collisions_knobs",.true.)
    
    call get_option_value &
         (ediff_scheme, eschemeopts, ediff_switch, &
         ierr, "ediff_scheme in collisions_knobs",.true.)

    has_lorentz = .false. ; has_diffuse = .false.
    select case (collision_model_switch)
    case (collision_model_full)
       has_lorentz = .true. ; has_diffuse = .true.
    case (collision_model_lorentz,collision_model_lorentz_test)
       has_lorentz = .true.
    case (collision_model_ediffuse)
       has_diffuse = .true.
    end select
    
    drag = has_lorentz .and. resistivity .and. (beta > epsilon(0.0)) &
         .and. (nspec > 1) .and. (fapar.gt.0)
  end subroutine read_parameters

  !> A wrapper to sqrt which replaces -ve values with 0.0 to avoid
  !> NaNs arising from slight floating point discrepancies. We could
  !> consider adding a debug check to abort/warn if the passed value
  !> is too negative (i.e. if it looks like an error rather than small
  !> round off).
  pure real function safe_sqrt(arg)
    implicit none
    real, intent(in) :: arg
    safe_sqrt = sqrt(max(0.0, arg))
  end function safe_sqrt

  !> FIXME : Add documentation
  subroutine init_arrays
    use species, only: nspec
    use le_grids, only: negrid, init_map
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: c_rate
    use array_utils, only: zero_array
    implicit none
    real, dimension (negrid) :: hee
    logical :: use_lz_layout, use_e_layout
! lowflow terms include higher-order corrections to GK equation
! such as parallel nonlinearity that require derivatives in v-space.
! most efficient way to take these derivatives is to go from g_lo to le_lo,
! i.e., bring all energies and lambdas onto each processor
!<DD>Note the user can still disable use_le_layout in the input file
!    this just changes the default.
!>@note, not clear the above is true -- init_map will always be called with
!le layout when lowflow active here. Given use_le_layout = T is default we
!might be able to get rid of this conditional compilation now?
# ifdef LOWFLOW
    use_le_layout = .true.
# endif
    use_lz_layout = .false. ; use_e_layout = .false.

    !<DD>This (potentially large) array should only be needed if
    !write_heating (gs2_diagnostics) is true.
    if (heating .and. .not. allocated(c_rate)) then
       allocate (c_rate(-ntgrid:ntgrid, ntheta0, naky, nspec, 3))
       call zero_array(c_rate)
    end if

    !<DD>Should probably put this first
    if (collision_model_switch == collision_model_none) then
       colls = .false.
       return
    end if

    if( .not. use_le_layout ) then
       if (has_lorentz) use_lz_layout = .true.
       if (has_diffuse) use_e_layout = .true.
    end if

    call init_map (use_lz_layout, use_e_layout, use_le_layout, test)

    call init_vnew (hee)
    if (all(abs(vnew(:,1,:)) <= 2.0*epsilon(0.0)) .and. .not. force_collisions) then
       collision_model_switch = collision_model_none
       colls = .false.
       return
    end if

    if (has_lorentz) then
       call init_lorentz
       if (conserve_moments) call init_lorentz_conserve
    end if

    if (has_diffuse) then
       call init_ediffuse
       if (conserve_moments) call init_diffuse_conserve
    end if

    if (use_le_layout .and. (conserve_moments .or. drag)) call init_le_bessel
  end subroutine init_arrays

  !> Communicate Bessel functions from g_lo to le_lo
  subroutine init_le_bessel
    use gs2_layouts, only: g_lo, le_lo, ig_idx
    use dist_fn_arrays, only: aj0, aj1, vpa
    use le_grids, only: negrid, nxi, g2le, ixi_to_il, energy => energy_maxw, al
    use theta_grid, only: ntgrid
    use redistribute, only: gather, scatter
    use array_utils, only: zero_array
    implicit none
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
    integer :: ile, ig, ie, ixi, il

    allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    ! We need to initialise ctmp as it is used as receiving buffer in
    ! g2le redistribute, which doesn't populate all elements
    call zero_array(ctmp)

    ! next set aj0le & aj1l
    z_big(:,1,:) = cmplx(aj0,aj1)
    z_big(:,2,:) = z_big(:,1,:)

    call gather (g2le, z_big, ctmp)

    if (.not. allocated(aj0le)) then
       allocate (aj0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       allocate (vperp_aj1le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if

    aj0le = real(ctmp)
    vperp_aj1le = aimag(ctmp)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, ig, ixi, ie, il) &
    !$OMP SHARED(le_lo, negrid, nxi, ixi_to_il, vperp_aj1le, energy, al) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ig = ig_idx(le_lo, ile)
       do ie = 1, negrid
          do ixi = 1, nxi + 1
             il = ixi_to_il(ig, ixi)
             vperp_aj1le(ixi, ie, ile) = vperp_aj1le(ixi, ie, ile) * energy(ie) * al(il)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    z_big = vpa
    call gather (g2le, z_big, ctmp)
    deallocate(z_big)
    if (.not. allocated(vpa_aj0_le)) then
       allocate (vpa_aj0_le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    end if
    vpa_aj0_le = real(ctmp) * aj0le
  end subroutine init_le_bessel

  !> Precompute three quantities needed for momentum and energy conservation:
  !> z0, w0, s0 (z0le, w0le, s0le if use_le_layout chosen in the input file defined)
  subroutine init_lorentz_conserve
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec, is_electron_species
    use kt_grids, only: kperp2, naky, ntheta0, kwork_filter
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy => energy_maxw, speed => speed_maxw, al, &
         integrate_moment, negrid, forbid
    use gs2_time, only: code_dt, tunits
    use dist_fn_arrays, only: aj0, aj1, vpa
    use le_grids, only: g2le, nxi
    use gs2_layouts, only: le_lo
    use redistribute, only: gather, scatter
    use array_utils, only: zero_array
    implicit none
    complex, dimension (1,1,1) :: dum1 = 0., dum2 = 0.
    real, dimension (:,:,:), allocatable :: gtmp
    real, dimension (:,:,:,:), allocatable :: duinv, dtmp, vns
    integer :: ie, il, ik, is, isgn, iglo,  it, ig
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
    complex, dimension(:,:,:), allocatable :: s0tmp, w0tmp, z0tmp
    logical :: all = .true.

    if(use_le_layout) then
       allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       ! We need to initialise ctmp as it is used as receiving buffer in
       ! g2le redistribute, which doesn't populate all elements
       call zero_array(ctmp)
    end if

    allocate(s0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(w0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(z0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,3))

    !This initialisation is needed in case kwork_filter is true anywhere
    if (any(kwork_filter)) then
       call zero_array(duinv)
       call zero_array(dtmp)
    end if

    vns(:,:,:,1) = vnmult(1)*vnew_D
    vns(:,:,:,2) = vnmult(1)*vnew_s
    vns(:,:,:,3) = 0.0

    if (drag) then
       do is = 1, nspec
          if (.not. is_electron_species(spec(is))) cycle
          do ik = 1, naky
             vns(ik,:,is,3) = vnmult(1)*spec(is)%vnewk*tunits(ik)/energy**1.5
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get z0 (first form)
    if(conserve_forbid_zero) then ! MRH
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, conservative, z0tmp, code_dt, vns, vpdiff, speed, aj0, vpa) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u0 = -2 nu_D^{ei} vpa J0 dt f0
                if ( .not. forbid(ig,il) ) then
                   if (conservative) then
                      z0tmp(ig,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpdiff(ig,isgn,il) &
                           * speed(ie)*aj0(ig,iglo)
                   else
                      z0tmp(ig,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpa(ig,isgn,iglo)*aj0(ig,iglo)
                   end if
                else
                   z0tmp(ig,isgn,iglo) = 0.
                endif
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
       !$OMP SHARED(g_lo, conservative, z0tmp, code_dt, vns, vpdiff, speed, aj0, vpa) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u0 = -2 nu_D^{ei} vpa J0 dt f0
             if (conservative) then
                z0tmp(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpdiff(:,isgn,il) &
                     * speed(ie)*aj0(:,iglo)
             else
                z0tmp(:,isgn,iglo) = -2.0*code_dt*vns(ik,ie,is,3)*vpa(:,isgn,iglo)*aj0(:,iglo)
             end if
          end do
       end do
       !$OMP END PARALLEL DO
    endif

    call zero_out_passing_hybrid_electrons(z0tmp)

    if(use_le_layout) then
       call gather (g2le, z0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, z0tmp)   ! z0 is redefined below
    else
       call solfp_lorentz (z0tmp,dum1,dum2)   ! z0 is redefined below
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, isgn) &
    !$OMP SHARED(g_lo, vpa, gtmp, aj0, z0tmp) &
    !$OMP COLLAPSE(2) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo) * aj0(:,iglo) * real(z0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all) ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, z0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          z0tmp(:,isgn,iglo) = z0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    all = .true.

    ! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
    ! duinv = 1/du
    if (conservative) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
       !$OMP SHARED(g_lo, gtmp, vns, vpa, vpdiff, speed) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                  * vpdiff(:,isgn,il)*speed(ie)
             !             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*energy(ie)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
       !$OMP SHARED(g_lo, gtmp, vpa, vns) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo)**2
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet

    ! Could replace this with OpenMP using an explicit loop. TAG
    where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
       !duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get s0 (first form)
    if(conserve_forbid_zero) then ! MRH
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, conservative, vns, vpdiff, speed, aj0, code_dt, duinv, vpa, s0tmp) &       
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u1 = -3 nu_s vpa dt J0 f_0 / du
                if ( .not. forbid(ig,il) ) then
                   if (conservative) then
                      s0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*vpdiff(ig,isgn,il)*speed(ie) &
                           * aj0(ig,iglo)*code_dt*duinv(ig,it,ik,is)
                   else
                      s0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*vpa(ig,isgn,iglo) &
                           * aj0(ig,iglo)*code_dt*duinv(ig,it,ik,is)
                   end if
                else
                   s0tmp(ig,isgn,iglo) = 0.
                endif
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, conservative, s0tmp, vns, vpdiff, speed, aj0, code_dt, duinv, vpa) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u1 = -3 nu_s vpa dt J0 f_0 / du
             if (conservative) then
                s0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpdiff(:,isgn,il)*speed(ie) &
                     * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
             else
                s0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                     * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
             end if
          end do
       end do
       !$OMP END PARALLEL DO
    endif

    call zero_out_passing_hybrid_electrons(s0tmp)

    if(use_le_layout) then
       call gather (g2le, s0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, s0tmp)   ! s0
    else
       call solfp_lorentz (s0tmp,dum1,dum2)   ! s0
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, isgn) &
    !$OMP SHARED(g_lo, gtmp, vpa, aj0, s0tmp) &
    !$OMP COLLAPSE(2) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo) * aj0(:,iglo) * real(s0tmp(:,isgn,iglo))
       end do
    end do
    !OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v0s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 - v0s0 * z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, s0tmp, dtmp, z0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          s0tmp(:,isgn,iglo) = s0tmp(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * z0tmp(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1s0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, conservative, gtmp, vns, speed, vpdiff, aj0, s0tmp, vpa) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nu_D vpa J0
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*speed(ie)*vpdiff(:,isgn,il) &
                  * aj0(:,iglo) * real(s0tmp(:,isgn,iglo))
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * real(s0tmp(:,isgn,iglo))
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v1s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 / (1 + v0s0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, s0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          s0tmp(:,isgn,iglo) = s0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get w0
    if(conserve_forbid_zero) then ! MRH
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, conservative, w0tmp, vns, energy, al , aj1, code_dt, &
       !$OMP spec, kperp2, duinv, bmag) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u2 = -3 dt J1 vperp vus a f0 / du
                if ( .not. forbid(ig,il) ) then
                   if (conservative) then
                      w0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(ig,iglo) &
                           * code_dt*spec(is)%smz**2*kperp2(ig,it,ik)*duinv(ig,it,ik,is) &
                           / bmag(ig)
                   else
                      !This is identical to the other branch, should it be?
                      w0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(ig,iglo) &
                           * code_dt*spec(is)%smz**2*kperp2(ig,it,ik)*duinv(ig,it,ik,is) &
                           / bmag(ig)
                   end if
                else
                   w0tmp(ig,isgn,iglo) = 0.
                endif
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, conservative, w0tmp, vns, energy, al, aj1, code_dt, spec, kperp2, duinv, bmag) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u2 = -3 dt J1 vperp vus a f0 / du
             if (conservative) then
                w0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
                     * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                     / bmag
             else
                !This is identical to the other branch, should it be?
                w0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
                     * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                     / bmag
             end if
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    call zero_out_passing_hybrid_electrons(w0tmp)

    if(use_le_layout) then
       call gather (g2le, w0tmp, ctmp)
       call solfp_lorentz (ctmp)
       call scatter (g2le, ctmp, w0tmp)
    else
       call solfp_lorentz (w0tmp,dum1,dum2)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, isgn) &
    !$OMP SHARED(g_lo, gtmp, vpa, aj0, w0tmp) &
    !$OMP COLLAPSE(2) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          ! v0 = vpa J0 f0
          gtmp(:,isgn,iglo) = vpa(:,isgn,iglo) * aj0(:,iglo) * real(w0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v0w0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w0 - v0w0 * z0 / (1 + v0z0) (this is w1 from MAB notes)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, w0tmp, z0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) - z0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v1w1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, conservative, gtmp, vns, speed, vpdiff, aj0, w0tmp, vpa) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nud vpa J0 f0
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*speed(ie)*vpdiff(:,isgn,il) &
                  * aj0(:,iglo) * real(w0tmp(:,isgn,iglo))
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * real(w0tmp(:,isgn,iglo))
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v1w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w1 - v1w1 * s1 / (1 + v1s1) (this is w2 from MAB notes)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, w0tmp, s0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) - s0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v2w2

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vns, energy, al, aj1, w0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = nud vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
               * real(w0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)   ! v2w2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn, ig) &
    !$OMP SHARED(g_lo, w0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          w0tmp(:,isgn,iglo) = w0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate (gtmp, duinv, dtmp, vns)

    if(use_le_layout) then
       allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))

       ! first set s0le, w0le & z0le
       z_big = cmplx(real(s0tmp), real(w0tmp))

       ! get rid of z0, s0, w0 now that we've converted to z0le, s0le, w0le
       if (allocated(s0tmp)) deallocate(s0tmp)
       if (allocated(w0tmp)) deallocate(w0tmp)

       call gather (g2le, z_big, ctmp)

       if (.not. allocated(s0le)) then
          allocate (s0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (w0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       end if

       s0le = real(ctmp)
       w0le = aimag(ctmp)

       ! set z0le
       call gather (g2le, z0tmp, ctmp)
       if (allocated(z0tmp)) deallocate(z0tmp)
       if (.not. allocated(z0le)) allocate (z0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       z0le = real(ctmp)

       deallocate (ctmp, z_big)
    else
       !Only need the real components (imaginary part should be zero, just
       !use complex arrays to allow reuse of existing integrate routines etc.)
       if (.not. allocated(s0)) then
          allocate (s0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          s0=real(s0tmp)
          deallocate(s0tmp)
       endif

       if (.not. allocated(w0)) then
          allocate (w0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          w0=real(w0tmp)
          deallocate(w0tmp)
       endif

       if (.not. allocated(z0)) then
          allocate (z0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          z0=real(z0tmp)
          deallocate(z0tmp)
       endif
    end if

  end subroutine init_lorentz_conserve

  !> Precompute three quantities needed for momentum and energy conservation:
  !> bz0, bw0, bs0
  subroutine init_diffuse_conserve
    use gs2_layouts, only: g_lo, ie_idx, is_idx, ik_idx, il_idx, it_idx
    use species, only: nspec, spec
    use kt_grids, only: naky, ntheta0, kperp2
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy => energy_maxw, al, integrate_moment, negrid, forbid
    use gs2_time, only: code_dt
    use dist_fn_arrays, only: aj0, aj1, vpa
    use le_grids, only: g2le, nxi
    use gs2_layouts, only: le_lo
    use redistribute, only: gather, scatter
    use array_utils, only: zero_array
    implicit none
    real, dimension (:,:,:), allocatable :: gtmp
    real, dimension (:,:,:,:), allocatable :: duinv, dtmp, vns
    integer :: ie, il, ik, is, isgn, iglo, it, ig
    complex, dimension (:,:,:), allocatable :: ctmp, z_big
    complex, dimension (:,:,:), allocatable :: bs0tmp, bw0tmp, bz0tmp
    logical :: all

    if(use_le_layout) then
       allocate (ctmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       ! We need to initialise ctmp as it is used as receiving buffer in
       ! g2le redistribute, which doesn't populate all elements
       call zero_array(ctmp)
    end  if

    allocate(bs0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(bw0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate(bz0tmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (duinv(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (dtmp(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (vns(naky,negrid,nspec,2))

    ! Following might only be needed if any kwork_filter
    call zero_array(duinv)
    call zero_array(dtmp)

    vns(:,:,:,1) = vnmult(2)*delvnew
    vns(:,:,:,2) = vnmult(2)*vnew_s

    ! first obtain 1/du
    !    if (conservative) then

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, energy, vnmult, vnew_E) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          gtmp(:,isgn,iglo) = energy(ie)*vnmult(2)*vnew_E(ik,ie,is)
       end do
    end do
    !$OMP END PARALLEL DO

    all = .true.
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
!    else
!       do is = 1, nspec
!          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
!       end do
!    end if

    ! Could replace this with OpenMP using an explicit loop. TAG
    where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
       !  duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now get z0 (first form)
    if(conserve_forbid_zero) then ! MRH
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, bz0tmp, code_dt, vnmult, vnew_E, aj0, duinv) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u0 = -nu_E E dt J0 f_0 / du
                if ( .not. forbid(ig,il) ) then
                   bz0tmp(ig,isgn,iglo) = -code_dt*vnmult(2)*vnew_E(ik,ie,is) &
                        * aj0(ig,iglo)*duinv(ig,it,ik,is)
                else
                   bz0tmp(ig,isgn,iglo) = 0.
                endif
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, is, isgn) &
       !$OMP SHARED(g_lo, bz0tmp, code_dt, vnmult, vnew_E, aj0, duinv) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u0 = -nu_E E dt J0 f_0 / du
             bz0tmp(:,isgn,iglo) = -code_dt*vnmult(2)*vnew_E(ik,ie,is) &
                  * aj0(:,iglo)*duinv(:,it,ik,is)
          end do
       end do
       !$OMP END PARALLEL DO
    endif

    call zero_out_passing_hybrid_electrons(bz0tmp)

    if(use_le_layout) then
       call gather (g2le, bz0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bz0tmp)   ! bz0 is redefined below
    else
       call solfp_ediffuse_standard_layout (bz0tmp)   ! bz0 is redefined below
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0z0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vnmult, vnew_E, aj0, bz0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0 f_0
          gtmp(:,isgn,iglo) = vnmult(2) * vnew_E(ik,ie,is) * aj0(:,iglo) &
               * real(bz0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all) ! v0z0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine z0 = z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bz0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          bz0tmp(:,isgn,iglo) = bz0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! redefine dq = du (for momentum-conserving terms)
    ! du == int (E nu_s f_0);  du = du(z, kx, ky, s)
    ! duinv = 1/du
    !    if (conservative) then

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vns, vpa) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          gtmp(:,isgn,iglo)  = vns(ik,ie,is,1)*vpa(:,isgn,iglo)**2
          !             gtmp(:,isgn,iglo)  = vns(ik,ie,is,2)*energy(ie)

       end do
    end do
    !$OMP END PARALLEL DO

    all = .true.
    call integrate_moment (gtmp, duinv, all)  ! not 1/du yet
    !    else
    !       do is = 1, nspec
    !          duinv(:,:,:,is) = vnmult(2)*spec(is)%vnewk*sqrt(2./pi)
    !       end do
    !    end if


    ! Could replace this with OpenMP using an explicit loop. TAG
    where (abs(duinv) > epsilon(0.0))  ! necessary b/c some species may have vnewk=0
       !  duinv=0 iff vnew=0 so ok to keep duinv=0.
       duinv = 1./duinv  ! now it is 1/du
    end where

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get s0 (first form)
    if(conserve_forbid_zero) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, bs0tmp, vns, vpa, aj0, code_dt, duinv) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u1 = -3 nu_s vpa dt J0 f_0 / du
                !          if (conservative) then
                if ( .not. forbid(ig,il) ) then
                   bs0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*vpa(ig,isgn,iglo) &
                        * aj0(ig,iglo)*code_dt*duinv(ig,it,ik,is)
                else
                   bs0tmp(ig,isgn,iglo) = 0.
                endif
             end do
             !          else
             !             bs0(:,isgn,iglo) = -3.0*vns(ik,ie,is,2)*vpa(:,isgn,iglo) &
             !                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
             !          end if
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, is, isgn) &
       !$OMP SHARED(g_lo, bs0tmp, vns, vpa, aj0, code_dt, duinv) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u1 = -3 nu_s vpa dt J0 f_0 / du
             !          if (conservative) then
             bs0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*vpa(:,isgn,iglo) &
                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
             !          else
             !             bs0(:,isgn,iglo) = -3.0*vns(ik,ie,is,2)*vpa(:,isgn,iglo) &
             !                  * aj0(:,iglo)*code_dt*duinv(:,it,ik,is)
             !          end if
          end do
       end do
       !$OMP END PARALLEL DO
    endif

    call zero_out_passing_hybrid_electrons(bs0tmp)

    if(use_le_layout) then
       call gather (g2le, bs0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bs0tmp)   ! bs0
    else
       call solfp_ediffuse_standard_layout (bs0tmp)    ! s0
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0s0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vnmult, vnew_E, aj0, bs0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * real(bs0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v0s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 - v0s0 * z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bs0tmp, dtmp, bz0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bs0tmp(:,isgn,iglo) = bs0tmp(:,isgn,iglo) - dtmp(:,it,ik,is) &
               * bz0tmp(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1s0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vns, vpa, aj0, bs0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nu_s - nu_D) vpa J0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * real(bs0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v1s0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine s0 = s0 / (1 + v0s0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bs0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bs0tmp(:,isgn,iglo) = bs0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get w0
    if (conserve_forbid_zero) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, forbid, bw0tmp, vns, energy, al, aj1, code_dt, &
       !$OMP spec, kperp2, duinv, bmag) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid,ntgrid
                ! u0 = -3 dt J1 vperp vus a f0 / du
                if ( .not. forbid(ig,il) ) then
                   bw0tmp(ig,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(ig,iglo) &
                        * code_dt*spec(is)%smz**2*kperp2(ig,it,ik)*duinv(ig,it,ik,is) &
                        / bmag(ig)
                else
                   bw0tmp(ig,isgn,iglo) = 0.
                endif
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, it, ie, il, is, isgn) &
       !$OMP SHARED(g_lo, bw0tmp, vns, energy, al, aj1, code_dt, spec, kperp2, duinv, bmag) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             ! u0 = -3 dt J1 vperp vus a f0 / du
             bw0tmp(:,isgn,iglo) = -vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
                  * code_dt*spec(is)%smz**2*kperp2(:,it,ik)*duinv(:,it,ik,is) &
                  / bmag
          end do
       end do
       !$OMP END PARALLEL DO
    endif

    call zero_out_passing_hybrid_electrons(bw0tmp)

    if(use_le_layout) then
       call gather (g2le, bw0tmp, ctmp)
       call solfp_ediffuse_le_layout (ctmp)
       call scatter (g2le, ctmp, bw0tmp)
    else
       call solfp_ediffuse_standard_layout (bw0tmp)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v0w0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vnmult, vnew_E, aj0, bw0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * real(bw0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v0w0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w0 - v0w0 * z0 / (1 + v0z0) (this is w1 from MAB notes)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bw0tmp, bz0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) - bz0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v1w1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vns, vpa, aj0, bw0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nus-nud) vpa J0 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * real(bw0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)    ! v1w1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w1 - v1w1 * s1 / (1 + v1s1) (this is w2 from MAB notes)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bw0tmp, bs0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) - bs0tmp(:,isgn,iglo)*dtmp(:,it,ik,is)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get v2w2

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, gtmp, vns, energy, al, aj1, bw0tmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = (nus-nud) vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is,1)*energy(ie)*al(il)*aj1(:,iglo) &
               * real(bw0tmp(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, dtmp, all)   ! v2w2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Redefine w0 = w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, is, isgn) &
    !$OMP SHARED(g_lo, bw0tmp, dtmp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn=1,2
          bw0tmp(:,isgn,iglo) = bw0tmp(:,isgn,iglo) / (1.0 + dtmp(:,it,ik,is))
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate (gtmp, duinv, dtmp, vns)

    if(use_le_layout) then
       allocate (z_big(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))
       z_big=cmplx(real(bs0tmp),real(bw0tmp))
       deallocate (bs0tmp)
       deallocate (bw0tmp)
       call gather (g2le, z_big, ctmp)
       deallocate(z_big)
       if (.not. allocated(bs0le)) allocate (bs0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       if (.not. allocated(bw0le)) allocate (bw0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       bs0le = real(ctmp)
       bw0le = aimag(ctmp)

       call gather (g2le, bz0tmp, ctmp)
       deallocate (bz0tmp)
       if (.not. allocated(bz0le)) allocate (bz0le(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       bz0le = real(ctmp)
       deallocate (ctmp)
    else
       !Only need the real components (imaginary part should be zero, just
       !use complex arrays to allow reuse of existing integrate routines etc.)
       if (.not. allocated(bs0)) then
          allocate (bs0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bs0=real(bs0tmp)
          deallocate(bs0tmp)
       endif

       if (.not. allocated(bw0)) then
          allocate (bw0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bw0=real(bw0tmp)
          deallocate(bw0tmp)
       endif

       if (.not. allocated(bz0)) then
          allocate (bz0(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          bz0=real(bz0tmp)
          deallocate(bz0tmp)
       endif
    end if

  end subroutine init_diffuse_conserve
  
  !> If we have hybrid electrons then we need to remove the
  !> contribution to the conservation terms from the adiabatic
  !> passing part. This routine does this for us.
  subroutine zero_out_passing_hybrid_electrons(arr)
    use gs2_layouts, only: g_lo, ik_idx, il_idx, is_idx
    use species, only: has_hybrid_electron_species, spec
    use le_grids, only: is_passing_hybrid_electron
    implicit none
    complex, dimension(:, :, g_lo%llim_proc:), intent(in out) :: arr
    integer :: iglo, is, ik, il
    if (has_hybrid_electron_species(spec)) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, ik, il, is) &
       !$OMP SHARED(g_lo, arr) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ik = ik_idx(g_lo,iglo)
          il = il_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          if (is_passing_hybrid_electron(is, ik, il)) arr(:, :, iglo) = 0.0
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine zero_out_passing_hybrid_electrons

  !> FIXME : Add documentation
  subroutine init_vnew (hee)
    use species, only: nspec, spec, is_electron_species
    use le_grids, only: negrid, energy => energy_maxw, w => w_maxw
    use kt_grids, only: naky, ntheta0, kperp2
    use theta_grid, only: ntgrid
    use run_parameters, only: zeff, delt_option_switch, delt_option_auto
    use gs2_time, only: tunits
    use constants, only: pi
    use gs2_save, only: init_vnm
    use spfunc, only: erf => erf_ext
    real, dimension (:), intent (out) :: hee
    real,dimension (negrid)::heevth, hsg, hsgvth
    integer :: ik, ie, is, it, ig
    integer :: istatus
    real :: k4max
    real :: vl, vr, dv2l, dv2r
    real, dimension (2) :: vnm_saved

    if (delt_option_switch == delt_option_auto) then
       call init_vnm (vnm_saved, istatus)
       if (istatus == 0) vnm_init = vnm_saved
    endif

    ! Initialise vnmult from the values in run_parameters. This is
    ! either what has been restored from the restart file  if
    ! `delt_option == 'check_restart'` or 1 otherwise.
    if (all(vnmult == -1)) then
       vnmult = vnm_init
    end if

    do ie = 1, negrid
       hee(ie) = exp(-energy(ie))/sqrt(pi*energy(ie)) &
            + (1.0 - 0.5/energy(ie))*erf(sqrt(energy(ie)))

!>MAB
! hsg is the G of Hirshman and Sigmar
! added to allow for momentum conservation with energy diffusion
       hsg(ie) = hsg_func(sqrt(energy(ie)))
!<MAB
    end do

!heevth is hee but using only thermal velocity (energy independent)

    do ie = 1, negrid
       heevth(ie) = exp(-1.0)/sqrt(pi) + 0.5*erf(1.0)

!>MAB
! hsg is the G of Helander and Sigmar
! added to allow for momentum conservation with energy diffusion
       hsgvth(ie) = hsg_func(1.0)
!<MAB
    end do                                                                  

    if (.not.allocated(vnew)) then
       ! Here the ky dependence is just to allow us to scale by tunits.
       ! This saves some operations in later use at expense of more memory
       allocate (vnew(naky,negrid,nspec))
       allocate (vnew_s(naky,negrid,nspec))
       allocate (vnew_D(naky,negrid,nspec))
       allocate (vnew_E(naky,negrid,nspec))
       allocate (delvnew(naky,negrid,nspec))
    end if

    do is = 1, nspec
       ! Don't include following lines, as haven't yet accounted for possibility of timestep changing.
       ! Correct collision frequency if applying collisions only every nth timestep:
       !spec(is)%vnewk = spec(is)%vnewk * float(timesteps_between_collisions)
       if (is_electron_species(spec(is))) then
          do ie = 1, negrid
             do ik = 1, naky
                if (ei_coll_only) then
                   if (const_v) then
                      vnew(ik,ie,is) = spec(is)%vnewk &
                           *zeff*0.5*tunits(ik)
                   else
                      vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                           *zeff*0.5*tunits(ik)
                   endif

                   vnew_s(ik,ie,is)=0.
                   vnew_D(ik,ie,is)=0.
                   vnew_E(ik,ie,is)=0.
                   delvnew(ik,ie,is)=0.
                else
                   if (const_v) then
                      vnew(ik,ie,is) = spec(is)%vnewk &
                           *(zeff+heevth(ie))*0.5*tunits(ik)
                      vnew_s(ik,ie,is) = spec(is)%vnewk &
                           *hsgvth(ie)*4.0*tunits(ik)
                      vnew_D(ik,ie,is) = spec(is)%vnewk &
                           *heevth(ie)*tunits(ik)
                      if (.not. conservative) then
                         vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                              - 2.0*vnew_D(ik,ie,is)
                         delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                      end if
                   else
                      vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                           *(zeff + hee(ie))*0.5*tunits(ik)
                      vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(energy(ie)) &
                           *hsg(ie)*4.0*tunits(ik)
                      vnew_D(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                           *hee(ie)*tunits(ik)
                      if (.not. conservative) then
                         vnew_E(ik,ie,is) = energy(ie)*(vnew_s(ik,ie,is)*(2.0-0.5/energy(ie)) &
                              - 2.0*vnew_D(ik,ie,is))
                         delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                      end if
                   end if
                endif
             end do
          end do
       else
          do ie = 1, negrid
             do ik = 1, naky
                if (const_v) then
                   vnew(ik,ie,is) = spec(is)%vnewk &
                        *heevth(ie)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk &
                        *hsgvth(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = vnew_s(ik,ie,is)*1.5 &
                           - 2.0*vnew_D(ik,ie,is)
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                else
                   vnew(ik,ie,is) = spec(is)%vnewk/energy(ie)**1.5 &
                        *hee(ie)*0.5*tunits(ik)
                   vnew_s(ik,ie,is) = spec(is)%vnewk/sqrt(energy(ie)) &
                        *hsg(ie)*4.0*tunits(ik)
                   vnew_D(ik,ie,is) = 2.0*vnew(ik,ie,is)
                   if (.not. conservative) then
                      vnew_E(ik,ie,is) = energy(ie)*(vnew_s(ik,ie,is)*(2.0-0.5/energy(ie)) &
                           - 2.0*vnew_D(ik,ie,is))
                      delvnew(ik,ie,is) = vnew_s(ik,ie,is)-vnew_D(ik,ie,is)
                   end if
                end if
                if (ei_coll_only) then
                   vnew(ik,ie,is) = 0.
                   vnew_s(ik,ie,is)=0.
                   vnew_D(ik,ie,is)=0.
                   vnew_E(ik,ie,is)=0.
                   delvnew(ik,ie,is)=0.
                endif
             end do
          end do
       end if

       if (conservative) then

          do ie = 2, negrid-1
             vr = 0.5*(sqrt(energy(ie+1)) + sqrt(energy(ie)))
             vl = 0.5*(sqrt(energy(ie  )) + sqrt(energy(ie-1)))
             dv2r = (energy(ie+1) - energy(ie)) / (sqrt(energy(ie+1)) - sqrt(energy(ie)))
             dv2l = (energy(ie) - energy(ie-1)) / (sqrt(energy(ie)) - sqrt(energy(ie-1)))
          
             vnew_E(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*dv2l*hsg_func(vl) &
                  - vr*exp(-vr**2)*dv2r*hsg_func(vr)) / (sqrt(pi)*w(ie))
             delvnew(:,ie,is) = spec(is)%vnewk*tunits*(vl*exp(-vl**2)*hsg_func(vl) &
                     - vr*exp(-vr**2)*hsg_func(vr)) / (sqrt(pi*energy(ie))*w(ie))
          end do

          ! boundary at v = 0
          vr = 0.5*(sqrt(energy(2)) + sqrt(energy(1)))
          dv2r = (energy(2) - energy(1)) / (sqrt(energy(2)) - sqrt(energy(1)))

          vnew_E(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr)*dv2r &
               / (sqrt(pi)*w(1))
          delvnew(:,1,is) = -spec(is)%vnewk*tunits*vr*exp(-vr**2)*hsg_func(vr) &
               / (sqrt(pi*energy(1))*w(1))

          ! boundary at v -> infinity
          vl = 0.5*(sqrt(energy(negrid)) + sqrt(energy(negrid-1)))
          dv2l = (energy(negrid) - energy(negrid-1)) / (sqrt(energy(negrid)) - sqrt(energy(negrid-1)))

          vnew_E(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl)*dv2l &
               / (sqrt(pi)*w(negrid))
          delvnew(:,negrid,is) = spec(is)%vnewk*tunits*vl*exp(-vl**2)*hsg_func(vl) &
               / (sqrt(pi*energy(negrid))*w(negrid))

       end if

       ! add hyper-terms inside collision operator
!BD: Warning!
!BD: For finite magnetic shear, this is different in form from what appears in hyper.f90 
!BD: because kperp2 /= akx**2 + aky**2;  there are cross terms that are dropped in hyper.f90
!BD: Warning!
!BD: Also: there is no "grid_norm" option here and the exponent is fixed to 4 for now
       if (hyper_colls) then
          if(.not.allocated(vnewh)) allocate (vnewh(-ntgrid:ntgrid,ntheta0,naky,nspec))
          k4max = (maxval(kperp2))**2 
          do ik = 1, naky
             do it = 1, ntheta0
                do ig=-ntgrid,ntgrid
                   vnewh(ig,it,ik,is) = spec(is)%nu_h * kperp2(ig,it,ik)**2/k4max
                end do
             end do
          end do
       end if
    end do
    
  end subroutine init_vnew

  !> FIXME : Add documentation  
  function hsg_func (vel)

    use constants, only: pi
    use spfunc, only: erf => erf_ext

    implicit none

    real, intent (in) :: vel
    real :: hsg_func

    hsg_func = 0.5*erf(vel)/vel**2-exp(-vel**2)/(sqrt(pi)*vel)

  end function hsg_func

  !> FIXME : Add documentation
  subroutine init_ediffuse (vnmult_target)
    use le_grids, only: negrid, nxi, ng2
    use le_grids, only: forbid, ixi_to_il
    use egrid, only: zeroes => zeroes_maxw, x0 => x0_maxw
    use gs2_layouts, only: le_lo, e_lo, il_idx
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx
    use array_utils, only: zero_array

    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ie, is, ik, il, ig, it
    real, dimension (:), allocatable :: aa, bb, cc, xe
    integer :: ile, ixi
    integer :: ielo
    integer :: nxi_lim

    allocate (aa(negrid), bb(negrid), cc(negrid))
    allocate (xe(negrid))

    ! The nxi > 2 * ng2 check appears to be checking if we have
    ! trapped particles or not so could be replaced with grid_has_trapped_particles()
    if(le_nxilim_fix .and. nxi > 2*ng2 ) then
       nxi_lim = nxi + 1
    else
       nxi_lim = nxi
    endif

    ! want to use x variables instead of e because we want conservative form
    ! for the x-integration
    xe(1:negrid-1) = zeroes
    xe(negrid) = x0

    if (present(vnmult_target)) then
       vnmult(2) = max (vnmult_target, 1.0)
    end if


    if(use_le_layout) then

       if (.not.allocated(ec1le)) then
          allocate (ec1le   (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (ebetaale(nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (eqle    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          vnmult(2) = max(1.0, vnmult(2))
       end if
       call zero_array(ec1le)
       call zero_array(ebetaale)
       call zero_array(eqle)

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, ik, it, is, ig, ie, ixi, il, aa, bb, cc) &
       !$OMP SHARED(le_lo, nxi_lim, forbid, ec1le, ebetaale, negrid, eqle, ixi_to_il, xe) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ik = ik_idx(le_lo, ile)
          it = it_idx(le_lo, ile)
          is = is_idx(le_lo, ile)
          ig = ig_idx(le_lo, ile)
          do ixi=1, nxi_lim
             il = ixi_to_il(ig,ixi)
             if (forbid(ig,il)) cycle
             call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)
             ec1le(ixi,:,ile) = cc
             ebetaale(ixi,1,ile) = 1.0 / bb(1)
             do ie = 1, negrid-1
                eqle(ixi,ie+1,ile) = aa(ie+1) * ebetaale(ixi,ie,ile)
                ebetaale(ixi,ie+1,ile) = 1.0 / ( bb(ie+1) - eqle(ixi,ie+1,ile) &
                     * ec1le(ixi,ie,ile) )
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else

       if (.not.allocated(ec1)) then
          allocate (ec1    (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
          allocate (ebetaa (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
          allocate (eql    (negrid,e_lo%llim_proc:e_lo%ulim_alloc))
          vnmult(2) = max(1.0, vnmult(2))
       endif
       call zero_array(ec1)
       call zero_array(ebetaa)
       call zero_array(eql)

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ielo, ik, it, is, ig, ie, ixi, il, aa, bb, cc) &
       !$OMP SHARED(e_lo, forbid, xe, ec1, ebetaa, negrid, eql) &
       !$OMP SCHEDULE(static)
       do ielo = e_lo%llim_proc, e_lo%ulim_proc
          il = il_idx(e_lo, ielo)
          ig = ig_idx(e_lo, ielo)
          if (forbid(ig,il)) cycle
          is = is_idx(e_lo, ielo)
          ik = ik_idx(e_lo, ielo)
          it = it_idx(e_lo, ielo)
          call get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)
          ec1(:,ielo) = cc

          ! fill in the arrays for the tridiagonal
          ebetaa(1,ielo) = 1.0/bb(1)
          do ie = 1, negrid-1
             eql(ie+1,ielo) = aa(ie+1)*ebetaa(ie,ielo)
             ebetaa(ie+1,ielo) = 1.0/(bb(ie+1)-eql(ie+1,ielo)*ec1(ie,ielo))
          end do

       end do
       !$OMP END PARALLEL DO
    end if

    deallocate(aa, bb, cc, xe)

  end subroutine init_ediffuse

  !> FIXME : Add documentation
  subroutine get_ediffuse_matrix (aa, bb, cc, ig, ik, it, il, is, xe)
    use species, only: spec
    use theta_grid, only: bmag
    use le_grids, only: al, negrid, w=>w_maxw
    use spfunc, only: erf => erf_ext
    use constants, only: pi
    use kt_grids, only: kperp2
    use gs2_time, only: code_dt, tunits

    implicit none

    integer, intent (in) :: ig, ik, it, il, is
    real, dimension (:), intent (in) :: xe
    real, dimension (:), intent (out) :: aa, bb, cc

    integer :: ie
    real :: vn, slb1, xe0, xe1, xe2, xel, xer
    real :: capgr, capgl, ee

    vn = vnmult(2)*spec(is)%vnewk*tunits(ik)

    slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))     ! xi_j

    select case (ediff_switch)

    case (ediff_scheme_default)
       
       do ie = 2, negrid-1
          xe0 = xe(ie-1)
          xe1 = xe(ie)
          xe2 = xe(ie+1)
          
          xel = (xe0 + xe1)*0.5
          xer = (xe1 + xe2)*0.5
         
          capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/(xer*sqrt(pi))
          capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/(xel*sqrt(pi))
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -0.25*vn*code_dt*capgr/(w(ie)*(xe2 - xe1))
          aa(ie) = -0.25*vn*code_dt*capgl/(w(ie)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
             ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at v = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       capgr = 2.0*exp(-xer**2)*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer/sqrt(pi)
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -0.25*vn*code_dt*capgr/(w(1)*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!          end if

       ! boundary at v = infinity
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       capgl = 2.0*exp(-xel**2)*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel/sqrt(pi)

       ee = 0.125*(1.-slb1**2)*vnew_s(ik,negrid,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(negrid) = 0.0
       aa(negrid) = -0.25*vn*code_dt*capgl/(w(negrid)*(xe1 - xe0))
       bb(negrid) = 1.0 - aa(negrid) + ee*code_dt
       
!          if (heating) then
!             dd(negrid) =vnc*ee
!             hh(negrid) =vnh*ee
!          end if

    case (ediff_scheme_old)
       
       ! non-conservative scheme
       
       do ie = 2, negrid-1
          xe0 = xe(ie-1)
          xe1 = xe(ie)
          xe2 = xe(ie+1)
          
          xel = (xe0 + xe1)*0.5
          xer = (xe1 + xe2)*0.5
          
          capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
          capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
          
          ee = 0.125*(1.-slb1**2)*vnew_s(ik,ie,is) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(ie) = -vn*code_dt*capgr/((xer-xel)*(xe2 - xe1))
          aa(ie) = -vn*code_dt*capgl/((xer-xel)*(xe1 - xe0))
          bb(ie) = 1.0 - (aa(ie) + cc(ie)) + ee*code_dt
          
          ! coefficients for entropy heating calculation
!             if (heating) then
!                dd(ie) =vnc*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!                hh(ie) =vnh*(0.25*capgr/(w(ie)*(xe2-xe1)) + ee)
!             end if
       end do

       ! boundary at xe = 0
       xe1 = xe(1)
       xe2 = xe(2)
       
       xer = (xe1 + xe2)*0.5
       
       capgr = 0.5*exp(xe1**2-xer**2)/xe1**2*(erf(xer)-2.*xer*exp(-xer**2)/sqrt(pi))/xer
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,1,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -vn*code_dt*capgr/(xer*(xe2 - xe1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*code_dt
       
!          if (heating) then
!             dd(1) =vnc*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!             hh(1) =vnh*(0.25*capgr/(w(1)*(xe2-xe1)) + ee)
!          end if

       ! boundary at xe = 1
       xe0 = xe(negrid-1)
       xe1 = xe(negrid)
       
       xel = (xe1 + xe0)*0.5
       
       capgl = 0.5*exp(xe1**2-xel**2)/xe1**2*(erf(xel)-2.*xel*exp(-xel**2)/sqrt(pi))/xel
       
       ee = 0.125*(1.-slb1**2)*vnew_s(ik,negrid,is) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(negrid) = 0.0
       aa(negrid) = -vn*code_dt*capgl/((1.0-xel)*(xe1 - xe0))
       bb(negrid) = 1.0 - aa(negrid) + ee*code_dt

!          if (heating) then
!             dd(negrid) =vnc*ee
!             hh(negrid) =vnh*ee
!          end if

    end select

  end subroutine get_ediffuse_matrix

  !> FIXME : Add documentation
  subroutine init_lorentz (vnmult_target)
    use le_grids, only: negrid, jend, ng2, nxi, g2le, nlambda, al
    use le_grids, only: setup_trapped_lambda_grids_old_finite_difference, setup_passing_lambda_grids, &
         il_is_passing, il_is_wfb
    use gs2_layouts, only: le_lo
    use gs2_layouts, only: lz_lo, g_lo
    use gs2_layouts, only: ig_idx, ik_idx, ie_idx, is_idx, it_idx
    use theta_grid, only: ntgrid
    use species, only: has_hybrid_electron_species, spec
    use array_utils, only: zero_array
    implicit none

    real, intent (in), optional :: vnmult_target

    integer :: ig, il, it, ik, ie, is, je, te2
    real, dimension (:), allocatable :: aa, bb, cc
    real, dimension (:), allocatable :: dd, hh
    integer :: ile
    integer :: ilz
    complex, dimension (:,:,:), allocatable :: ctmp, z_big

    allocate (aa(nxi+1), bb(nxi+1), cc(nxi+1), dd(nxi+1), hh(nxi+1))

    if (.not. allocated(pitch_weights)) then
       ! Start by just copying the existing weights
       allocate(pitch_weights(-ntgrid:ntgrid, nlambda))
       pitch_weights = 0.0
       ! We have to recaculate the passing weights here as when Radau-Gauss
       ! grids are used the WFB (il=ng2+1) is treated as both passing and
       ! trapped. Otherwise we could just copy the passing weights from wl.
       call setup_passing_lambda_grids(al, pitch_weights)
       call setup_trapped_lambda_grids_old_finite_difference(al, pitch_weights)
    end if

    call init_vpdiff

    if(use_le_layout) then
       if (.not.allocated(c1le)) then
          allocate (c1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (betaale (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          allocate (qle     (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
          if (heating) then
             allocate (d1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
             allocate (h1le    (nxi+1, negrid, le_lo%llim_proc:le_lo%ulim_alloc))
             call zero_array(d1le)
             call zero_array(h1le)
          end if
          vnmult(1) = max(1.0, vnmult(1))
       endif
       call zero_array(c1le)
       call zero_array(betaale)
       call zero_array(qle)

       if (present(vnmult_target)) then
          vnmult(1) = max (vnmult_target, 1.0)
       end if


       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, ik, it, is, ig, je, te2, ie, aa, bb, cc, dd, hh, il) &
       !$OMP SHARED(le_lo, jend, special_wfb_lorentz, ng2, negrid, spec, c1le, d1le, h1le, qle, betaale) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ig = ig_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          is = is_idx(le_lo,ile)
          je = jend(ig)
          !if (je <= ng2+1) then
          ! MRH the above line is likely the original cause of the wfb bug in collisions
          ! by treating collisions arrays here as if there are only passing particles
          ! when there are only passing particles plus wfb (and no other trapped particles)
          ! introduced an inconsistency in the tri-diagonal solve for theta =+/- pi 
          ! Fixed below
          if (il_is_passing(je) .or. (il_is_wfb(je) .and. special_wfb_lorentz) ) then ! MRH
             te2 = 2*ng2
          else
             te2 = 2*je-1
          end if
          do ie = 1, negrid

             call get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)

             if (has_hybrid_electron_species(spec)) &
                  call set_hzero_lorentz_collisions_matrix(aa, bb, cc, je, ik, is)
             
             c1le(:,ie,ile) = cc
             if (allocated(d1le)) then
                d1le(:,ie,ile) = dd
                h1le(:,ie,ile) = hh
             end if

             qle(1,ie,ile) = 0.0
             betaale(1,ie,ile) = 1.0/bb(1)
             do il = 1, te2-1
                qle(il+1,ie,ile) = aa(il+1) * betaale(il,ie,ile)
                betaale(il+1,ie,ile) = 1.0 / ( bb(il+1) - qle(il+1,ie,ile) &
                     * c1le(il,ie,ile) )
             end do
             qle(te2+1,ie,ile) = 0.0
             betaale(te2+1,ie,ile) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    else

       if (.not.allocated(c1)) then
          allocate (c1(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
          allocate (betaa(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
          allocate (ql(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
          if (heating) then
             allocate (d1   (nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
             allocate (h1   (nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
             call zero_array(d1)
             call zero_array(h1)
          end if
          vnmult(1) = max(1.0, vnmult(1))
       end if

       call zero_array(c1)
       call zero_array(betaa)
       call zero_array(ql)

       if (present(vnmult_target)) then
          vnmult(1) = max (vnmult_target, 1.0)
       end if

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ilz, ik, it, is, ig, je, te2, ie, aa, bb, cc, dd, hh, il) &
       !$OMP SHARED(lz_lo, jend, special_wfb_lorentz, spec, c1, d1, h1, betaa, ql, ng2) &
       !$OMP SCHEDULE(static)
       do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
          is = is_idx(lz_lo,ilz)
          ik = ik_idx(lz_lo,ilz)
          it = it_idx(lz_lo,ilz)
          ie = ie_idx(lz_lo,ilz)
          ig = ig_idx(lz_lo,ilz)
          je = jend(ig)
          !if (je <= ng2+1) then
          ! MRH the above line is likely the original cause of the wfb bug in collisions
          ! by treating collisions arrays here as if there are only passing particles
          ! when there are only passing particles plus wfb (and no other trapped particles)
          ! introduced an inconsistency in the tri-diagonal solve for theta =+/- pi
          ! Fixed below
          if (il_is_passing(je) .or. (il_is_wfb(je) .and. special_wfb_lorentz) ) then ! MRH
             te2 = 2*ng2
          else
             te2 = 2*je-1
          end if

          call get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)

          if (has_hybrid_electron_species(spec)) &
               call set_hzero_lorentz_collisions_matrix(aa, bb, cc, je, ik, is)

          c1(:,ilz) = cc
          if (allocated(d1)) then
             d1(:,ilz) = dd
             h1(:,ilz) = hh
          end if

          betaa(1,ilz) = 1.0/bb(1)
          do il = 1, te2-1
             ql(il+1,ilz) = aa(il+1)*betaa(il,ilz)
             betaa(il+1,ilz) = 1.0/(bb(il+1)-ql(il+1,ilz)*c1(il,ilz))
          end do

          ql(1,ilz) = 0.0
          ql(te2+1:,ilz) = 0.0
          c1(te2+1:,ilz) = 0.0
          betaa(te2+1:,ilz) = 0.0
       end do
       !$OMP END PARALLEL DO
    end if

    deallocate (aa, bb, cc, dd, hh)

  end subroutine init_lorentz

  !> Special behaviour when h=0 for passing non-zonal electrons
  !>
  !> The effect of these changes is to exclude passing electrons
  !> From pitch angle scattering, and to enforce 0 passing as a
  !> boundary condition For the trapped particle pitch angle
  !> scattering.
  subroutine set_hzero_lorentz_collisions_matrix(aa, bb, cc, je, ik, is)
    use species, only: is_hybrid_electron_species, spec
    use kt_grids, only: aky
    use le_grids, only: ng2
    implicit none
    real, dimension(:), intent(in out) :: aa, bb, cc
    integer, intent(in) :: je, ik, is
    !> il index of 1st non-wfb trapped particle
    integer :: il_llim
    !> il index of last non-wfb trapped particle
    integer :: il_ulim

    il_llim = ng2 + 2
    il_ulim = 2*je-1 - (ng2+ 1)

    if ( is_hybrid_electron_species(spec(is)) .and. aky(ik) /= 0.) then
       aa(:il_llim) = 0.
       bb(:il_llim-1) = 1.
       cc(:il_llim-1) = 0.
       aa(il_ulim+1:) = 0.
       bb(il_ulim+1:) = 1.
       cc(il_ulim:) = 0.
    endif
  end subroutine set_hzero_lorentz_collisions_matrix

  !> FIXME : Add documentation
  subroutine get_lorentz_matrix (aa, bb, cc, dd, hh, ig, ik, it, ie, is)
    use species, only: spec
    use le_grids, only: al, energy => energy_maxw, xi, ng2
    use le_grids, only: jend, al, il_is_passing, il_is_wfb
    use gs2_time, only: code_dt, tunits
    use kt_grids, only: kperp2
    use theta_grid, only: bmag
    implicit none
    real, dimension (:), intent (out) :: aa, bb, cc, dd, hh
    integer, intent (in) :: ig, ik, it, ie, is
    integer :: il, je, te, te2, teh
    real :: slb0, slb1, slb2, slbl, slbr, vhyp, vn, vnh, vnc, ee, deltaxi

    je = jend(ig)
!
!CMR, 17/2/2014:
!         te, te2, teh: indices in xi, which runs from +1 -> -1.
!         te   :  index of minimum xi value >= 0.
!         te2  :  total #xi values = index of minimum xi value (= -1) 
!         teh  :  index of minimum xi value > 0.  
!                 teh = te if no bouncing particle at this location
!              OR teh = te-1 if there is a bouncing particle
!
    if (special_wfb_lorentz) then
!CMRDDGC, 17/2/2014:
!   This clause is appropriate for Lorentz collisons with 
!         SPECIAL (unphysical) treatment of wfb at its bounce point
       if (il_is_passing(je) .or. il_is_wfb(je)) then
          te = ng2
          te2 = 2*ng2
          teh = ng2
       else
          te = je
          te2 = 2*je-1
          teh = je-1
       end if
    else
!CMRDDGC, 17/2/2014:
!   This clause is appropriate for Lorentz collisons with 
!         STANDARD treatment of wfb at its bounce point
       if (il_is_passing(je)) then
          te = ng2
          te2 = 2*ng2
          teh = ng2
       else
          te = je
          te2 = 2*je-1
          teh = je-1
       end if
    endif 
    if (collision_model_switch == collision_model_lorentz_test) then
       vn = vnmult(1)*abs(spec(is)%vnewk)*tunits(ik)
       vnc = 0.
       vnh = 0.
    else
       if (hyper_colls) then
          vhyp = vnewh(ig, it, ik, is)
       else
          vhyp = 0.0
       end if
       ! vnc and vnh only needed when heating is true
       vnc = vnmult(1) * vnew(ik, ie, is)
       if (hypermult) then
          vn = vnmult(1) * vnew(ik, ie, is) * (1 + vhyp)
          vnh = vhyp * vnc
       else
          vn = vnmult(1) * vnew(ik, ie, is) + vhyp
          vnh = vhyp
       end if
    end if

    aa = 0.0 ; bb = 0.0 ; cc = 0.0 ; dd = 0.0 ; hh = 0.0

    select case (lorentz_switch)
       
    case (lorentz_scheme_default)          
       
       do il = 2, te-1
          slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
          slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
          slb2 = safe_sqrt(1.0 - bmag(ig)*al(il+1))
          
          slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
          slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
          
          ee = 0.5*energy(ie)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(il) = 2.0*vn*code_dt*(1.0 - slbr**2)/(pitch_weights(ig,il)*(slb2 - slb1))
          aa(il) = 2.0*vn*code_dt*(1.0 - slbl**2)/(pitch_weights(ig,il)*(slb1 - slb0))
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
          
          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*(-2.0*(1.0-slbr*slbr)/(pitch_weights(ig,il)*(slb2-slb1)) + ee)
             hh(il) =vnh*(-2.0*(1.0-slbr*slbr)/(pitch_weights(ig,il)*(slb2-slb1)) + ee)
          end if
       end do
       
       ! boundary at xi = 1
       slb1 = safe_sqrt(1.0-bmag(ig)*al(1))
       slb2 = safe_sqrt(1.0-bmag(ig)*al(2))
       
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = 2.0*vn*code_dt*(1.0 - slbr**2)/(pitch_weights(ig,1)*(slb2 - slb1))
       aa(1) = 0.0
       bb(1) = 1.0 - cc(1) + ee*vn*code_dt
       
       if (heating) then
          dd(1) =vnc*(-2.0*(1.0-slbr**2)/(pitch_weights(ig,1)*(slb2-slb1)) + ee)
          hh(1) =vnh*(-2.0*(1.0-slbr**2)/(pitch_weights(ig,1)*(slb2-slb1)) + ee)
       end if
       
       ! boundary at xi = 0
       il = te
       slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
       if (te == ng2) then
          slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
          slb2 = -slb1
       else
          slb1 = 0.0
          slb2 = -slb0
       end if
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
!
!CMR, 6/3/2014:  
! STANDARD treatment of pitch angle scattering must resolve T-P boundary.
! NEED special_wfb= .false. to resolve T-P boundary at wfb bounce point 
!     (special_wfb= .true. AVOIDS TP boundary at wfb bounce point)
!
! Original code (pre-r2766) used eq.(42) Barnes et al, Phys Plasmas 16, 072107 
! (2009), with pitch angle weights to enhance accuracy in derivatives.
! NB THIS FAILS at wfb bounce point, giving aa=cc=infinite, 
!    because weight wl(ig,il)=0 at wfb bounce point.
!    UNPHYSICAL as d/dxi ((1-xi^2)g) IS NOT resolved numerically for wl=0.
! MUST accept limitations of the grid resolution and USE FINITE coefficients! 
! FIX here by setting a FINITE width of the trapped region at wfb B-P 
!              deltaxi=xi(ig,ng2)-xi(ig,ng2+2)
! ASIDE: NB    deltaxi=wl is actually 2*spacing in xi !!!
!              which explains upfront factor 2 in definition of aa, cc
       deltaxi=pitch_weights(ig,il)
       if (te == je) deltaxi = 2.*deltaxi ! MRH vpar =0 fix 
       ! MRH appropriate when endpoint (te) is a vpar = 0 point (je)
       ! factor of 2 required above because xi=0 is a repeated point 
       ! on vpar > 0, vpar < 0 grids, and hence has half the weight associated 
       ! to it at each appearance compared to what the weight should be  
       ! as calculated in the continuum of points xi = (-1,1)
       if (.not. special_wfb_lorentz .and. deltaxi .eq. 0. .and. il_is_wfb(il)) then
          deltaxi=xi(ig,ng2)-xi(ig,ng2+2)
       endif
       cc(il) = 2.0*vn*code_dt*(1.0 - slbr**2)/(deltaxi*(slb2 - slb1))
       aa(il) = 2.0*vn*code_dt*(1.0 - slbl**2)/(deltaxi*(slb1 - slb0))
       bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
       
       if (heating) then
          dd(il) =vnc*(-2.0*(1.0-slbr**2)/(deltaxi*(slb2-slb1)) + ee)
          hh(il) =vnh*(-2.0*(1.0-slbr**2)/(deltaxi*(slb2-slb1)) + ee)
       end if
!CMRend

    case (lorentz_scheme_old)
       
       do il = 2, te-1
          slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
          slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
          slb2 = safe_sqrt(1.0 - bmag(ig)*al(il+1))
          
          slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
          slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
          
          ee = 0.5*energy(ie)*(1+slb1**2) &
               / (bmag(ig)*spec(is)%zstm)**2 &
               * kperp2(ig,it,ik)*cfac
          
          ! coefficients for tridiagonal matrix:
          cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
          aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
          
          ! coefficients for entropy heating calculation
          if (heating) then
             dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
             hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          end if
       end do
       
       ! boundary at xi = 1
       slb0 = 1.0
       slb1 = safe_sqrt(1.0-bmag(ig)*al(1))
       slb2 = safe_sqrt(1.0-bmag(ig)*al(2))
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(1) = -vn*code_dt*(-1.0 - slbr)/(slb2-slb1)
       aa(1) = 0.0
       bb(1) = 1.0 - (aa(1) + cc(1)) + ee*vn*code_dt
       
       if (heating) then
          dd(1) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(1) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
       end if
       
       ! boundary at xi = 0
       il = te
       slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
       if (te == ng2) then
          slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
          slb2 = -slb1
       else
          slb1 = 0.0
          slb2 = -slb0
       end if
       
       slbl = (slb1 + slb0)*0.5
       slbr = (slb1 + slb2)*0.5
       
       ee = 0.5*energy(ie)*(1+slb1**2) &
            / (bmag(ig)*spec(is)%zstm)**2 &
            * kperp2(ig,it,ik)*cfac
       
       cc(il) = -vn*code_dt*(1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
       aa(il) = -vn*code_dt*(1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
       bb(il) = 1.0 - (aa(il) + cc(il)) + ee*vn*code_dt
       
       if (heating) then
          dd(il) =vnc*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
          hh(il) =vnh*((1.0-slbr*slbr)/(slbr-slbl)/(slb2-slb1) + ee)
       end if
       
    end select

    ! assuming symmetry in xi, fill in the rest of the arrays.
    aa(te+1:te2) = cc(teh:1:-1)
    bb(te+1:te2) = bb(teh:1:-1)
    cc(te+1:te2) = aa(teh:1:-1)
       
    if (heating) then
       dd(te+1:te2) = dd(teh:1:-1)
       hh(te+1:te2) = hh(teh:1:-1)
    end if

  end subroutine get_lorentz_matrix

  !> FIXME : Add documentation  
  subroutine init_lorentz_error

    use le_grids, only: jend, al, ng2, nlambda, &
         il_is_wfb, il_is_passing, il_is_trapped
    use theta_grid, only: ntgrid, bmag
    use array_utils, only: zero_array
    implicit none
    
    integer :: je, ig, il, ip, ij, im
    real :: slb0, slb1, slb2, slbr, slbl
    real, dimension (:), allocatable :: slb
    real, dimension (:,:), allocatable :: dprod
    real, dimension (:,:,:), allocatable :: dlcoef, d2lcoef

    allocate(slb(2*nlambda))
    allocate (dprod(nlambda,5))

    allocate (dlcoef(-ntgrid:ntgrid,nlambda,5))
    allocate (d2lcoef(-ntgrid:ntgrid,nlambda,5))
    allocate (dtot(-ntgrid:ntgrid,nlambda,5))
    allocate (fdf(-ntgrid:ntgrid,nlambda), fdb(-ntgrid:ntgrid,nlambda))

    dlcoef = 1.0; call zero_array(d2lcoef); call zero_array(dtot)
    call zero_array(fdf); call zero_array(fdb); slb = 0.0

    do ig=-ntgrid,ntgrid
       je = jend(ig)
       
       if (il_is_passing(je) .or. il_is_wfb(je)) then            ! no trapped particles

! calculation of xi and finite difference coefficients for non-boundary points
          do il=2,ng2-1
             slb(il) = safe_sqrt(1.0-al(il)*bmag(ig))   ! xi_{j}
             
             slb2 = safe_sqrt(1.0-al(il+1)*bmag(ig))    ! xi_{j+1}
             slb1 = slb(il)
             slb0 = safe_sqrt(1.0-al(il-1)*bmag(ig))    ! xi_{j-1}
             
             slbr = (slb2+slb1)*0.5                     ! xi_{j+1/2}
             slbl = (slb1+slb0)*0.5                     ! xi_{j-1/2}

! finite difference coefficients
             fdf(ig,il) = (1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             fdb(ig,il) = (1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          end do

! boundary at xi = 1
          slb(1) = safe_sqrt(1.0-al(1)*bmag(ig))
          slb0 = 1.0
          slb1 = slb(1)
          slb2 = slb(2)

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

! derivative of [(1-xi**2)*df/dxi] at xi_{j=1} is centered, with upper xi=1 and
! lower xi = xi_{j+1/2}
          fdf(ig,1) = (-1.0-slbr)/(slb2-slb1)
          fdb(ig,1) = 0.0

! boundary at xi = 0
          il = ng2
          slb(il) = safe_sqrt(1.0 - al(il)*bmag(ig))
          slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
          slb1 = slb(il)
          slb2 = -slb1

          slbl = (slb1 + slb0)/2.0
          slbr = (slb1 + slb2)/2.0

          fdf(ig,il) = (1.0 - slbr*slbr)/(slbr-slbl)/(slb2-slb1)
          fdb(ig,il) = (1.0 - slbl*slbl)/(slbr-slbl)/(slb1-slb0)

          slb(ng2+1:) = -slb(ng2:1:-1)
       else          ! run with trapped particles
          do il=2,je-1
             slb(il) = safe_sqrt(1.0-al(il)*bmag(ig))
             
             slb2 = safe_sqrt(1.0-al(il+1)*bmag(ig))
             slb1 = slb(il)
             slb0 = safe_sqrt(1.0-al(il-1)*bmag(ig))
             
             slbr = (slb2+slb1)*0.5
             slbl = (slb1+slb0)*0.5

             fdf(ig,il) = (1.0 - slbr*slbr)/(slbr - slbl)/(slb2 - slb1)
             fdb(ig,il) = (1.0 - slbl*slbl)/(slbr - slbl)/(slb1 - slb0)
          end do

! boundary at xi = 1
          slb(1) = safe_sqrt(1.0-bmag(ig)*al(1))
          slb0 = 1.0
          slb1 = slb(1)
          slb2 = slb(2)

          slbr = (slb1 + slb2)/2.0

          fdf(ig,1) = (-1.0 - slbr)/(slb2-slb1)
          fdb(ig,1) = 0.0

! boundary at xi = 0
          il = je
          slb(il) = safe_sqrt(1.0-bmag(ig)*al(il))
          slb0 = slb(je-1)
          slb1 = 0.
          slb2 = -slb0                                                        

          slbl = (slb1 + slb0)/2.0

          fdf(ig,il) = (1.0 - slbl*slbl)/slb0/slb0
          fdb(ig,il) = fdf(ig,il)

          slb(je+1:2*je-1) = -slb(je-1:1:-1)
       end if

! compute coefficients (dlcoef) multipyling first derivative of h
       do il=3,ng2
          do ip=il-2,il+2
             if (il == ip) then
                dlcoef(ig,il,ip-il+3) = 0.0
                do ij=il-2,il+2
                   if (ij /= ip) dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3) + 1/(slb(il)-slb(ij))
                end do
             else
                do ij=il-2,il+2
                   if (ij /= ip .and. ij /= il) then
                      dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do
                dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
             end if
             dlcoef(ig,il,ip-il+3) = -2.0*slb(il)*dlcoef(ig,il,ip-il+3)
          end do
       end do

       il = 1
       do ip=il,il+2
          if (il == ip) then
             dlcoef(ig,il,ip) = 0.0
             do ij=il,il+2
                if (ij /= ip) dlcoef(ig,il,ip) = dlcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
             end do
          else
             do ij=il,il+2
                if (ij /= ip .and. ij /= il) then
                   dlcoef(ig,il,ip) = dlcoef(ig,il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do
             dlcoef(ig,il,ip) = dlcoef(ig,il,ip)/(slb(ip)-slb(il))
          end if
          dlcoef(ig,il,ip) = -2.0*slb(il)*dlcoef(ig,il,ip)
       end do

       il = 2
       do ip=il-1,il+1
          if (il == ip) then
             dlcoef(ig,il,ip-il+2) = 0.0
             do ij=il-1,il+1
                if (ij /= ip) dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2) + 1/(slb(il)-slb(ij))
             end do
          else
             do ij=il-1,il+1
                if (ij /= ip .and. ij /= il) then
                   dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do
             dlcoef(ig,il,ip-il+2) = dlcoef(ig,il,ip-il+2)/(slb(ip)-slb(il))
          end if
          dlcoef(ig,il,ip-il+2) = -2.0*slb(il)*dlcoef(ig,il,ip-il+2)
       end do

       dprod = 2.0

! compute coefficients (d2lcoef) multiplying second derivative of h
       do il=3,ng2
          do ip=il-2,il+2
             if (il == ip) then
                do ij=il-2,il+2
                   if (ij /= ip) then
                      do im=il-2,il+2
                         if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+3) = &
                              d2lcoef(ig,il,ip-il+3) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                      end do
                   end if
                end do
             else
                do ij=il-2,il+2
                   if (ij /= il .and. ij /= ip) then
                      dprod(il,ip-il+3) = dprod(il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                   end if
                end do

                do ij=il-2,il+2
                   if (ij /= ip .and. ij /= il) then
                      d2lcoef(ig,il,ip-il+3) = d2lcoef(ig,il,ip-il+3) + 1./(slb(il)-slb(ij))
                   end if
                end do
                d2lcoef(ig,il,ip-il+3) = dprod(il,ip-il+3) &
                     *d2lcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
             end if
             d2lcoef(ig,il,ip-il+3) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+3)
          end do
       end do

       il = 1
       do ip=il,il+2
          if (il == ip) then
             do ij=il,il+2
                if (ij /= ip) then
                   do im=il,il+2
                      if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                   end do
                end if
             end do
          else
             do ij=il,il+2
                if (ij /= il .and. ij /= ip) then
                   dprod(il,ip) = dprod(il,ip)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do

             do ij=il,il+2
                if (ij /= ip .and. ij /= il) then
                   d2lcoef(ig,il,ip) = d2lcoef(ig,il,ip) + 1./(slb(il)-slb(ij))
                end if
             end do
             d2lcoef(ig,il,ip) = dprod(il,ip)*d2lcoef(ig,il,ip)/(slb(ip)-slb(il))
          end if
          d2lcoef(ig,il,ip) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip)
       end do

       il = 2
       do ip=il-1,il+1
          if (il == ip) then
             do ij=il-1,il+1
                if (ij /= ip) then
                   do im=il-1,il+1
                      if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+2) &
                           = d2lcoef(ig,il,ip-il+2) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                   end do
                end if
             end do
          else
             do ij=il-1,il+1
                if (ij /= il .and. ij /= ip) then
                   dprod(il,ip-il+2) = dprod(il,ip-il+2)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                end if
             end do

             do ij=il-1,il+1
                if (ij /= ip .and. ij /= il) then
                   d2lcoef(ig,il,ip-il+2) = d2lcoef(ig,il,ip-il+2) + 1./(slb(il)-slb(ij))
                end if
             end do
             d2lcoef(ig,il,ip-il+2) = dprod(il,ip-il+2)*d2lcoef(ig,il,ip-il+2)/(slb(ip)-slb(il))
          end if
          d2lcoef(ig,il,ip-il+2) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+2)
       end do
       
       if (il_is_trapped(je)) then      ! have to handle trapped particles

          do il=ng2+1,je
             do ip=il-2,il+2
                if (il == ip) then
                   dlcoef(ig,il,ip-il+3) = 0.0
                   do ij=il-2,il+2
                      if (ij /= ip) dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3) + 1/(slb(il)-slb(ij))
                   end do
                else
                   do ij=il-2,il+2
                      if (ij /= ip .and. ij /= il) then
                         dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                      end if
                   end do
                   dlcoef(ig,il,ip-il+3) = dlcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
                end if
                dlcoef(ig,il,ip-il+3) = -2.0*slb(il)*dlcoef(ig,il,ip-il+3)
             end do
          end do

          do il=ng2+1,je
             do ip=il-2,il+2
                if (il == ip) then
                   do ij=il-2,il+2
                      if (ij /= ip) then
                         do im=il-2,il+2
                            if (im /= ip .and. im /= ij) d2lcoef(ig,il,ip-il+3) = &
                                 d2lcoef(ig,il,ip-il+3) + 1./((slb(il)-slb(im))*(slb(il)-slb(ij)))
                         end do
                      end if
                   end do
                else
                   do ij=il-2,il+2
                      if (ij /= il .and. ij /= ip) then
                         dprod(il,ip-il+3) = dprod(il,ip-il+3)*(slb(il)-slb(ij))/(slb(ip)-slb(ij))
                      end if
                   end do
                   
                   do ij=il-2,il+2
                      if (ij /= ip .and. ij /= il) then
                         d2lcoef(ig,il,ip-il+3) = d2lcoef(ig,il,ip-il+3) + 1./(slb(il)-slb(ij))
                      end if
                   end do
                   d2lcoef(ig,il,ip-il+3) = dprod(il,ip-il+3) &
                        *d2lcoef(ig,il,ip-il+3)/(slb(ip)-slb(il))
                end if
                d2lcoef(ig,il,ip-il+3) = (1.0-slb(il)**2)*d2lcoef(ig,il,ip-il+3)
             end do
          end do
          
       end if
    end do

    dtot = dlcoef + d2lcoef

    deallocate (slb, dprod, dlcoef, d2lcoef)
  end subroutine init_lorentz_error

  !> FIXME : Add documentation  
  subroutine solfp1_standard_layout (g, g1, gc1, gc2, diagnostics, gtoc, ctog)
    use gs2_layouts, only: g_lo, it_idx, ik_idx, ie_idx, is_idx
    use theta_grid, only: ntgrid
    use run_parameters, only: beta
    use gs2_time, only: code_dt
    use le_grids, only: energy => energy_maxw
    use species, only: spec, is_electron_species
    use dist_fn_arrays, only: vpa, aj0
    use fields_arrays, only: aparnew
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter, kperp2
    use optionals, only: get_option_with_default
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1, gc1, gc2
    integer, optional, intent (in) :: diagnostics
    logical, optional, intent (in) :: gtoc, ctog
    integer :: isgn, it, ik, ie, is, iglo
!CMR, 12/9/2013: 
!CMR   New logical optional input parameters gtoc, ctog used to set
!CMR   flags (g_to_c and c_to_g) to control whether redistributes required
!CMR   to map g_lo to collision_lo, and collision_lo to g_lo.
!CMR   All redistributes are performed by default.
!CMR  
    logical :: g_to_c, c_to_g

    g_to_c = get_option_with_default(gtoc, .true.)
    c_to_g = get_option_with_default(ctog, .true.)

    if (has_diffuse) then
       call solfp_ediffuse_standard_layout (g)
       if (conserve_moments) call conserve_diffuse (g, g1)
    end if

    if (has_lorentz) then
       if (drag) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(iglo, is, it, ik, ie, isgn) &
          !$OMP SHARED(g_lo, spec, kwork_filter, g, ieqzip, vnmult, code_dt, vpa, &
          !$OMP kperp2, aparnew, aj0, beta, energy) &
          !$OMP SCHEDULE(static)
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             is = is_idx(g_lo,iglo)
             if (.not. is_electron_species(spec(is))) cycle
             it = it_idx(g_lo,iglo)
             ik = ik_idx(g_lo,iglo)
             if(kwork_filter(it,ik)) cycle
             if(ieqzip(it,ik)) cycle
             ie = ie_idx(g_lo,iglo)
             do isgn = 1, 2
                g(:, isgn, iglo) = g(:, isgn, iglo) + vnmult(1)*spec(is)%vnewk*code_dt &
                     * vpa(:,isgn,iglo)*kperp2(:,it,ik)*aparnew(:,it,ik)*aj0(:,iglo) &
                     / ((-spec(is)%z*spec(is)%dens)*beta*spec(is)%stm*energy(ie)**1.5)
                ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
                ! This has been implemented as 1/-spec(electron)%z*spec(electron)%dens
                ! in an attempt handle the multi-ion species case.
             end do
          end do
          !$OMP END PARALLEL DO
       end if

       call solfp_lorentz (g, gc1, gc2, diagnostics)
       if (conserve_moments) call conserve_lorentz (g, g1)
    end if
  end subroutine solfp1_standard_layout

  !> FIXME : Add documentation  
  subroutine solfp1_le_layout (gle, diagnostics)
    use gs2_layouts, only: le_lo, it_idx, ik_idx, ig_idx, is_idx
    use run_parameters, only: beta, ieqzip
    use gs2_time, only: code_dt
    use le_grids, only: energy => energy_maxw, negrid, nxi, ixi_to_il, ixi_to_isgn, ng2
    use species, only: spec, is_electron_species
    use fields_arrays, only: aparnew
    use kt_grids, only: kwork_filter, kperp2
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    integer, optional, intent (in) :: diagnostics
    complex :: tmp
    integer :: ig, it, ik, il, ie, is, ile, ixi, isgn
    integer :: nxi_lim

    ! The nxi > 2 * ng2 check appears to be checking if we have
    ! trapped particles or not so could be replaced with grid_has_trapped_particles()
    if(le_nxilim_fix .and. nxi > 2*ng2 ) then
        nxi_lim = nxi + 1
    else
        nxi_lim = nxi
    endif

    if (has_diffuse) then
       call solfp_ediffuse_le_layout (gle)
       if (conserve_moments) call conserve_diffuse (gle)
    end if

    if (has_lorentz) then
       if (drag) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ile, is, it, ik, ie, ig, ixi, tmp) &
          !$OMP SHARED(le_lo, spec, kwork_filter, negrid, ieqzip, vnmult, code_dt, kperp2, &
          !$OMP aparnew, beta, energy, nxi_lim, gle, vpa_aj0_le) &
          !$OMP SCHEDULE(static)
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             is = is_idx(le_lo,ile)
             if (.not. is_electron_species(spec(is))) cycle
             it = it_idx(le_lo,ile)
             ik = ik_idx(le_lo,ile)
             if(kwork_filter(it,ik)) cycle
             if(ieqzip(it,ik)) cycle
             ig = ig_idx(le_lo,ile)
             do ie = 1, negrid
                tmp = vnmult(1)*spec(is)%vnewk*code_dt &
                     * kperp2(ig,it,ik)*aparnew(ig,it,ik) &
                     / ((-spec(is)%z*spec(is)%dens)*beta*spec(is)%stm*energy(ie)**1.5)
                do ixi = 1, nxi_lim
                   gle(ixi, ie, ile) = gle(ixi, ie, ile) + tmp * vpa_aj0_le(ixi, ie, ile)
                ! probably need 1/(spec(is_ion)%z*spec(is_ion)%dens) above
                ! This has been implemented as 1/-spec(electron)%z*spec(electron)%dens
                ! in an attempt handle the multi-ion species case.
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if

       call solfp_lorentz (gle, diagnostics)
       if (conserve_moments) call conserve_lorentz (gle)
    end if
  end subroutine solfp1_le_layout

  !> FIXME : Add documentation
  subroutine conserve_lorentz_standard_layout (g, g1)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy => energy_maxw, speed => speed_maxw, al, &
         integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    use array_utils, only: copy
    
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2

    integer :: isgn, iglo, ik, ie, il, is, it
    logical :: all = .true. 

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (vns(naky,negrid,nspec))
    vns = vnmult(1)*vnew_D

    if (drag) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, it, ik, isgn) &
       !$OMP SHARED(g_lo, kwork_filter, gtmp, vpa, aj0, g) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          if(kwork_filter(it,ik))cycle
          do isgn = 1, 2
             ! v0 = vpa J0 f0, y0 = g
             gtmp(:,isgn,iglo) = vpa(:,isgn,iglo)*aj0(:,iglo)*g(:,isgn,iglo)
          end do
       end do
       !$OMP END PARALLEL DO
       call integrate_moment (gtmp, v0y0, all)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, it, ik, is, isgn) &
       !$OMP SHARED(g_lo, kwork_filter, g1, g, ieqzip, v0y0, z0) &
       !$OMP SCHEDULE(static)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          if(kwork_filter(it,ik)) cycle
          if(ieqzip(it,ik)) cycle
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g1(:,isgn,iglo) = g(:,isgn,iglo) - v0y0(:,it,ik,is) * z0(:,isgn,iglo)
          end do
       end do
       !$OMP END PARALLEL DO
    else

       call copy(g, g1)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, conservative, gtmp, vns, speed, vpdiff, aj0, g1, vpa) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = nud vpa J0 f0, y1 = g1
          if (conservative) then
             gtmp(:,isgn,iglo) = vns(ik,ie,is)*speed(ie)*vpdiff(:,isgn,il) &
                  * aj0(:,iglo)*g1(:,isgn,iglo)
          else
             gtmp(:,isgn,iglo) = vns(ik,ie,is)*vpa(:,isgn,iglo)*aj0(:,iglo) &
                  * g1(:,isgn,iglo)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, v1y1, all)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, g1, ieqzip, v1y1, s0) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g1(:,isgn,iglo) - v1y1(:,it,ik,is) * s0(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, gtmp, vns, energy, al, aj1, g1) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = nud vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*energy(ie)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, g, g1, ieqzip, v2y2, w0) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - v2y2(:,it,ik,is) * w0(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate (vns, v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_lorentz_standard_layout

  !> FIXME : Add documentation
  subroutine conserve_lorentz_le_layout (gle)
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: speed => speed_maxw, nxi, ng2
    use le_grids, only: ixi_to_il, ixi_to_isgn, integrate_moment, negrid
    use run_parameters, only: ieqzip
    use gs2_layouts, only: le_lo, ig_idx
    use redistribute, only: scatter
    use kt_grids, only: kwork_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    complex, dimension (:,:,:), allocatable :: gtmp

    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2
    real :: vpa_nud, nud
    integer :: ig, isgn, ik, ie, il, is, it
    integer :: ile, ixi
    integer :: nxi_lim

    ! The nxi > 2 * ng2 check appears to be checking if we have
    ! trapped particles or not so could be replaced with grid_has_trapped_particles()
    if(le_nxilim_fix .and. nxi > 2*ng2 ) then
        nxi_lim = nxi + 1
    else
        nxi_lim = nxi
    endif

    ! We only need one of these at a time so could just reuse one allocation.
    ! We'll keep two for now to help with code clarity and note that these
    ! should not be large arrays.
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    ! Let's work on gle directly instead of g for the moment
    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
    call zero_array(gtmp)

    if (drag) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

       allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))

       ! v0 = vpa J0 f0, y0 = gle
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, it, ik, ixi, ie) &
       !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gtmp, vpa_aj0_le, gle) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          if(kwork_filter(it, ik)) cycle
          do ie = 1, negrid
             do ixi = 1, nxi_lim
                ! Should we use vpdiff here if conservative?
                gtmp(ixi, ie, ile) =  vpa_aj0_le(ixi, ie, ile) * gle(ixi, ie, ile)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, it, ik, ixi, ie) &
       !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, z0le, v0y0) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          if(kwork_filter(it,ik)) cycle
          if(ieqzip(it,ik)) cycle
          do ie = 1, negrid
             do ixi = 1, nxi_lim
                gle(ixi,ie,ile) = gle(ixi,ie,ile) - z0le(ixi,ie,ile) * v0y0(ile)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       deallocate(v0y0)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    ! v1 = nud vpa J0 f0, y1 = gle
    if (conservative) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, is, it, ik, ig, ixi, il, ie, isgn, vpa_nud) &
       !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, ixi_to_il, ixi_to_isgn, vpdiff, &
       !$OMP speed, vnmult, vnew_D, gtmp, aj0le, gle) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ik = ik_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          if(kwork_filter(it,ik)) cycle
          is = is_idx(le_lo,ile)
          ig = ig_idx(le_lo,ile)
          do ie = 1, negrid
             do ixi = 1, nxi_lim
                il = ixi_to_il(ig,ixi)
                isgn = ixi_to_isgn(ig,ixi)
                vpa_nud = vpdiff(ig,isgn,il) * speed(ie) * vnmult(1) * vnew_D(ik, ie, is)
                gtmp(ixi,ie,ile) = vpa_nud * aj0le(ixi,ie,ile) * gle(ixi,ie,ile)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, is, it, ik, ig, ixi, ie, nud) &
       !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gtmp, vpa_aj0_le, vnmult, &
       !$OMP vnew_D, gle) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ik = ik_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          if(kwork_filter(it,ik)) cycle
          is = is_idx(le_lo,ile)
          do ie = 1, negrid
             nud =  vnmult(1) * vnew_D(ik, ie, is)
             do ixi = 1, nxi_lim
                gtmp(ixi,ie,ile) = vpa_aj0_le(ixi, ie, ile) * nud * gle(ixi, ie, ile)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, s0le, v1y1) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gle(ixi,ie,ile) = gle(ixi,ie,ile) - s0le(ixi,ie,ile) * v1y1(ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    deallocate(v1y1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, is, ig, il, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, ixi_to_il, gtmp, vnmult, vnew_D, vperp_aj1le, gle) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       do ie = 1, negrid
          do ixi=1, nxi_lim
             il = ixi_to_il(ig,ixi)
             ! aj1vp2 = 2 * J1(arg)/arg * vperp^2
             gtmp(ixi,ie,ile) = vnmult(1) * vnew_D(ik,ie,is) * vperp_aj1le(ixi,ie,ile) &
                  * gle(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, w0le, v2y2) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gle(ixi,ie,ile) = gle(ixi,ie,ile) - w0le(ixi,ie,ile) * v2y2(ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    deallocate (gtmp, v2y2)

  end subroutine conserve_lorentz_le_layout

  !> FIXME : Add documentation
  subroutine conserve_diffuse_standard_layout (g, g1)
    use theta_grid, only: ntgrid
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use gs2_layouts, only: g_lo, ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: energy => energy_maxw, al, integrate_moment, negrid
    use dist_fn_arrays, only: aj0, aj1, vpa
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, g1
    complex, dimension (:,:,:), allocatable :: gtmp

    real, dimension (:,:,:), allocatable :: vns

    integer :: isgn, iglo, ik, ie, il, is, it 
    complex, dimension (:,:,:,:), allocatable :: v0y0, v1y1, v2y2    
    logical :: all = .true.

    allocate (v0y0(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v1y1(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (v2y2(-ntgrid:ntgrid, ntheta0, naky, nspec))

    allocate (vns(naky,negrid,nspec))
    allocate (gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))

    vns = vnmult(2)*delvnew

    !This is needed to to ensure the it,ik values we don't set aren't included
    !in the integral (can also be enforced in integrate_moment routine)
    if(any(kwork_filter)) call zero_array(gtmp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, gtmp, vnmult, vnew_E, aj0, g) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v0 = nu_E E J0 f0
          gtmp(:,isgn,iglo) = vnmult(2)*vnew_E(ik,ie,is)*aj0(:,iglo) &
               * g(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, v0y0, all)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, g1, g, ieqzip, v0y0, bz0) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g(:,isgn,iglo) - v0y0(:,it,ik,is) * bz0(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, ie, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, gtmp, vns, vpa, aj0, g1) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v1 = (nus-nud) vpa J0 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*vpa(:,isgn,iglo)*aj0(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, v1y1, all)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, g1, ieqzip, v1y1, bs0) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g1(:,isgn,iglo) = g1(:,isgn,iglo) - v1y1(:,it,ik,is) * bs0(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, ie, il, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, gtmp, vns, energy, al, aj1, g1) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if(kwork_filter(it,ik))cycle
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          ! v2 = (nus-nud) vperp J1 f0
          gtmp(:,isgn,iglo) = vns(ik,ie,is)*energy(ie)*al(il)*aj1(:,iglo) &
               * g1(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (gtmp, v2y2, all)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, g, g1, ieqzip, v2y2, bw0) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g(:,isgn,iglo) = g1(:,isgn,iglo) - v2y2(:,it,ik,is) * bw0(:,isgn,iglo)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate (vns, v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_diffuse_standard_layout

  !> FIXME : Add documentation
  subroutine conserve_diffuse_le_layout (gle)
    use gs2_layouts, only: ik_idx, it_idx, ie_idx, il_idx, is_idx
    use le_grids, only: integrate_moment, negrid, ixi_to_il, nxi, ng2
    use run_parameters, only: ieqzip
    use gs2_layouts, only: le_lo, ig_idx, ik_idx, it_idx
    use redistribute, only: scatter
    use kt_grids, only: kwork_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    complex, dimension (:,:,:), allocatable :: gtmp

    real :: delnu
    integer :: ig, ik, ie, il, is, it
    integer :: ile, ixi
    integer :: nxi_lim
    complex, dimension (:), allocatable :: v0y0, v1y1, v2y2

    ! We only use one of these at a time so could just use a single array.
    ! We keep all three here for clarity and note that these should be small.
    allocate (v0y0(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v1y1(le_lo%llim_proc:le_lo%ulim_alloc))
    allocate (v2y2(le_lo%llim_proc:le_lo%ulim_alloc))

    allocate (gtmp(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))

    !This is needed to to ensure the it,ik values we don't set aren't included
    !in the integral (can also be enforced in integrate_moment routine)
    !if(any(kwork_filter)) gtmp=0.
    call zero_array(gtmp)

    ! The nxi > 2 * ng2 check appears to be checking if we have
    ! trapped particles or not so could be replaced with grid_has_trapped_particles()
    if(le_nxilim_fix .and. nxi > 2*ng2 ) then
        nxi_lim = nxi + 1
    else
        nxi_lim = nxi
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First get v0y0

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, is, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gtmp, vnmult, vnew_E, aj0le, gle) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gtmp(ixi,ie,ile) = vnmult(2) * vnew_E(ik,ie,is) * aj0le(ixi,ie,ile) * gle(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (le_lo, gtmp, v0y0)    ! v0y0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y1 = y0 - v0y0 * z0 / (1 + v0z0)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, v0y0, bz0le) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gle(ixi,ie,ile) = gle(ixi,ie,ile) - v0y0(ile) * bz0le(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v1y1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, is, ig, delnu, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, vnmult, delvnew, nxi_lim, gtmp, &
    !$OMP vpa_aj0_le, gle) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       do ie = 1, negrid
          delnu = vnmult(2) * delvnew(ik, ie, is)
          do ixi = 1, nxi_lim
             gtmp (ixi,ie,ile) = vpa_aj0_le(ixi, ie, ile) * delnu * gle(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (le_lo, gtmp, v1y1)    ! v1y1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get y2 = y1 - v1y1 * s1 / (1 + v1s1)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, bs0le, v1y1) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gle(ixi,ie,ile) = gle(ixi,ie,ile) - bs0le(ixi,ie,ile) * v1y1(ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now get v2y2
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, is, ig, il, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, ixi_to_il, gtmp, &
    !$OMP delvnew, vperp_aj1le, gle, vnmult) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       ik = ik_idx(le_lo,ile)
       it = it_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       is = is_idx(le_lo,ile)
       ig = ig_idx(le_lo,ile)
       do ie=1, negrid
          do ixi = 1, nxi_lim
             il = ixi_to_il(ig,ixi)
             gtmp(ixi,ie,ile) = vnmult(2) * delvnew(ik,ie,is) * vperp_aj1le(ixi,ie,ile) * gle(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (le_lo, gtmp, v2y2)    ! v2y2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finally get x = y2 - v2y2 * w2 / (1 + v2w2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, it, ik, ie, ixi) &
    !$OMP SHARED(le_lo, kwork_filter, negrid, nxi_lim, gle, ieqzip, bw0le, v2y2) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       it = it_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if(kwork_filter(it,ik)) cycle
       if(ieqzip(it,ik)) cycle
       do ie=1, negrid
          do ixi = 1, nxi_lim
             gle(ixi,ie,ile) = gle(ixi,ie,ile) - bw0le(ixi,ie,ile) * v2y2(ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    deallocate (v0y0, v1y1, v2y2, gtmp)

  end subroutine conserve_diffuse_le_layout

  !> FIXME : Add documentation
  subroutine solfp_lorentz_standard_layout (g, gc, gh, diagnostics)
    use theta_grid, only: ntgrid
    use le_grids, only: nlambda, jend, ng2, lambda_map, nxi, il_is_wfb
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: ig_idx, ik_idx, il_idx, is_idx, it_idx, ie_idx
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    use gs2_layouts, only: lz_lo
    use array_utils, only: zero_array

    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g, gc, gh
    integer, optional, intent (in) :: diagnostics
    complex, dimension (:,:), allocatable :: glz, glzc
    complex, dimension (nxi+1) :: delta
    complex :: fac, gwfb
    integer :: ig, ik, il, is, it, je, ie, nxi_scatt, ilz
    logical :: is_wfb

    allocate (glz(nxi+1,lz_lo%llim_proc:lz_lo%ulim_alloc))
    call zero_array(glz)

    if (heating .and. present(diagnostics)) then
       allocate (glzc(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))
       call zero_array(glzc)
    end if

    call gather (lambda_map, g, glz)

    if (heating .and. present(diagnostics)) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ilz, ig, je, il, fac) &
       !$OMP SHARED(lz_lo, jend, ng2, glz, glzc, d1) &
       !$OMP SCHEDULE(static)
       do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
          ig = ig_idx(lz_lo,ilz)

          je = 2*jend(ig)          
          if (je == 0) then
             je = 2*ng2 
          end if

! when il=je-1 below, and we have trapped particles, glz is evaluated at glz(2*jend(ig),ilz).
! this seems like a bug, since there are only 2*jend(ig)-1 grid points and
! the value glz(2*jend(ig),ilz) corresponds to the value of g at xi = 0...this
! doesn't make any sense...MAB

          do il = 1, je-1
             fac = glz(il+1,ilz)-glz(il,ilz)
             glzc(il,ilz) = conjg(fac)*fac*d1(il,ilz)  ! d1 accounts for hC(h) entropy
          end do
       end do
       !$OMP END PARALLEL DO

       call scatter (lambda_map, glzc, gc)

       if (hyper_colls) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ilz, ig, je, il, fac) &
          !$OMP SHARED(lz_lo, jend, ng2, glz, glzc, h1) &
          !$OMP SCHEDULE(static)
          do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
             ig = ig_idx(lz_lo,ilz)
             
             je = 2*jend(ig)          
             if (je == 0) then
                je = 2*ng2 
             end if
             
             do il = 1, je-1
                fac = glz(il+1,ilz)-glz(il,ilz)
                glzc(il,ilz) = conjg(fac)*fac*h1(il,ilz)  ! h1 accounts for hH(h) entropy
             end do
          end do
          !$OMP END PARALLEL DO

          call scatter (lambda_map, glzc, gh)
       end if
    end if

    ! solve for glz row by row
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ilz, ik, it, is, ie, ig, je, nxi_scatt, &
    !$OMP il, gwfb, delta, is_wfb) &
    !$OMP SHARED(lz_lo, kwork_filter, ieqzip, vnew, force_collisions, jend, ng2, vpar_zero_mean, &
    !$OMP glz, special_wfb_lorentz, ql, c1, betaa) &
    !$OMP SCHEDULE(static)
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       is = is_idx(lz_lo,ilz)
       ik = ik_idx(lz_lo,ilz)
       if ( (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) .and. .not. force_collisions) cycle
       it = it_idx(lz_lo,ilz)
       if(kwork_filter(it,ik))cycle
       if (ieqzip(it,ik)) cycle
       ie = ie_idx(lz_lo,ilz)
       ig = ig_idx(lz_lo,ilz)

!CMRDDGC, 10/2/2014: 
! Fixes for wfb treatment below, use same je definition in ALL cases
!   je  = #physical xi values + 1 
!       NB +1 above WITH TRAPPED is duplicate xi=vpar=0 point with isign=2
!          +1 above WITHOUT TRAPPED is entirely unphysical extra grid point
!  je-1 = #physical xi values removing unphysical/duplicate extra point
       je=max(2*jend(ig),2*ng2+1)
       nxi_scatt=je-1
       is_wfb = il_is_wfb(jend(ig))
       
       if (je==2*jend(ig) .and. vpar_zero_mean) then 
          ! MRH if we have trapped particles and hence 2 vpar=0 points on the grid
          ! use the average of both vpar = 0 points in the lorentz diffusion in pitch angle
          ! this turns out to be necessary to suppress numerical instability
          ! when we handle pitch angle scattering physically at theta = +/- pi
          ! by setting special_wfb_lorentz =.false.
          glz(jend(ig),ilz)=(glz(jend(ig),ilz)+glz(2*jend(ig),ilz))/2.0
       endif 
       
       glz(je,ilz)=0.0  ! zero unphysical/duplicate extra point
       if (is_wfb .and. special_wfb_lorentz) then
!CMRDDGC:  special_wfb_lorentz = t  => unphysical handling of wfb at bounce pt: 
!          remove wfb from collisions, reinsert later
!
! first save gwfb for reinsertion later
          gwfb = glz(ng2+1,ilz) 
! then remove vpa = 0 point, weight 0: (CMR confused by this comment!)  
          glz(ng2+1:je-2,ilz) = glz(ng2+2:je-1,ilz)
          nxi_scatt=nxi_scatt-1
       end if
!CMRDDGCend

       ! right and left sweeps for tridiagonal solve:

       delta(1) = glz(1,ilz)
       do il = 1, nxi_scatt
          delta(il+1) = glz(il+1,ilz) - ql(il+1,ilz)*delta(il)
       end do
       
!CMR, 5/9/14:
! Correction to next line fixes a bug introduced at r2699 which
! affected Lorentz scattering operator at wfb bounce point (theta=pi)
! with default setting of special_wfb_lorentz=t        
       glz(nxi_scatt+1,ilz) = delta(nxi_scatt+1)*betaa(nxi_scatt+1,ilz)
       do il = nxi_scatt, 1, -1
          glz(il,ilz) = (delta(il) - c1(il,ilz)*glz(il+1,ilz))*betaa(il,ilz)
       end do

!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
       ! interpolation described above mysteriously causing numerical instability
       ! stabilized by using old (pre-collision) value of g for wfb
       if (is_wfb .and. special_wfb_lorentz) then
          glz(ng2+2:je-1,ilz) = glz(ng2+1:je-2,ilz)
          glz(ng2+1,ilz) = gwfb
       end if

!
! bug fixed 4.14.03
! was only a problem with collisions but with no trapped particles
! because of an array index out of bounds 
! Overall, this was a rare bug.
!

       ! update the duplicate vpar=0 point with the post pitch angle scattering value 
       if (jend(ig) /= 0) glz(2*jend(ig),ilz) = glz(jend(ig),ilz)

    end do
    !$OMP END PARALLEL DO

    call scatter (lambda_map, glz, g)

    deallocate (glz)
    if (allocated(glzc)) deallocate (glzc)

  end subroutine solfp_lorentz_standard_layout

  !> FIXME : Add documentation
  subroutine solfp_lorentz_le_layout (gle, diagnostics)
    use le_grids, only: jend, ng2, negrid, nxi, integrate_moment, il_is_wfb
    use gs2_layouts, only: ig_idx, ik_idx, is_idx, it_idx
    use run_parameters, only: ieqzip
    use gs2_layouts, only: le_lo
    use dist_fn_arrays, only: c_rate
    use kt_grids, only: kwork_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle
    integer, optional, intent (in) :: diagnostics
    complex, dimension (:), allocatable :: tmp
    complex, dimension (:,:,:), allocatable :: glec
    complex, dimension (nxi+1) :: delta
    complex :: fac, gwfb
    integer :: ig, ik, il, is, je, it, ie, nxi_scatt, ile
    logical :: is_wfb

    if (heating .and. present(diagnostics)) then
       allocate (tmp(le_lo%llim_proc:le_lo%ulim_alloc)) ; call zero_array(tmp)
       allocate (glec(nxi+1, negrid+1, le_lo%llim_proc:le_lo%ulim_alloc))
       call zero_array(glec)
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, ig, je, ie, il, fac) &
       !$OMP SHARED(le_lo, jend, ng2, negrid, gle, glec, d1le) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ig = ig_idx(le_lo,ile)

          je = 2*jend(ig)
          if (je == 0) then
             je = 2*ng2 
          end if

! when il=je-1 below, and we have trapped particles, gle is evaluated at gle(2*jend(ig),ie,ile).
! this seems like a bug, since there are only 2*jend(ig)-1 grid points and
! the value gle(2*jend(ig),ie,ile) corresponds to the value of g at xi = 0...this
! doesn't make any sense...MAB

          do ie = 1, negrid
             do il = 1, je-1
                fac = gle(il+1,ie,ile)-gle(il,ie,ile)
                glec(il,ie,ile) = conjg(fac)*fac*d1le(il,ie,ile)  ! d1le accounts for hC(h) entropy
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       call integrate_moment (le_lo, glec, tmp)
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ile, ig, it, ik, is) &
       !$OMP SHARED(le_lo, c_rate, tmp) &
       !$OMP SCHEDULE(static)
       do ile = le_lo%llim_proc, le_lo%ulim_proc
          ig = ig_idx(le_lo,ile)
          it = it_idx(le_lo,ile)
          ik = ik_idx(le_lo,ile)
          is = is_idx(le_lo,ile)
          c_rate(ig,it,ik,is,1) = tmp(ile)
       end do
       !$OMP END PARALLEL DO

       if (hyper_colls) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ile, ig, je, ie, il, fac) &
          !$OMP SHARED(le_lo, jend, ng2, negrid, gle, glec, h1le) &
          !$OMP SCHEDULE(static)
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             ig = ig_idx(le_lo,ile)

             je = 2*jend(ig)
             if (je == 0) then
                je = 2*ng2
             end if

             do ie = 1, negrid
                do il = 1, je-1
                   fac = gle(il+1,ie,ile)-gle(il,ie,ile)
                   glec(il,ie,ile) = conjg(fac)*fac*h1le(il,ie,ile)  ! h1le accounts for hH(h) entropy
                end do
             end do
          end do
          !$OMP END PARALLEL DO

          call integrate_moment (le_lo, glec, tmp)

          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ile, ig, it, ik, is) &
          !$OMP SHARED(le_lo, c_rate, tmp) &
          !$OMP SCHEDULE(static)
          do ile = le_lo%llim_proc, le_lo%ulim_proc
             ig = ig_idx(le_lo,ile)
             it = it_idx(le_lo,ile)
             ik = ik_idx(le_lo,ile)
             is = is_idx(le_lo,ile)
             c_rate(ig,it,ik,is,2) = tmp(ile)
          end do
          !$OMP END PARALLEL DO
       end if
       deallocate(tmp, glec)
    end if

    ! solve for gle row by row
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, is, it, ik, ig, il, ie, je, nxi_scatt, &
    !$OMP gwfb, delta, is_wfb) &
    !$OMP SHARED(le_lo, kwork_filter, ieqzip, vnew, force_collisions, jend, ng2, special_wfb_lorentz, &
    !$OMP vpar_zero_mean, negrid, gle, qle, betaale, c1le) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if ( (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) .and. .not. force_collisions) cycle
       it = it_idx(le_lo,ile)
       if (kwork_filter(it,ik)) cycle
       if (ieqzip(it,ik)) cycle
       ig = ig_idx(le_lo,ile)

!CMRDDGC, 10/2/1014: 
! Fixes for wfb treatment below, use same je definition in ALL cases
!   je  = #physical xi values at location, includes duplicate point at vpar=0
!  je-1 = #physical xi values removing duplicate vpar=0 point
       je=max(2*jend(ig),2*ng2+1)
       nxi_scatt=je-1
       is_wfb = il_is_wfb(jend(ig))
       if (is_wfb .and. special_wfb_lorentz) nxi_scatt = nxi_scatt - 1

! These fixes address previous comment by CMR: 
! "surely wfb's bounce point should be handled like any other trapped bp?"
!CMR, 1/11/2013:
! Numerical instability referred to below more likely arises from wfb failing
! to satisfy trapping condition after invert_rhs.
!
       do ie = 1, negrid
          if (je==2*jend(ig) .and. vpar_zero_mean) then
             ! MRH if we have trapped particles and hence 2 vpar=0 points on the grid
             ! use the average of both vpar = 0 points in the lorentz diffusion in pitch angle
             ! this turns out to be necessary to suppress numerical instability
             ! when we handle pitch angle scattering physically at theta = +/- pi
             ! by setting special_wfb_lorentz =.false.
             gle(jend(ig),ie,ile)=(gle(jend(ig),ie,ile)+gle(2*jend(ig),ie,ile))/2.0
          endif 
          gle(je,ie,ile)=0.0d0  ! zero redundant duplicate xi, isign=2 for vpar=0!
          if (is_wfb .and. special_wfb_lorentz) then
!CMRDDGC:  special_wfb_lorentz = t  => unphysical handling of wfb at bounce pt: 
!          remove wfb from collisions, reinsert later
!
! first save gwfb for reinsertion later
             gwfb = gle(ng2+1,ie,ile)
! then remove vpa = 0 point, weight 0: (CMR confused by this comment!)  
             gle(ng2+1:je-2,ie,ile) = gle(ng2+2:je-1,ie,ile)
          end if
!CMRDDGCend

          ! right and left sweeps for tridiagonal solve:
          
          delta(1) = gle(1,ie,ile)
          do il = 1, nxi_scatt
             delta(il+1) = gle(il+1,ie,ile) - qle(il+1,ie,ile)*delta(il)
          end do
!CMR, 5/9/14:
! Correction to next line fixes a bug introduced at r2699 which
! affected Lorentz scattering operator at wfb bounce point (theta=pi)
! with default setting of special_wfb_lorentz=t        
          gle(nxi_scatt+1,ie,ile) = delta(nxi_scatt+1)*betaale(nxi_scatt+1,ie,ile)
          do il = nxi_scatt, 1, -1
             gle(il,ie,ile) = (delta(il) - c1le(il,ie,ile)*gle(il+1,ie,ile))*betaale(il,ie,ile)
          end do
!       ! interpolate to obtain glz(vpa = 0) point for wfb
!       ! and insert this point into glz
          ! interpolation described above mysteriously causing numerical instability
          ! stabilized by using old (pre-collision) value of g for wfb
          if (is_wfb .and. special_wfb_lorentz) then
             gle(ng2+2:je-1,ie,ile)=gle(ng2+1:je-2,ie,ile)
             gle(ng2+1,ie,ile)=gwfb
!          gle(ng2+1,ie,ile) = 0.5*(gle(ng2,ie,ile)+gle(ng2+1,ie,ile))
          end if

!CMR, 1/11/2013:
! next line ensures bounce condition is satisfied after lorentz collisions
! this is right thing to do, but it would mask any prior bug in trapping condition!
          ! update the duplicate vpar=0 point with the post pitch angle scattering value 
          if (jend(ig) /= 0) gle(2*jend(ig),ie,ile) = gle(jend(ig),ie,ile)

       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine solfp_lorentz_le_layout

  !> Energy diffusion subroutine used with energy layout (not le_layout)
  !> this is always the case when initializing the conserving terms,
  !> otherwise is the case if use_le_layout is no specified in the input file.
  subroutine solfp_ediffuse_standard_layout (g)
    use theta_grid, only: ntgrid
    use le_grids, only: negrid, forbid
    use le_grids, only: energy_map
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, il_idx, is_idx, e_lo, g_lo
    use redistribute, only: gather, scatter
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g

    integer :: ie, is, ig, il, it, ik
    complex, dimension (negrid) :: delta
    complex, dimension (:,:), allocatable :: ged
    integer :: ielo

    allocate (ged(negrid+1,e_lo%llim_proc:e_lo%ulim_alloc)) ; call zero_array(ged)

    call gather (energy_map, g, ged)

    ! solve for ged row by row
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ielo, it, ik, is, ig, il, delta, ie) &
    !$OMP SHARED(e_lo, kwork_filter, ieqzip, vnew, force_collisions, forbid, negrid, ged, eql, ec1, ebetaa) &
    !$OMP SCHEDULE(static)
    do ielo = e_lo%llim_proc, e_lo%ulim_proc
       is = is_idx(e_lo,ielo)
       ik = ik_idx(e_lo,ielo)
       if ( (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) .and. .not. force_collisions) cycle
       it = it_idx(e_lo,ielo)
       if(kwork_filter(it,ik))cycle
       if (ieqzip(it,ik)) cycle
       ig = ig_idx(e_lo,ielo)
       il = il_idx(e_lo,ielo)
       if (forbid(ig,il)) cycle

       delta(1) = ged(1,ielo)
       do ie = 1, negrid-1
          delta(ie+1) = ged(ie+1,ielo) - eql(ie+1,ielo)*delta(ie)
       end do
       
       ged(negrid+1,ielo) = 0.0
       do ie = negrid, 1, -1
          ged(ie,ielo) = (delta(ie) - ec1(ie,ielo)*ged(ie+1,ielo))*ebetaa(ie,ielo)
       end do

    end do
    !$OMP END PARALLEL DO

    call scatter (energy_map, ged, g)

    deallocate (ged)

  end subroutine solfp_ediffuse_standard_layout

  !> FIXME : Add documentation
  subroutine solfp_ediffuse_le_layout (gle)
    use le_grids, only: nxi, negrid, forbid, ixi_to_il, ng2
    use gs2_layouts, only: ig_idx, it_idx, ik_idx, is_idx, le_lo
    use run_parameters, only: ieqzip
    use kt_grids, only: kwork_filter
    implicit none

    complex, dimension (:,:,le_lo%llim_proc:), intent (in out) :: gle

    integer :: ie, is, ig, il
    complex, dimension (negrid) :: delta
    integer :: ile, ixi, ik, it
    integer :: nxi_lim
    
    ! The nxi > 2 * ng2 check appears to be checking if we have
    ! trapped particles or not so could be replaced with grid_has_trapped_particles()
    if(le_nxilim_fix .and. nxi > 2*ng2 ) then
        nxi_lim = nxi + 1
    else
        nxi_lim = nxi
    endif
    
    ! solve for gle row by row
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ile, is, it, ik, ig, ixi, il, ie, delta) &
    !$OMP SHARED(le_lo, vnew, force_collisions, kwork_filter, ieqzip, ixi_to_il, forbid, gle, negrid, &
    !$OMP eqle, ec1le, ebetaale, nxi_lim) &
    !$OMP SCHEDULE(static)
    do ile = le_lo%llim_proc, le_lo%ulim_proc
       is = is_idx(le_lo,ile)
       ik = ik_idx(le_lo,ile)
       if ( (abs(vnew(ik,1,is)) < 2.0*epsilon(0.0)) .and. .not. force_collisions) cycle
       it=it_idx(le_lo,ile)
       if (kwork_filter(it,ik)) cycle
       if (ieqzip(it,ik)) cycle
       ig = ig_idx(le_lo,ile)

       do ixi = 1, nxi_lim
          il = ixi_to_il(ig,ixi)
          if (forbid(ig,il)) cycle

          delta(1) = gle(ixi,1,ile)
          do ie = 1, negrid-1
             delta(ie+1) = gle(ixi,ie+1,ile) - eqle(ixi,ie+1,ile)*delta(ie)
          end do

          gle(ixi,negrid+1,ile) = 0.0
          do ie = negrid, 1, -1
             gle(ixi,ie,ile) = (delta(ie) - ec1le(ixi,ie,ile)*gle(ixi,ie+1,ile)) &
                  * ebetaale(ixi,ie,ile)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine solfp_ediffuse_le_layout

  !> FIXME : Add documentation
  subroutine init_vpdiff
    use le_grids, only: al, nlambda, jend, ng2, il_is_wfb, il_is_passing
    use theta_grid, only: ntgrid, bmag
    use array_utils, only: zero_array
    implicit none

    integer :: il, ig, je, te
    real :: slb0, slb1, slb2, slbl, slbr

    if (.not. allocated(vpdiff) .and. conservative) then

       allocate (vpdiff(-ntgrid:ntgrid,2,nlambda)) ; call zero_array(vpdiff)
       
       do ig = -ntgrid, ntgrid
          je = jend(ig)
          if (il_is_passing(je) .or. il_is_wfb(je)) then
             te = ng2
          else
             te = je
          end if
          do il = 2, te-1
             slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
             slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
             slb2 = safe_sqrt(1.0 - bmag(ig)*al(il+1))
             
             slbl = (slb1 + slb0)*0.5  ! xi(j-1/2)
             slbr = (slb1 + slb2)*0.5  ! xi(j+1/2)
             
             vpdiff(ig,1,il) = (slbl**2 - slbr**2)/pitch_weights(ig,il)
          end do

          ! boundary at xi = 1
          slb1 = safe_sqrt(1.0 - bmag(ig)*al(1))
          slb2 = safe_sqrt(1.0 - bmag(ig)*al(2))
          slbr = 0.5*(slb1 + slb2)
          vpdiff(ig,1,1) = (1.0 - slbr**2)/pitch_weights(ig,1)
          
          ! boundary at xi = 0
          il = te
          slb0 = safe_sqrt(1.0 - bmag(ig)*al(il-1))
          if (te == ng2) then
             slb1 = safe_sqrt(1.0 - bmag(ig)*al(il))
             slb2 = -slb1
          else
             slb1 = 0.0
             slb2 = -slb0
          end if
          
          slbl = (slb1 + slb0)*0.5
          slbr = (slb1 + slb2)*0.5
          vpdiff(ig,1,il) = (slbl**2 - slbr**2)/pitch_weights(ig,il)
          
          vpdiff(ig,2,:) = -vpdiff(ig,1,:)
          
       end do
       
    end if
    
  end subroutine init_vpdiff

  !> Forces recalculation of coefficients in collision operator
  !! when timestep changes. Currently just calls [[finish_collisions]]
  subroutine reset_init
    call finish_collisions
  end subroutine reset_init

  !> Forces recalculation of coefficients in collision operator
  !! when timestep changes.
  subroutine finish_collisions
    use dist_fn_arrays, only: c_rate
    implicit none

    ! This is a bit of a hack to make sure that we get the correct vnmult
    ! value during a timestep change. Whilst we do restore vnmult from the
    ! restart file, this unfortunately doesn't happen until after we've initialised
    ! the collision operator so we choose to store this in vnm_init, which is
    ! where we get our initial value of vnmult from anyway.
    vnm_init = vnmult
    vnmult = -1.0

    accelerated_x = .false. ; accelerated_v = .false.
    ediffinit = .false. ; lzinit = .false. ; leinit = .false.
    initialized = .false.  

    if (allocated(c_rate)) deallocate (c_rate)
    if (allocated(z0)) deallocate (z0, w0, s0)
    if (allocated(bz0)) deallocate (bz0, bw0, bs0)
    if (allocated(vnew)) deallocate (vnew, vnew_s, vnew_D, vnew_E, delvnew)
    if (allocated(vnewh)) deallocate (vnewh)
    if (allocated(pitch_weights)) deallocate(pitch_weights)

    if(use_le_layout) then
      if (allocated(c1le)) then
         deallocate (c1le, betaale, qle)
         if (heating) deallocate (d1le, h1le)
      end if
      if (allocated(ec1le)) deallocate (ec1le, ebetaale, eqle)
      if (allocated(s0le)) deallocate(s0le)
      if (allocated(w0le)) deallocate(w0le)
      if (allocated(z0le)) deallocate(z0le)
      if (allocated(aj0le)) deallocate(aj0le)
      if (allocated(vperp_aj1le)) deallocate(vperp_aj1le)
      if (allocated(vpa_aj0_le)) deallocate(vpa_aj0_le)
      if (allocated(bs0le)) deallocate(bs0le)
      if (allocated(bw0le)) deallocate(bw0le)
      if (allocated(bz0le)) deallocate(bz0le)
    else
      if (allocated(c1)) then
         deallocate (c1, betaa, ql)
         if (heating) deallocate (d1, h1)
      end if
      if (allocated(ec1)) deallocate (ec1, ebetaa, eql)
    end if
    if (allocated(vpdiff)) deallocate (vpdiff)
    if (allocated(dtot)) deallocate (dtot, fdf, fdb)

    call collisions_config%reset()
  end subroutine finish_collisions

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_collisions_config(collisions_config_in)
    use mp, only: mp_abort
    type(collisions_config_type), intent(in), optional :: collisions_config_in
    if (initialized) then
       call mp_abort("Trying to set collisions_config when already initialized.", to_screen = .true.)
    end if
    if (present(collisions_config_in)) then
       collisions_config = collisions_config_in
    end if
  end subroutine set_collisions_config

  !> Get the module level config instance
  function get_collisions_config()
    type(collisions_config_type) :: get_collisions_config
    get_collisions_config = collisions_config
  end function get_collisions_config

  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the collisions_knobs namelist and populates the member variables
  subroutine read_collisions_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(collisions_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: collision_model, ediff_scheme, lorentz_scheme
    integer :: ncheck, timesteps_between_collisions
    logical :: adjust, conservative, conserve_forbid_zero, conserve_moments, const_v, ei_coll_only, force_collisions, heating, hypermult, le_nxilim_fix, resistivity, special_wfb_lorentz, split_collisions, test, use_le_layout
    logical :: vary_vnew, vpar_zero_mean
    real :: cfac, etol, etola, ewindow, ewindowa, vnfac, vnslow
    
    namelist /collisions_knobs/ adjust, cfac, collision_model, conservative, conserve_forbid_zero, conserve_moments, const_v, ediff_scheme, ei_coll_only, etol, etola, ewindow, ewindowa, force_collisions, heating, hypermult, le_nxilim_fix, &
         lorentz_scheme, ncheck, resistivity, special_wfb_lorentz, split_collisions, test, timesteps_between_collisions, use_le_layout, vary_vnew, vnfac, vnslow, vpar_zero_mean

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    adjust = self%adjust
    cfac = self%cfac
    collision_model = self%collision_model
    conservative = self%conservative
    conserve_forbid_zero = self%conserve_forbid_zero
    conserve_moments = self%conserve_moments
    const_v = self%const_v
    ediff_scheme = self%ediff_scheme
    ei_coll_only = self%ei_coll_only
    etol = self%etol
    etola = self%etola
    ewindow = self%ewindow
    ewindowa = self%ewindowa
    force_collisions = self%force_collisions
    heating = self%heating
    hypermult = self%hypermult
    le_nxilim_fix = self%le_nxilim_fix
    lorentz_scheme = self%lorentz_scheme
    ncheck = self%ncheck
    resistivity = self%resistivity
    special_wfb_lorentz = self%special_wfb_lorentz
    split_collisions = self%split_collisions
    test = self%test
    timesteps_between_collisions = self%timesteps_between_collisions
    use_le_layout = self%use_le_layout
    vary_vnew = self%vary_vnew
    vnfac = self%vnfac
    vnslow = self%vnslow
    vpar_zero_mean = self%vpar_zero_mean

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = collisions_knobs)

    ! Now copy from local variables into type members
    self%adjust = adjust
    self%cfac = cfac
    self%collision_model = collision_model
    self%conservative = conservative
    self%conserve_forbid_zero = conserve_forbid_zero
    self%conserve_moments = conserve_moments
    self%const_v = const_v
    self%ediff_scheme = ediff_scheme
    self%ei_coll_only = ei_coll_only
    self%etol = etol
    self%etola = etola
    self%ewindow = ewindow
    self%ewindowa = ewindowa
    self%force_collisions = force_collisions
    self%heating = heating
    self%hypermult = hypermult
    self%le_nxilim_fix = le_nxilim_fix
    self%lorentz_scheme = lorentz_scheme
    self%ncheck = ncheck
    self%resistivity = resistivity
    self%special_wfb_lorentz = special_wfb_lorentz
    self%split_collisions = split_collisions
    self%test = test
    self%timesteps_between_collisions = timesteps_between_collisions
    self%use_le_layout = use_le_layout
    self%vary_vnew = vary_vnew
    self%vnfac = vnfac
    self%vnslow = vnslow
    self%vpar_zero_mean = vpar_zero_mean

    self%exist = exist
  end subroutine read_collisions_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_collisions_config(self, unit)
    implicit none
    class(collisions_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("adjust", self%adjust, unit_internal)
    call self%write_key_val("cfac", self%cfac, unit_internal)
    call self%write_key_val("collision_model", self%collision_model, unit_internal)
    call self%write_key_val("conservative", self%conservative, unit_internal)
    call self%write_key_val("conserve_forbid_zero", self%conserve_forbid_zero, unit_internal)
    call self%write_key_val("conserve_moments", self%conserve_moments, unit_internal)
    call self%write_key_val("const_v", self%const_v, unit_internal)
    call self%write_key_val("ediff_scheme", self%ediff_scheme, unit_internal)
    call self%write_key_val("ei_coll_only", self%ei_coll_only, unit_internal)
    call self%write_key_val("etol", self%etol, unit_internal)
    call self%write_key_val("etola", self%etola, unit_internal)
    call self%write_key_val("ewindow", self%ewindow, unit_internal)
    call self%write_key_val("ewindowa", self%ewindowa, unit_internal)
    call self%write_key_val("force_collisions", self%force_collisions, unit_internal)
    call self%write_key_val("heating", self%heating, unit_internal)
    call self%write_key_val("hypermult", self%hypermult, unit_internal)
    call self%write_key_val("le_nxilim_fix", self%le_nxilim_fix, unit_internal)
    call self%write_key_val("lorentz_scheme", self%lorentz_scheme, unit_internal)
    call self%write_key_val("ncheck", self%ncheck, unit_internal)
    call self%write_key_val("resistivity", self%resistivity, unit_internal)
    call self%write_key_val("special_wfb_lorentz", self%special_wfb_lorentz, unit_internal)
    call self%write_key_val("split_collisions", self%split_collisions, unit_internal)
    call self%write_key_val("test", self%test, unit_internal)
    call self%write_key_val("timesteps_between_collisions", self%timesteps_between_collisions, unit_internal)
    call self%write_key_val("use_le_layout", self%use_le_layout, unit_internal)
    call self%write_key_val("vary_vnew", self%vary_vnew, unit_internal)
    call self%write_key_val("vnfac", self%vnfac, unit_internal)
    call self%write_key_val("vnslow", self%vnslow, unit_internal)
    call self%write_key_val("vpar_zero_mean", self%vpar_zero_mean, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_collisions_config

  !> Resets the config object to the initial empty state
  subroutine reset_collisions_config(self)
    class(collisions_config_type), intent(in out) :: self
    type(collisions_config_type) :: empty
    select type (self)
    type is (collisions_config_type)
       self = empty
    end select
  end subroutine reset_collisions_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_collisions_config(self)
    use mp, only: broadcast
    implicit none
    class(collisions_config_type), intent(in out) :: self
    call broadcast(self%adjust)
    call broadcast(self%cfac)
    call broadcast(self%collision_model)
    call broadcast(self%conservative)
    call broadcast(self%conserve_forbid_zero)
    call broadcast(self%conserve_moments)
    call broadcast(self%const_v)
    call broadcast(self%ediff_scheme)
    call broadcast(self%ei_coll_only)
    call broadcast(self%etol)
    call broadcast(self%etola)
    call broadcast(self%ewindow)
    call broadcast(self%ewindowa)
    call broadcast(self%force_collisions)
    call broadcast(self%heating)
    call broadcast(self%hypermult)
    call broadcast(self%le_nxilim_fix)
    call broadcast(self%lorentz_scheme)
    call broadcast(self%ncheck)
    call broadcast(self%resistivity)
    call broadcast(self%special_wfb_lorentz)
    call broadcast(self%split_collisions)
    call broadcast(self%test)
    call broadcast(self%timesteps_between_collisions)
    call broadcast(self%use_le_layout)
    call broadcast(self%vary_vnew)
    call broadcast(self%vnfac)
    call broadcast(self%vnslow)
    call broadcast(self%vpar_zero_mean)

    call broadcast(self%exist)
  end subroutine broadcast_collisions_config
  
  !> Gets the default name for this namelist
  function get_default_name_collisions_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_collisions_config
    get_default_name_collisions_config = "collisions_knobs"
  end function get_default_name_collisions_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_collisions_config()
    implicit none
    logical :: get_default_requires_index_collisions_config
    get_default_requires_index_collisions_config = .false.
  end function get_default_requires_index_collisions_config
end module collisions

