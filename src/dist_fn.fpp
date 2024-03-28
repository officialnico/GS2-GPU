! Notes from BD, 7/2011:
!
! Need to extend the verr tools to include delta B_parallel integrals
! There are new factors of 1/B here and there which I do not understand.  

!> The principal function of this module is to evolve the distribution
!> function, that is, to advance the discrete gyrokinetic equation.
!> This involves calculating the source and dealing with the
!> complexities of the parallel boundary conditions.  In addition it
!> contains a routine for implementing perpendicular velocity shear
!> and calculating the right-hand side of the field equation, as well
!> as a host of other functions.
module dist_fn
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use redistribute, only: redist_type
  implicit none

  private

  public :: init_dist_fn, finish_dist_fn
  
  !> Initializes a limited selection of arrays,
  !> for example g, gnew, vperp2, typically those which
  !> are needed by other modules that don't need
  !> dist_fn to be fully initialized (e.g. nonlinear_terms)
  !> This initialization level depends on grid sizes.
  public :: init_dist_fn_arrays

  !> Deallocates arrays allocated in init_dist_fn_arrays
  public :: finish_dist_fn_arrays

  !> Reads the dist_fn_knobs, source knobs, and
  !> dist_fn_species_knobs namelists. This has be done independently
  !> of initializing the distribution function because
  !> it needs to be possible to override parameters such
  !> as g_exb.
  public :: init_dist_fn_parameters
  public :: finish_dist_fn_parameters

  !> Initializes parallel boundary conditions. This level
  !> depends on geometry.
  public :: init_dist_fn_level_1
  public :: finish_dist_fn_level_1

  !> Initializes bessel functions and field_eq. Note 
  !> that this level depends on species paramters.
  public :: init_dist_fn_level_2
  public :: finish_dist_fn_level_2

  !> Fully initialize the dist_fn module. Note that
  !> this level depends on the size of the timestep.
  public :: init_dist_fn_level_3
  public :: finish_dist_fn_level_3

  public :: set_overrides

  interface set_overrides
    module procedure set_profiles_overrides
    module procedure set_optimisations_overrides
  end interface set_overrides

  public :: read_parameters, wnml_dist_fn, wnml_dist_fn_species, check_dist_fn
  public :: timeadv, exb_shear, g_exb, g_exbfac
  public :: collisions_advance
  public :: getfieldeq, getan, getmoms, getemoms, getmoms_gryfx_dist
  public :: getfieldeq_nogath
  public :: flux_dist
  public :: get_heat
  public :: omega0, gamma0
  public :: write_f
  public :: M_class, N_class, i_class
  public :: l_links, r_links, itright, itleft, boundary
  public :: supercell_labels, n_supercells, n_cells
  public :: enforce_linked_boundary_conditions
  public :: get_verr, get_gtran, write_fyx, collision_error
  public :: get_init_field, getan_nogath
  public :: flux_vs_theta_vs_vpa
  public :: pflux_vs_theta_vs_vpa
  public :: wdrift_func
  
  public :: gamtot,gamtot1,gamtot2
  public :: getmoms_notgc, init_mom_coeff

!CMR, 25/1/13:
!  add public variables below for init_g
  public :: boundary_option_switch, boundary_option_linked
  public :: pass_right, init_pass_ends
  public :: get_leftmost_it
  public :: gridfac1, fl_avg, apfac
  public :: adiabatic_option_switch, adiabatic_option_fieldlineavg

  public :: def_parity, even

  public :: gf_lo_integrate
  public :: calculate_flux_surface_average

  public :: calculate_potentials_from_nonadiabatic_dfn
  public :: flux_tube_field_like_to_ballooning_space
  public :: check_linked_boundaries_are_satisfied

  public :: dist_fn_config_type
  public :: dist_fn_species_config_type
  public :: source_config_type
  public :: set_dist_fn_config
  public :: get_dist_fn_config
  public :: get_dist_fn_species_config
  public :: get_source_config

  ! The following routines are sometimes used for debugging or testing purposes
  public :: dump_homogeneous_solution, dump_current_source_term, get_field_inconsistency
  public :: check_getan, get_mom_glm

  ! knobs
  complex, dimension (:), allocatable :: fexp ! (nspec)
  real, dimension (:), allocatable :: bkdiff  ! (nspec)
  real :: gridfac, apfac, driftknob, tpdriftknob, poisfac, vparknob
  real :: t0, omega0, gamma0, source0
  real :: phi_ext, afilter
  real :: wfb, g_exb, g_exbfac, omprimfac, btor_slab, mach
  real :: g_exb_start_time, g_exb_error_limit
  integer :: g_exb_start_timestep
  logical :: hyper_in_initialisation
  logical :: start_from_previous_solution
  logical :: dfexist, skexist, nonad_zero, lf_default, lf_decompose, esv, opt_source

  integer :: adiabatic_option_switch
  integer, parameter :: adiabatic_option_default = 1, &
       adiabatic_option_zero = 2, &
       adiabatic_option_fieldlineavg = 3, &
       adiabatic_option_yavg = 4

  integer :: source_option_switch
  integer, parameter :: source_option_full = 1, &
       source_option_phiext_full = 5, source_option_homogeneous = 6
  
  integer :: boundary_option_switch
  integer, parameter :: boundary_option_default = 0, &
       boundary_option_zero = 1, &
       boundary_option_self_periodic = 2, &
       boundary_option_alternate_zero = 3, &
       boundary_option_linked = 4

  logical :: def_parity, even
  logical :: zero_forbid
  logical :: gf_lo_integrate
  logical :: mult_imp
!Not sure why these are module level as only used in one routine
  logical :: increase = .true., decrease = .false. 

#ifdef LOWFLOW
  complex, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, 2, -g-layout-)
#else
  real, dimension (:,:,:), allocatable :: wdrift
  ! (-ntgrid:ntgrid, 2, -g-layout-)
#endif

  real, dimension (:,:,:), allocatable :: wdriftttp
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  real, dimension (:,:,:), allocatable :: wstar
  ! (naky,negrid,nspec) replicated

  ! fieldeq
  real, dimension (:,:,:), allocatable :: gamtot, gamtot1, gamtot2, gamtot3
  ! (-ntgrid:ntgrid,ntheta0,naky) replicated

  complex, dimension (:,:,:), allocatable :: a, b, r, ainv
  ! (-ntgrid:ntgrid, 2, -g-layout-)

  real, dimension (:,:,:), allocatable :: gridfac1
  ! (-ntgrid:ntgrid,ntheta0,naky)

  real, dimension(:, :, :), allocatable :: vpar

#ifndef SHMEM
  complex, dimension (:,:,:), allocatable :: g_h
#else
  complex, dimension (:,:,:), allocatable, target :: g_h
#endif
  ! (-ntgrid:ntgrid,2, -g-layout-)

  complex, dimension (:,:,:), allocatable :: g_adj
  ! (N(links), 2, -g-layout-)

  ! exb shear
  integer, dimension(:), allocatable :: jump, ikx_indexed

  ! set_source
  real, dimension(:,:), allocatable :: ufac

  ! set_source_opt
  complex, dimension (:,:,:,:), allocatable :: source_coeffs

  ! getfieldeq1
  real, allocatable, dimension(:,:) :: awgt
  complex, allocatable, dimension(:,:) :: fl_avg !Changed

  ! get_verr
  real, dimension (:,:), allocatable :: kmax

  ! connected bc
  integer, dimension (:,:), allocatable :: itleft, itright
  ! (naky,ntheta0)

  !> Labels the supercells as a function of {it, ik}. Primarily
  !> to help with post-processing, although can also help with
  !> conversion from flux tube to ballooning space. This probably
  !> really belongs in kt_grids, but placing here alongside related
  !> itleft, itright etc.
  integer, dimension(:, :), allocatable :: supercell_labels

  !> Records the number of unique supercells for each ky,
  integer, dimension(:), allocatable :: n_supercells

  !> Records the number of cells which are members of the supercell
  !> which owns the particular {it, ik} point. One would expect this
  !> to be equivalent to l_links + r_links + 1.
  integer, dimension(:, :), allocatable :: n_cells

  ! For use in get_fields_direct_from_dfn
  real, dimension(:, :, :), allocatable :: inv_phi_denominator_g, inv_bpar_denominator_g

  !> A type to record information about connections between iglo points
  !> Store the iglo and iproc indices for connections to the left and
  !> right of a given point.
  !> Initialised to no connections, indicated by negative values
  type :: connections_type
     integer :: iproc_left = -1,  iglo_left = -1
     integer :: iproc_right = -1, iglo_right = -1
     logical :: neighbor = .false.
  end type connections_type

  type (connections_type), dimension (:), allocatable, save :: connections
  ! (-g-layout-)

  ! linked only
  type (redist_type), save :: gc_from_left, gc_from_right
  type (redist_type), save :: links_p, links_h
  type (redist_type), save :: wfb_p, wfb_h
  type (redist_type), save :: pass_right
  type (redist_type), save :: pass_left
  type (redist_type), save :: incoming_links
  type (redist_type), save :: parity_redist

  integer, dimension (:,:), allocatable :: l_links, r_links
  integer, dimension (:,:,:), allocatable :: n_links
  logical, dimension (:,:), allocatable :: save_h

  !> Indicates if there are connected theta domains in this simulation or not.
  !> Should only be used as a part of linked boundary conditions, but worth
  !> noting that the default given may suggest that there are connections in
  !> simulations with boundary conditions other than linked.
  logical :: no_connections = .false.

  integer, dimension(:), allocatable :: M_class, N_class
  integer :: i_class

  logical :: initialized = .false.

  logical :: exb_first = .true.

  logical :: initialized_dist_fn_parameters = .false.
  logical :: initialized_dist_fn_arrays = .false.
  logical :: initialized_dist_fn_level_1 = .false.
  logical :: initialized_dist_fn_level_2 = .false.
  logical :: initialized_dist_fn_level_3 = .false.
  !logical :: initializing = .true.
  logical :: readinit = .false.
  logical :: bessinit = .false.
  logical :: connectinit = .false.
  logical :: feqinit = .false.
  logical :: lpolinit = .false.
  logical :: fyxinit = .false.
  logical :: cerrinit = .false.

  !> Arbitrary limit on the maximum number of allowed links in connected boundary conditions
  integer, parameter :: max_allowed_links = 5000

#ifdef NETCDF_PARALLEL
  logical, parameter :: moment_to_allprocs = .true.
#else
  logical, parameter :: moment_to_allprocs = .false.
#endif 

  !> Used to represent the input configuration of dist_fn
  type, extends(abstract_config_type) :: dist_fn_config_type
     ! namelist : dist_fn_knobs
     ! indexed : false
     !> The form of the adiabatic response (if a species is being modeled as adiabatic). Ignored if there are electrons in the species list.
     !>
     !> - 'no-field-line-average-term'  Adiabatic species has n = Phi.  Appropriate for single-species ETG simulations.
     !> - 'default'  Same as 'no-field-line-average-term'
     !> - 'iphi00=0' Same as 'no-field-line-average-term'
     !> - 'iphi00=1' Same as 'no-field-line-average-term'
     !> - 'field-line-average-term'  Adiabatic species has n=Phi-< Phi >.  Appropriate for single-species ITG simulations.
     !> - 'iphi00=2' Same as field-line-average-term'
     !> - 'iphi00=3' Adiabatic species has n=Phi-< Phi >_y.  Incorrect implementation of field-line-average-term.
     !>
     !> @todo Remove 'iphi00=3'
     character(len = 30) :: adiabatic_option = 'default'
     !> For debugging only.
     !>
     !> @todo Improve documentation
     real :: afilter = 0.0
     !> Leave as unity.  For debugging.
     !>
     !> @todo Improve documentation
     real :: apfac = 1.0
     !> Sets the boundary condition along the field line (i.e. the
     !> boundary conditions at \(\theta = \pm \pi\) in flux-tube
     !> simulations or \(\theta = \pm (2*\textrm{nperiod}-1)*\pi\) in
     !> ballooning space). Possible values are:
     !>
     !> - 'default' - "Smart" default -- equivalent to 'zero' except if
     !>                kt_grids:grid_option='box' in which case equivalent to 'linked'.
     !>                In simulations with suffifienctly small shear default boundaries correspond
     !>                to 'periodic'.
     !> - 'zero', 'unconnected' - Use Kotschenreuther boundary condition at endpoints
     !>                           of theta grid
     !> - 'self-periodic', 'periodic', 'kperiod=1' - Each mode is periodic in theta with
     !>                                              itself.
     !> - 'linked' - Twist and shift boundary conditions (used for kt_grids:grid_option='box')
     !> - 'alternate-zero' - Ignored
     !>
     !> See also [[dist_fn_knobs:nonad_zero]]
     character(len = 20) :: boundary_option = 'default'
     !> Overrides the [[theta_grid:itor_over_B]] internal parameter,
     !> meant only for slab runs where it sets the angle between the
     !> magnetic field and the toroidal flow.
     real :: btor_slab = 0.0
     !> True only allows solutions of single parity as determined by
     !> the input [[dist_fn_knobs:even]].
     logical :: def_parity = .false.
     !> Leave as unity for physical runs can be used for
     !> debugging. Multiplies the passing particle drifts (also see
     !> [[dist_fn_knobs:tpdriftknob]]).
     real :: driftknob = 1.0
     !> If `esv=.true.` and `boundary_option='linked'` (i.e. flux tube
     !> simulation) then at every time step we ensure the fields are
     !> exactly single valued by replacing the field at one of the
     !> repeated boundary points with the value at the other boundary
     !> (i.e. of the two array locations which should be storing the
     !> field at the same point in space we pick one and set the other
     !> equal to the one we picked). This can be important in
     !> correctly tracking the amplitude of damped modes to very small
     !> values. Also see [[init_g_knobs:clean_init]].
     logical :: esv = .false.
     !> If `def_parity` is true, determines allowed parity.
     logical :: even = .true.
     !> \(\frac{\rho}{q} \frac{d\Omega^{\rm GS2}}{d\rho_n}\) where
     !> \(\Omega^{\rm GS2}\) is toroidal angular velocity normalised
     !> to the reference frequency \(v_{t}^{\rm ref}/a\) and
     !> \(\rho_n\) is the normalised flux label which has value
     !> \(\rho\) on the local surface.
     real :: g_exb = 0.0
     !> With flow shear in single mode, constrains relative error in `phi^2`.
     !>
     !> @todo Add more documentation
     real :: g_exb_error_limit = 0.0
     !> Flow shear switched on when simulation time exceeds this time.
     real :: g_exb_start_time = -1
     !> Flow shear switched on when simulation timestep exceeds this timestep.
     integer :: g_exb_start_timestep = -1
     !> Multiplies `g_exb` in the perpendicular shearing term *but
     !> not* in the parallel drive term. Can be used to construct simulations with
     !> purely parallel flow.
     real :: g_exbfac = 1.0
     !> Perform velocity space integration using `gf_lo`, rather than
     !> `g_lo`, data decomposition if true.  If used without
     !> `field_option = 'gf_local'` in [[fields_knobs]] it is likely to give
     !> poor performance. If `field_option = 'gf_local'` in
     !> [[fields_knobs]] then this needs to be present and set to `.true.`.
     !>
     !> @todo Consider if this should be a field input.
     logical :: gf_lo_integrate = .false.
     !> Affects boundary condition at end of theta grid.   Recommended value: 1.
     !>_
     !> @note This used to default to 5.e4
     !> @todo Improve documentation
     real :: gridfac = 1.0
     !> Determines if we include the hyperviscosity term during the
     !> initialisation of the response matrix or not.
     logical :: hyper_in_initialisation = .true.
     !> Set up lowflow term \(F_1/F_0 = xi*(v/v_{th})*(a + b*(v/v_{th})^2)\)
     logical :: lf_decompose = .false.
     !> Treatment of radial variation of temperature in NEO energy variable:
     !>
     !> - (true, default) take it into account when constructing distribution
     !>   function, i.e. `H(EGS2, r)`
     !> - (false) construct `H(ENEO)` and deal with it later by including a temp
     !>   gradient term in `wstarfac` in [[dist_fn.fpp]]
     logical :: lf_default = .true.
     !> Number multiplying the coriolis drift.
     !>
     !> @todo Expand documentation
     real :: mach = 0.0
     !> Allow different species to have different values of `bakdif` and
     !> `fexpr`.  Not allowed for nonlinear runs.
     !>
     !> @note The restriction to `.false.` for nonlinear runs isn't actually enforced.
     logical :: mult_imp = .false.
     !> If true switches on "new" parallel boundary condition where h=0 at incoming boundary instead of g=0.
     !> This is considered the correct boundary condition but old behaviour retained for experimentation.
     logical :: nonad_zero = .true.
     !> Factor multiplying the parallel shearing drive term when
     !> running with non-zero [[dist_fn_knobs:g_exb]]
     real :: omprimfac = 1.0
     !>  If true then use an optimised linear source calculation which
     !>  uses pre-calculated coefficients, calculates both sigma
     !>  together and skips work associated with empty fields. Can
     !>  contribute at least 10-25% savings for linear electrostatic
     !>  collisionless simulations. For more complicated runs the
     !>  savings may be less. If enabled memory usage will
     !>  increase due to using an additional array of size 2-4 times
     !>  `gnew`.
     logical :: opt_source = .false.
     !> If non-zero, quasineutrality is not enforced,
     !> `poisfac`=  \((\lambda_\textrm{Debye}/\rho)^2\)
     real :: poisfac = 0.0
     !> Determines if we set gnew = 0 (if flag is false) or gnew = g
     !> in [[invert_rhs_1]] or not. This is currently considered
     !> experimental, but the intention is to change the default to true.
     logical :: start_from_previous_solution = .false.
     !> For debugging only. Multiplies the trapped particle drifts
     !> (also see [[dist_fn_knobs:driftknob]]).
     real :: tpdriftknob = -9.9e9
     !> For debugging only. Scales the calculated vparallel.
     !>
     !> @note If set to zero then the homogeneous contribution, `g1`, will be
     !> exactly 1 everywhere it is defined. This can lead to a divide by zero
     !> in the trapped particle continuity calculation in [[invert_rhs_1]],
     !> leading to NaNs appearing in the solution.
     real :: vparknob = 1.0
     !> For debugging only. Sets the boundary value for the barely trapped/passing particle.
     real :: wfb = 1.0
     !> Deprecated option, moved to [[le_grids]]. Setting
     !> will lead to an error message and abort. To be removed
     !> entirely in future release
     character(len = 20) :: wfbbc_option = 'unset'
     !> If true then force `gnew=0` in the forbidden region at the end
     !> of [[invert_rhs_1]] (this is the original behaviour).
     !> Nothing should depend on the forbidden region so `g` should be 0 here
     !> and if it is not for some reason then it shouldn't impact on
     !> results. If the results of your simulation depend upon this
     !> flag then something has likely gone wrong somewhere.
     logical :: zero_forbid = .false.
   contains
     procedure, public :: read => read_dist_fn_config
     procedure, public :: write => write_dist_fn_config
     procedure, public :: reset => reset_dist_fn_config
     procedure, public :: broadcast => broadcast_dist_fn_config
     procedure, public, nopass :: get_default_name => get_default_name_dist_fn_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_dist_fn_config
  end type dist_fn_config_type
  
  type(dist_fn_config_type) :: dist_fn_config
  
  !> Used to represent the input configuration of source
  type, extends(abstract_config_type) :: source_config_type
     ! namelist : source_knobs
     ! indexed : false
     !> Growth rate of non-standard source used with `source_option = 'phiext_full'`.
     real :: gamma0 = 0.0
     !> Frequency of non-standard source used with `source_option = 'phiext_full'`.
     real :: omega0 = 0.0
     !> Amplitude of external phi added as source term with
     !> `source_option = 'phiext_full'`. If zero (default) then no
     !> extra term included in the source.
     real :: phi_ext = 0.0
     !> Amplitude of non-standard source used with `source_option =
     !> 'phiext_full'` when time >= t0.
     real :: source0 = 1.0
     !> Controls the source term used in the time advance. Should be
     !> one of:
     !>
     !> - 'source_option_full' : Solve GK equation in standard form (with no artificial sources)
     !> - 'default' : Same as 'source_option_full'
     !> - 'phiext_full' : Solve GK equation with additional source
     !>   proportional to `phi_ext*F_0`.
     !> - 'homogeneous' : Solve the homogeneous equation, i.e. no potential related sources.
     character(len = 20) :: source_option = 'default'
     !> Turn on any artificial sources after time = t0. Only used with
     !> `source_option = 'phiext_full'`.
     real :: t0 = 100.0
   contains
     procedure, public :: read => read_source_config
     procedure, public :: write => write_source_config
     procedure, public :: reset => reset_source_config
     procedure, public :: broadcast => broadcast_source_config
     procedure, public, nopass :: get_default_name => get_default_name_source_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_source_config
  end type source_config_type

  type(source_config_type) :: source_config

  !> Used to represent the input configuration of
  !> dist_fn_species. There should be one of this namelist for each
  !> species simulated.
  type, extends(abstract_config_type) :: dist_fn_species_config_type
     ! namelist : dist_fn_species_knobs
     ! indexed : true
     !> Spatial implicitness parameter. Any value greater than 0 adds
     !> numerical dissipation which is often necessary to avoid grid
     !> scale oscillations. When `bakdif = 0.0` (default) we use a
     !> \(2^\textrm{nd}\) order grid centered approach in the parallel
     !> direction. When `bakdif = 1.0` this becomes a fully upwinded
     !> method. Recommended value is 0.05 for electrostatic
     !> simulations and 0.0 for electromagnetic.
     !>
     !>
     !> @warning It is possible to have different values for the
     !> different species simulated, but in a number of places we only
     !> use the value given by the first species. It is strongly
     !> recommended to use the same value for all species.
     !>
     !> @todo Clarify the motivation for the electromagnetic
     !> recommendation.
     !>
     !> @todo Consider forcing this to be constant across all species.
     !>
     !> @todo Consider changing the default to the recommended value.     
     real :: bakdif = 0.0
     !> Sets the imaginary part of the temporal implicitness
     !> parameter. It is considered an error for this to be non-zero.
     !>
     !> @todo Consider removing this variable as it should always be
     !> zero. Is this used for numerical testing?
     real :: fexpi = 0.0
     !> Sets the real part of the temporal implicitness parameter. Any
     !> value smaller than 0.5 adds numerical dissipation. When `fexpr
     !> = 0.5` we have a \(2^\textrm{nd}\) order time centered
     !> approach. If `fexpr = 0.0` this reduces to a fully implicit
     !> backward Euler method. When `fexpr = 1.0` this instead becomes
     !> a fully explicit forward Euler method (not recommended). The
     !> recommended value is 0.48.
     !>
     !> @todo Consider changing the default to the recommended value.
     real :: fexpr = 0.4
   contains
     procedure, public :: read => read_dist_fn_species_config
     procedure, public :: write => write_dist_fn_species_config
     procedure, public :: reset => reset_dist_fn_species_config
     procedure, public :: broadcast => broadcast_dist_fn_species_config
     procedure, public, nopass :: get_default_name => get_default_name_dist_fn_species_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_dist_fn_species_config
  end type dist_fn_species_config_type

  type(dist_fn_species_config_type), dimension(:), allocatable :: dist_fn_species_config
  
contains

  !> FIXME : Add documentation
  subroutine check_dist_fn(report_unit)
    use kt_grids, only: grid_option, gridopt_box, gridopt_switch
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_on
    use species, only: spec, nspec, has_electron_species
    implicit none
    integer, intent(in) :: report_unit
    integer :: is 
    if (gridfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected gridfac = ',e11.4,' in dist_fn_knobs.')") gridfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is gridfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (apfac /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected apfac = ',e11.4,' in dist_fn_knobs.')") apfac
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('The normal choice is apfac = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (driftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected driftknob = ',e11.4,' in dist_fn_knobs.')") driftknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE DRIFTS.')") 
       write (report_unit, fmt="('The normal choice is driftknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (tpdriftknob /= 1.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected tpdriftknob = ',e11.4,' in dist_fn_knobs.')") tpdriftknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE TRAPPED PARTICLE DRIFTS (either via driftknob or via tpdriftknob).')") 
       write (report_unit, fmt="('The normal choice is tpdriftknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (vparknob /= 1.0) then
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('You selected vparknob = ',e11.4,' in dist_fn_knobs.')") vparknob
       write (report_unit, fmt="('THIS IS EITHER AN ERROR, or you are DELIBERATELY SCALING THE PARALLEL VELOCITY.')")
       write (report_unit, fmt="('The normal choice is vparknob = 1.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    end if

    select case (boundary_option_switch)
    case (boundary_option_linked)
       write (report_unit, *) 
       if (gridopt_switch /= gridopt_box) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Linked boundary conditions require a box for a simulation domain.')")
          write (report_unit, fmt="('You have grid_option = ',a)") trim(grid_option)
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       else
          write (report_unit, *) 
          write (report_unit, fmt="('Linked (twist and shift) boundary conditions will be used.')")
          write (report_unit, *) 
          if (esv) then 
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, fmt="('Single valued antot arrays will be enforced.')")
             write (report_unit, fmt="('This can significantly increase the cost of the run.')")
             write (report_unit, fmt="('################# WARNING #######################')")
             write (report_unit, *) 
          endif
       end if
    case (boundary_option_self_periodic)
       write (report_unit, *) 
       write (report_unit, fmt="('Periodic boundary conditions will be used.')")
       write (report_unit, fmt="('(No twist and shift.)')")
       write (report_unit, *) 
    case default
       write (report_unit, *) 
       write (report_unit, fmt="('Outgoing boundary conditions will be used.')")
    end select

    write (report_unit, fmt="('Parallel bc for passing particles at ends of the domain is:')")
    if (nonad_zero) then
       write (report_unit, fmt="(T20,'g_wesson = g_krt = 0')")
       write (report_unit, fmt="('ie NO incoming particles in the nonadiabatic piece of delta(f)')") 
    else
       write (report_unit, fmt="(T20,'g_gs2 = 0')")
       write (report_unit, fmt="('NB this ONLY gives NO incoming particles in the nonadiabatic piece of delta(f)')")
       write (report_unit, fmt="(T20,'if phi and bpar are zero at the ends of the domain')") 
    endif
    write (report_unit, *) 

    if (.not. has_electron_species(spec)) then
       select case (adiabatic_option_switch)
       case (adiabatic_option_default)
          write (report_unit, *) 
          write (report_unit, fmt="('The adiabatic electron response is of the form:')")
          write (report_unit, *) 
          write (report_unit, fmt="('             ne = Phi')")
          write (report_unit, *) 
          write (report_unit, fmt="('This is appropriate for an ETG simulation,')") 
          write (report_unit, fmt="('where the role of ions and electrons in GS2 is switched.')")
          write (report_unit, *) 

       case (adiabatic_option_fieldlineavg)
          write (report_unit, *) 
          write (report_unit, fmt="('The adiabatic electron response is of the form:')")
          write (report_unit, *) 
          write (report_unit, fmt="('             ne = Phi - <Phi>')")
          write (report_unit, *) 
          write (report_unit, fmt="('The angle brackets denote a proper field line average.')") 
          write (report_unit, fmt="('This is appropriate for an ITG simulation.')") 
          write (report_unit, *) 

       case (adiabatic_option_yavg)
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('The adiabatic electron response is of the form:')")
          write (report_unit, *) 
          write (report_unit, fmt="('             ne = Phi - <Phi>_y')")
          write (report_unit, *) 
          write (report_unit, fmt="('The angle brackets denote an average over y only.')") 
          write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
          write (report_unit, fmt="('Perhaps you want field-line-average-term for adiabatic_option.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end select
    end if

    if (poisfac /= 0.) then
       write (report_unit, *) 
       write (report_unit, fmt="('Quasineutrality is not enforced.  The ratio (lambda_Debye/rho)**2 = ',e11.4)") poisfac
       write (report_unit, *) 
    end if

    if (mult_imp .and. nonlinear_mode_switch == nonlinear_mode_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('For nonlinear runs, all species must use the same values of fexpr and bakdif')")
       write (report_unit, fmt="('in the dist_fn_species_knobs_x namelists.')")
       write (report_unit, fmt="('THIS IS AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (def_parity .and. nonlinear_mode_switch == nonlinear_mode_on) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('Choosing a definite parity for a nonlinear run has never been tested.')")
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if

    if (def_parity) then
       if (even) then
          write (report_unit, fmt="('Only eigenmodes of even parity will be included.')")
       else
          write (report_unit, fmt="('Only eigenmodes of odd parity will be included.')")
       end if
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 
    write (report_unit, fmt="('The ExB parameter is ',f7.4)") g_exb
    if (abs(g_exb) .gt. epsilon(0.0)) then
       write (report_unit, fmt="('Perp shear terms will be multiplied by factor',f7.4)") g_exbfac
       write (report_unit, fmt="('Parallel shear term will be multiplied by factor',f7.4)") omprimfac
    endif

    write (report_unit, *) 
    write (report_unit, fmt="('------------------------------------------------------------')")
    write (report_unit, *) 

    select case (source_option_switch)

    case (source_option_full)
       write (report_unit, *) 
       write (report_unit, fmt="('The standard GK equation will be solved.')")
       write (report_unit, *) 

    case(source_option_phiext_full)
       write (report_unit, *) 
       write (report_unit, fmt="('The standard GK equation will be solved,')")
       write (report_unit, fmt="('with an additional source proportional to Phi*F_0')")
       write (report_unit, fmt="('Together with phi_ext = -1., this is the usual way to &
            & calculate the Rosenbluth-Hinton response.')")
       write (report_unit, *)

    case (source_option_homogeneous)
       write (report_unit, *)
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The homogeneous GK equation will be solved.')")
       write (report_unit, fmt="('The calculated response matrices and resulting')")
       write (report_unit, fmt="('potentials will not be accurate.')")
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *)
    end select

    ! 
    ! implicitness parameters
    !
    do is = 1, nspec
       if (aimag(fexp(is)) /= 0.) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('Species ',i2,' has fexpi = ',e11.4)") is, aimag(fexp(is))
          write (report_unit, fmt="('THIS IS AN ERROR')")
          write (report_unit, fmt="('fexpi should be zero for all species.')")
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if

       write (report_unit, fmt="('Species ',i2,' has fexpr = ', e11.4)") is, real(fexp(is))
    end do
  end subroutine check_dist_fn

  !> FIXME : Add documentation  
  subroutine wnml_dist_fn(unit)
    use species, only: spec, has_electron_species
    implicit none
    integer, intent(in) :: unit
    if (dfexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "dist_fn_knobs"

       select case (boundary_option_switch)

       case (boundary_option_default)
          write (unit, fmt="(' boundary_option = ',a)") '"default"'

       case (boundary_option_zero)
          write (unit, fmt="(' boundary_option = ',a)") '"zero"'

       case (boundary_option_self_periodic)
          write (unit, fmt="(' boundary_option = ',a)") '"periodic"'

       case (boundary_option_linked)
          write (unit, fmt="(' boundary_option = ',a)") '"linked"'

       case (boundary_option_alternate_zero)
          write (unit, fmt="(' boundary_option = ',a)") '"alternate-zero"'

       end select

       write (unit, fmt="(' nonad_zero = ',L1)") nonad_zero
       write (unit, fmt="(' gridfac = ',e17.10)") gridfac
       write (unit, fmt="(' esv = ',L1)") esv
       write (unit, fmt="(' opt_source = ',L1)") opt_source

       if (.not. has_electron_species(spec)) then
          select case (adiabatic_option_switch)

          case (adiabatic_option_default)
             write (unit, *)
             write (unit, fmt="(' adiabatic_option = ',a)") &
                  & '"no-field-line-average-term"'

          case (adiabatic_option_fieldlineavg)
             write (unit, fmt="(' adiabatic_option = ',a)") '"field-line-average-term"'

          case (adiabatic_option_yavg)
             write (unit, fmt="(' adiabatic_option = ',a)") '"iphi00=3"'

          end select
       end if

       if (apfac /= 1.) write (unit, fmt="(' apfac = ',e17.10)") apfac
       if (driftknob /= 1.) write (unit, fmt="(' driftknob = ',e17.10)") driftknob
       if (tpdriftknob /= 1.) write (unit, fmt="(' tpdriftknob = ',e17.10)") tpdriftknob
       if (poisfac /= 0.) write (unit, fmt="(' poisfac = ',e17.10)") poisfac
       if (afilter /= 0.) write (unit, fmt="(' afilter = ',e17.10)") afilter
       if (mult_imp) write (unit, fmt="(' mult_imp = ',L1)") mult_imp
       if (def_parity) then
          write (unit, fmt="(' def_parity = ',L1)") def_parity
          if (even) write (unit, fmt="(' even = ',L1)") even
       end if
       write (unit, fmt="(' /')")
    endif
    if (skexist) then
       write (unit, *)
       write (unit, fmt="(' &',a)") "source_knobs"
       select case (source_option_switch)

       case (source_option_full)
          write (unit, fmt="(' source_option = ',a)") '"full"'

       case(source_option_phiext_full)
          write (unit, fmt="(' source_option = ',a)") '"phiext_full"'
          write (unit, fmt="(' source0 = ',e17.10)") source0
          write (unit, fmt="(' omega0 = ',e17.10)") omega0
          write (unit, fmt="(' gamma0 = ',e17.10)") gamma0
          write (unit, fmt="(' t0 = ',e17.10)") t0
          write (unit, fmt="(' phi_ext = ',e17.10)") phi_ext
       case(source_option_homogeneous)
          write (unit, fmt="(' source_option = ',a)") '"homogeneous"'
       end select
       write (unit, fmt="(' /')")
    endif
  end subroutine wnml_dist_fn

  !> FIXME : Add documentation
  subroutine wnml_dist_fn_species(unit)
    use species, only: nspec
    implicit none
    integer, intent(in) :: unit
    integer :: i
    character (100) :: line
    do i=1,nspec
       write (unit, *)
       write (line, *) i
       write (unit, fmt="(' &',a)") &
            & trim("dist_fn_species_knobs_"//trim(adjustl(line)))
       write (unit, fmt="(' fexpr = ',e13.6)") real(fexp(i))
       write (unit, fmt="(' bakdif = ',e13.6)") bkdiff(i)
    end do
  end subroutine wnml_dist_fn_species

  !> FIXME : Add documentation  
  subroutine init_dist_fn_parameters(dist_fn_config_in, source_config_in, dist_fn_species_config_in)
    use gs2_layouts, only: init_gs2_layouts
    use kt_grids, only: init_kt_grids
    use le_grids, only: init_le_grids
    use species, only: init_species
    use theta_grid, only: init_theta_grid
    implicit none
    type(dist_fn_config_type), intent(in), optional :: dist_fn_config_in
    type(source_config_type), intent(in), optional :: source_config_in
    type(dist_fn_species_config_type), intent(in), dimension(:), allocatable, optional :: dist_fn_species_config_in
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_parameters) return
    initialized_dist_fn_parameters  = .true.

    if (debug) write(6,*) "init_dist_fn: init_gs2_layouts"
    call init_gs2_layouts

    if (debug) write(6,*) "init_dist_fn: init_species"
    call init_species

    if (debug) write(6,*) "init_dist_fn: init_theta_grid"
    call init_theta_grid

    if (debug) write(6,*) "init_dist_fn: init_kt_grids"
    call init_kt_grids

    if (debug) write(6,*) "init_dist_fn: init_le_grids"
    call init_le_grids

    if (debug) write(6,*) "init_dist_fn: read_parameters"
    call read_parameters(dist_fn_config_in, source_config_in, dist_fn_species_config_in)

  end subroutine init_dist_fn_parameters

  !> FIXME : Add documentation
  subroutine init_dist_fn_arrays
    use mp, only: nproc, iproc
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    use le_grids, only: nlambda, negrid
    use gs2_time, only: init_gs2_time
    use theta_grid, only: ntgrid
    use gs2_layouts, only: init_dist_fn_layouts, init_gf_layouts
    use nonlinear_terms, only: init_nonlinear_terms
    use run_parameters, only: init_run_parameters
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_arrays) return
    initialized_dist_fn_arrays  = .true.

    call init_dist_fn_parameters

    if (debug) write(6,*) "init_dist_fn: run_parameters"
    call init_run_parameters

    if (debug) write(6,*) "init_dist_fn: gs2_time"
    call init_gs2_time

    if (debug) write(6,*) "init_dist_fn: dist_fn_layouts"
    call init_dist_fn_layouts (ntgrid, naky, ntheta0, nlambda, negrid, nspec, nproc, iproc)

    if (debug) write(6,*) "init_dist_fn: gf_layouts"
    call init_gf_layouts (ntgrid, naky, ntheta0, negrid, nlambda, nspec, nproc, iproc)

    ! TODO: I don't think we need to call init_nonlinear_terms
    ! before calling allocate_arrays here because allocate_arrays
    ! does not require any of the transform layouts. I think 
    ! we can take this line out and remove the dependency 
    ! dist_fn_arrays<=nonlinear_terms in gs2_init.rb. Anyone 
    ! see any problems with this?  EGH
    if (debug) write(6,*) "init_dist_fn: nonlinear_terms"
    call init_nonlinear_terms 

    if (debug) write(6,*) "init_dist_fn: allocate_arrays"
    call allocate_arrays
  end subroutine init_dist_fn_arrays

  !> FIXME : Add documentation
  subroutine init_dist_fn_level_1
    if (initialized_dist_fn_level_1) return
    initialized_dist_fn_level_1 = .true.

    call init_dist_fn_arrays
    call init_vperp2
    call init_vpa_vpac
    call init_bc
  end subroutine init_dist_fn_level_1

  !> FIXME : Add documentation  
  subroutine init_dist_fn_level_2
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_level_2) return
    initialized_dist_fn_level_2 = .true.

    call init_dist_fn_level_1

    if (debug) write(6,*) "init_dist_fn: init_bessel"
    call init_bessel

    if (debug) write(6,*) "init_dist_fn: init_fieldeq"
    call init_fieldeq
  end subroutine init_dist_fn_level_2

  !> FIXME : Add documentation  
  subroutine init_dist_fn
    call init_dist_fn_level_3
  end subroutine init_dist_fn

  !> FIXME : Add documentation  
  subroutine init_dist_fn_level_3
    use mp, only: finish_mp, mp_abort
    use species, only: init_species
    use collisions, only: init_collisions!, vnmult
    use hyper, only: init_hyper
    use le_grids, only: init_g2gf
#ifdef LOWFLOW
    use dist_fn_arrays, only: vpac
    use lowflow, only: init_lowflow
#endif
    implicit none
    logical,parameter:: debug=.false.

    if (initialized_dist_fn_level_3) return
    initialized_dist_fn_level_3 = .true.

    if (initialized) return
    initialized = .true.

    call init_dist_fn_level_2

    if (debug) write(6,*) "init_dist_fn: init_vpar"
    call init_vpar

    if (debug) write(6,*) "init_dist_fn: init_wdrift"
    call init_wdrift

    if (debug) write(6,*) "init_dist_fn: init_wstar"
    call init_wstar

    if (debug) write(6,*) "init_dist_fn: init_collisions"
    call init_collisions ! needs to be after init_run_parameters

#ifdef LOWFLOW
    if (debug) write(6,*) "init_dist_fn: init_lowflow"
    ! needs to be before init_invert_rhs but after
    ! init_collisions as use_le_layouts is set in
    ! collisions
    call init_lowflow(wdrift, wdriftttp, wstar, vpar, vpac, &
         lf_default, lf_decompose, driftknob, tpdriftknob)
#endif

    if (debug) write(6,*) "init_dist_fn: init_invert_rhs"
    call init_invert_rhs

    if (debug) write(6,*) "init_dist_fn: init_hyper"
    call init_hyper

    if (debug) write(6,*) "init_dist_fn: init_source_term"
    call init_source_term

    if (gf_lo_integrate) call init_g2gf(debug)

  end subroutine init_dist_fn_level_3

  !> FIXME : Add documentation
  subroutine finish_dist_fn
    call finish_dist_fn_parameters
  end subroutine finish_dist_fn

  !> FIXME : Add documentation  
  subroutine finish_dist_fn_parameters
    if (.not. initialized_dist_fn_parameters) return
    initialized_dist_fn_parameters = .false.
    call finish_dist_fn_arrays
    readinit = .false.
    if (allocated(fexp)) deallocate (fexp, bkdiff)
  end subroutine finish_dist_fn_parameters

  subroutine finish_dist_fn_arrays
    use dist_fn_arrays, only: g, gnew, g_work, kx_shift, theta0_shift
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3
    use dist_fn_arrays, only: antot, antota, antotp
    use dist_fn_arrays, only: fieldeq, fieldeqa, fieldeqp
#ifdef SHMEM
    use shm_mpi3, only : shm_free
#endif

    implicit none

    if (.not. initialized_dist_fn_arrays) return
    initialized_dist_fn_arrays = .false.

    call finish_dist_fn_level_1
    if (allocated(g)) deallocate (g, gnew, g_work)
    if (allocated(source_coeffs)) deallocate(source_coeffs)
#ifndef SHMEM
    if (allocated(gexp_1)) deallocate (gexp_1, gexp_2, gexp_3)
#else
    if (allocated(gexp_2)) deallocate (gexp_2, gexp_3)
    if (associated(gexp_1)) call shm_free(gexp_1)
#endif
    if (allocated(g_h)) deallocate (g_h)
    if (allocated(save_h)) deallocate (save_h)
    if (allocated(kx_shift)) deallocate (kx_shift)
    if (allocated(theta0_shift)) deallocate (theta0_shift)
    if (allocated(antot)) deallocate(antot)
    if (allocated(antota)) deallocate(antota)
    if (allocated(antotp)) deallocate(antotp)
    if (allocated(fieldeq)) deallocate(fieldeq)
    if (allocated(fieldeqa)) deallocate(fieldeqa)
    if (allocated(fieldeqp)) deallocate(fieldeqp)
    initialized_dist_fn_arrays = .false.
  end subroutine finish_dist_fn_arrays

  subroutine finish_dist_fn_level_1
    use redistribute, only: delete_redist
    use dist_fn_arrays, only: vperp2, vpa, vpac, vpa_gf, vperp2_gf
#ifdef LOWFLOW
    use lowflow, only: finish_lowflow_terms
#endif
    implicit none
    if (.not. initialized_dist_fn_level_1) return
    initialized_dist_fn_level_1 = .false.
    no_connections = .false.
    connectinit = .false.
    lpolinit = .false. ; fyxinit = .false. ; cerrinit = .false.
    increase = .true. ; decrease = .false.
    exb_first = .true.

    call finish_dist_fn_level_2

    if (allocated(vperp2)) deallocate (vperp2)
    if (allocated(vpa)) deallocate (vpa, vpac)
    if (allocated(vpa_gf)) deallocate (vpa_gf, vperp2_gf)

    if (allocated(l_links)) deallocate (l_links, r_links, n_links)
    if (allocated(M_class)) deallocate (M_class, N_class)
    if (allocated(itleft)) deallocate (itleft, itright)
    if (allocated(supercell_labels)) deallocate (supercell_labels)
    if (allocated(n_supercells)) deallocate (n_supercells)
    if (allocated(n_cells)) deallocate (n_cells)
    if (allocated(connections)) deallocate (connections)
    if (allocated(g_adj)) deallocate (g_adj)
    if (allocated(jump)) deallocate (jump)
    if (allocated(ikx_indexed)) deallocate (ikx_indexed)
    if (allocated(ufac)) deallocate (ufac)
    if (allocated(fl_avg)) deallocate (fl_avg)
    if (allocated(awgt)) deallocate (awgt)
    if (allocated(kmax)) deallocate (kmax)

    call delete_redist(gc_from_left)
    call delete_redist(gc_from_right)
    call delete_redist(links_p)
    call delete_redist(links_h)
    call delete_redist(wfb_p)
    call delete_redist(wfb_h)
    call delete_redist(pass_right)
    call delete_redist(pass_left)
    call delete_redist(incoming_links)
    call delete_redist(parity_redist)

#ifdef LOWFLOW
    call finish_lowflow_terms
#endif

  end subroutine finish_dist_fn_level_1

  !> FIXME : Add documentation  
  subroutine finish_dist_fn_level_2
    use dist_fn_arrays, only: aj0, aj1, aj0_gf, aj1_gf

    if (.not. initialized_dist_fn_level_2) return
    initialized_dist_fn_level_2 = .false.

    call finish_dist_fn_level_3
    bessinit = .false. 
    feqinit = .false.  
!AJ
    if (allocated(aj0)) deallocate (aj0, aj1)
    if (allocated(aj0_gf)) deallocate (aj0_gf, aj1_gf)
    if (allocated(gridfac1)) deallocate (gridfac1, gamtot, gamtot1, gamtot2)
    if (allocated(gamtot3)) deallocate (gamtot3)
    if (allocated(a)) deallocate (a, b, r, ainv)
    if (allocated(inv_phi_denominator_g)) deallocate(inv_phi_denominator_g)
    if (allocated(inv_bpar_denominator_g)) deallocate(inv_bpar_denominator_g)
  end subroutine finish_dist_fn_level_2

  !> FIXME : Add documentation  
  subroutine finish_dist_fn_level_3
    if (.not. initialized_dist_fn_level_3) return
    initialized_dist_fn_level_3 = .false.

    initialized = .false.

    if (allocated(vpar)) deallocate (vpar)
    if (allocated(wdrift)) deallocate (wdrift, wdriftttp)
    if (allocated(wstar)) deallocate (wstar)

    call dist_fn_config%reset()
    call source_config%reset()
    if(allocated(dist_fn_species_config)) deallocate(dist_fn_species_config)
  end subroutine finish_dist_fn_level_3

  !> FIXME : Add documentation  
  subroutine set_profiles_overrides(prof_ov)
    use overrides, only: profiles_overrides_type
    use unit_tests, only: debug_message
    type(profiles_overrides_type), intent(in) :: prof_ov
    call debug_message(3, 'dist_fn::set_overrides starting')
    if (prof_ov%override_g_exb) g_exb = prof_ov%g_exb
    if (prof_ov%override_mach) mach = prof_ov%mach
  end subroutine set_profiles_overrides

  !> FIXME : Add documentation  
  subroutine set_optimisations_overrides(opt_ov)
    use overrides, only: optimisations_overrides_type
    type(optimisations_overrides_type), intent(in) :: opt_ov
    if (opt_ov%override_opt_source) opt_source = opt_ov%opt_source
  end subroutine set_optimisations_overrides

  !> FIXME : Add documentation  
  subroutine read_parameters(dist_fn_config_in, source_config_in, dist_fn_species_config_in)
    use file_utils, only: error_unit
    use theta_grid, only: is_effectively_zero_shear
    use text_options, only: text_option, get_option_value
    use species, only: nspec
    use mp, only: proc0, broadcast, mp_abort
    use theta_grid, only: itor_over_B
    use kt_grids, only: gridopt_box, gridopt_switch
    implicit none
    type(dist_fn_config_type), intent(in), optional :: dist_fn_config_in
    type(source_config_type), intent(in), optional :: source_config_in
    type(dist_fn_species_config_type), intent(in), dimension(:), allocatable, optional :: dist_fn_species_config_in
    type (text_option), dimension (4), parameter :: sourceopts = &
         (/ text_option('default', source_option_full), &
            text_option('full', source_option_full), &
            text_option('phiext_full', source_option_phiext_full), &
            text_option('homogeneous', source_option_homogeneous) /)
    character(20) :: source_option

    type (text_option), dimension (8), parameter :: boundaryopts = &
         (/ text_option('default', boundary_option_default), &
            text_option('zero', boundary_option_zero), &
            text_option('unconnected', boundary_option_zero), &
            text_option('self-periodic', boundary_option_self_periodic), &
            text_option('periodic', boundary_option_self_periodic), &
            text_option('kperiod=1', boundary_option_self_periodic), &
            text_option('linked', boundary_option_linked), &
            text_option('alternate-zero', boundary_option_alternate_zero) /)
    character(20) :: boundary_option

    type (text_option), dimension (7), parameter :: adiabaticopts = &
         (/ text_option('default', adiabatic_option_default), &
            text_option('no-field-line-average-term', adiabatic_option_default), &
            text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
! eventually add in iphi00 = 0 option:
            text_option('iphi00=0', adiabatic_option_default), &
            text_option('iphi00=1', adiabatic_option_default), &
            text_option('iphi00=2', adiabatic_option_fieldlineavg), &
            text_option('iphi00=3', adiabatic_option_yavg)/)
    character(30) :: adiabatic_option
    character(20) :: wfbbc_option

    integer :: ierr, is
    real :: bd

    if (readinit) return
    readinit = .true.


    !Deal with source_knobs namelist
    if (present(source_config_in)) source_config = source_config_in

    call source_config%init(name = 'source_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    gamma0 = source_config%gamma0
    omega0 = source_config%omega0
    phi_ext = source_config%phi_ext
    source0 = source_config%source0
    source_option = source_config%source_option
    t0 = source_config%t0

    skexist = source_config%exist

    !Deal with dist_fn_knobs namelist    
    if (present(dist_fn_config_in)) dist_fn_config = dist_fn_config_in

    call dist_fn_config%init(name = 'dist_fn_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    adiabatic_option = dist_fn_config%adiabatic_option
    afilter = dist_fn_config%afilter
    apfac = dist_fn_config%apfac
    boundary_option = dist_fn_config%boundary_option
    btor_slab = dist_fn_config%btor_slab
    def_parity = dist_fn_config%def_parity
    driftknob = dist_fn_config%driftknob
    esv = dist_fn_config%esv
    even = dist_fn_config%even
    g_exb = dist_fn_config%g_exb
    g_exb_error_limit = dist_fn_config%g_exb_error_limit
    g_exb_start_time = dist_fn_config%g_exb_start_time
    g_exb_start_timestep = dist_fn_config%g_exb_start_timestep
    g_exbfac = dist_fn_config%g_exbfac
    gf_lo_integrate = dist_fn_config%gf_lo_integrate
    gridfac = dist_fn_config%gridfac
    hyper_in_initialisation = dist_fn_config%hyper_in_initialisation
    lf_decompose = dist_fn_config%lf_decompose
    lf_default = dist_fn_config%lf_default
    mach = dist_fn_config%mach
    mult_imp = dist_fn_config%mult_imp
    nonad_zero = dist_fn_config%nonad_zero
    omprimfac = dist_fn_config%omprimfac
    opt_source = dist_fn_config%opt_source
    poisfac = dist_fn_config%poisfac
    start_from_previous_solution = dist_fn_config%start_from_previous_solution
    tpdriftknob = dist_fn_config%tpdriftknob
    vparknob = dist_fn_config%vparknob
    wfb = dist_fn_config%wfb
    wfbbc_option = dist_fn_config%wfbbc_option
    zero_forbid = dist_fn_config%zero_forbid

    dfexist = dist_fn_config%exist

    if (trim(adjustl(wfbbc_option)) /= 'unset') then
       call mp_abort("The input wfbbc_option has been set in 'dist_fn_knobs' but this "// &
            "input has moved to 'le_grids_knobs'", .true.)
    end if
    
    if (tpdriftknob == -9.9e9) tpdriftknob=driftknob

    ierr = error_unit()

    call get_option_value &
         (boundary_option, boundaryopts, boundary_option_switch, &
         ierr, "boundary_option in dist_fn_knobs",.true.)

    ! Handle a request for default boundaries
    if (boundary_option_switch == boundary_option_default) then
       ! Linked for box, zero otherwise
       if (gridopt_switch == gridopt_box) then
          boundary_option_switch = boundary_option_linked
          boundary_option = 'linked'
       else
          boundary_option_switch = boundary_option_zero
          boundary_option = 'zero'
       end if

       ! Fully periodic at low enough shear
       if(is_effectively_zero_shear()) then
          boundary_option_switch = boundary_option_self_periodic
          boundary_option = 'periodic'
       end if
    end if

    call get_option_value &
         (source_option, sourceopts, source_option_switch, &
         ierr, "source_option in source_knobs",.true.)
    call get_option_value &
         (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
         ierr, "adiabatic_option in dist_fn_knobs",.true.)
       
    if (.not.allocated(fexp)) allocate (fexp(nspec), bkdiff(nspec))

    call read_species_knobs(dist_fn_species_config_in)

    call broadcast (fexp)
    call broadcast (bkdiff)
  
    !<DD>Turn off esv if not using linked boundaries
    esv=esv.and.(boundary_option_switch.eq.boundary_option_linked)

    if (mult_imp) then
       ! nothing -- fine for linear runs, but not implemented nonlinearly
    else
! consistency check for bkdiff
       bd = bkdiff(1)
       do is = 1, nspec
          if (bkdiff(is) /= bd) then
             if (proc0) write(*,*) 'Forcing bkdiff for species ',is,' equal to ',bd
             if (proc0) write(*,*) 'If this is a linear run, and you want unequal bkdiff'
             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
             bkdiff(is) = bd
          endif
       end do
! consistency check for fexp
!       fe = fexp(1)
!       do is = 1, nspec
!          if (fexp(is) /= fe) then
!             if (proc0) write(*,*) 'Forcing fexp for species ',is,' equal to ',fe
!             if (proc0) write(*,*) 'If this is a linear run, and you want unequal fexp'
!             if (proc0) write(*,*) 'for different species, specify mult_imp = .true.'
!             if (proc0) write(*,*) 'in the dist_fn_knobs namelist.'
!             fexp(is) = fe
!          endif
!       end do
    end if

! consistency check for afilter
!    if (afilter /= 0.0) then
!       if (proc0) write(*,*) 'Forcing afilter = 0.0'
!       afilter = 0.0
!    end if


!CMR, 15/2/11:  Move following lines to here from init_dist_fn, so read_parameters 
!               sets up itor_over_B
!!CMR, 19/10/10:
!! Override itor_over_B, if "dist_fn_knobs" parameter btor_slab ne 0
!! Not ideal to set geometry quantity here, but its historical! 
       if (abs(btor_slab) > epsilon(0.0)) itor_over_B = btor_slab
!! Done for slab, where itor_over_B is determined by angle between B-field 
!! and toroidal flow: itor_over_B = (d(u_z)/dx) / (d(u_y)/dx) = Btor / Bpol
!! u = u0 (phihat) = x d(u0)/dx (phihat) = x d(uy)/dx (yhat + Btor/Bpol zhat)
!! g_exb = d(uy)/dx => d(uz)/dx = g_exb * Btor/Bpol = g_exb * itor_over_B

  end subroutine read_parameters 

  !> FIXME : Add documentation  
  subroutine read_species_knobs(dist_fn_species_config_in)
    use species, only: nspec
    use mp, only: proc0
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    type(dist_fn_species_config_type), intent(in), dimension(:), allocatable, optional :: dist_fn_species_config_in
    integer :: is
    logical :: exist
    real :: fexpr, fexpi

    if(present(dist_fn_species_config_in)) dist_fn_species_config = dist_fn_species_config_in
    if(.not.allocated(dist_fn_species_config)) allocate(dist_fn_species_config(nspec))
    if(size(dist_fn_species_config).ne.nspec) then
       if(proc0) print*,"Warning: inconsistent config size in read_species_knobs"
    endif
    
    do is = 1, nspec
       call dist_fn_species_config(is)%init(name = 'dist_fn_species_knobs', requires_index = .true., index = is)

       ! Copy out internal values into module level parameters
       bkdiff(is) = dist_fn_species_config(is)%bakdif
       fexpi = dist_fn_species_config(is)%fexpi
       fexpr = dist_fn_species_config(is)%fexpr
       
       exist = dist_fn_species_config(is)%exist
       
       fexp(is) = cmplx(fexpr,fexpi)
    end do
  end subroutine read_species_knobs

  !> FIXME : Add documentation
  subroutine init_wdrift
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, il_is_wfb, il_is_trapped, can_be_ttp, is_ttp
    use gs2_layouts, only: g_lo, il_idx
    use array_utils, only: zero_array
    implicit none
    integer :: ig, il, iglo
    real :: scaling_knob, wdrift_value, wcoriolis_value
    logical, parameter :: debug = .false.

    if (.not. allocated(wdrift)) then
       ! allocate wdrift with sign(vpa) dependence because will contain
       ! Coriolis as well as magnetic drifts
       allocate (wdrift(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (wdriftttp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    end if
    call zero_array(wdrift) ; call zero_array(wdriftttp)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ig, scaling_knob, wdrift_value, wcoriolis_value) &
    !$OMP SHARED(g_lo, tpdriftknob, driftknob, ntgrid, forbid, wdrift) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il=il_idx(g_lo,iglo)

       ! Note here we always treat the wfb as trapped
       if (il_is_wfb(il) .or. il_is_trapped(il)) then
          scaling_knob = tpdriftknob
       else
          scaling_knob = driftknob
       end if

       do ig = -ntgrid, ntgrid
          if (forbid(ig, il)) cycle
          wdrift_value = wdrift_func(ig, iglo)
          wcoriolis_value = wcoriolis_func(ig, iglo)

          ! Add Coriolis drift to magnetic drifts. Is there a reason we don't scale
          ! this by scaling_knob?
          wdrift(ig, 1, iglo) = wdrift_value * scaling_knob + wcoriolis_value
          wdrift(ig, 2, iglo) = wdrift_value * scaling_knob - wcoriolis_value
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ig, scaling_knob, wdrift_value, wcoriolis_value) &
    !$OMP SHARED(g_lo, tpdriftknob, ntgrid, can_be_ttp, wdriftttp, is_ttp) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il=il_idx(g_lo,iglo)
       ! Cycle if this il is never ttp as we only care about ttp here.
       if (.not. can_be_ttp(il)) cycle
       do ig = -ntgrid, ntgrid
          ! The below cycle will activate regularly, indicating wdriftttp
          ! may be a relatively inefficent way of storing this
          ! information, requiring a lot of memory to hold mostly 0.
          if (.not. is_ttp(ig, il)) cycle

          !GWH+JAB: should this be calculated only at ittp? or for
          !each totally trapped pitch angle? (Orig logic: there was
          !only one totally trapped pitch angle; now multiple ttp are
          !allowed). Previously passed ittp(ig) instead of il. Now we pass
          !iglo instead of il (although this is equivalent here).
          wdrift_value = wdrift_func(ig, iglo)
          wcoriolis_value = wcoriolis_func(ig, iglo)

          !CMR:  totally trapped particle drifts also scaled by tpdriftknob
          ! add Coriolis drift to magnetic drifts. Is there a reason we don't scale
          ! this by scaling_knob?
          wdriftttp(ig, 1, iglo) = wdrift_value * tpdriftknob + wcoriolis_value
          wdriftttp(ig, 2, iglo) = wdrift_value * tpdriftknob - wcoriolis_value
       end do
    end do
    !$OMP END PARALLEL DO

    ! This should be weighted by bakdif to be completely consistent
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ig) &
    !$OMP SHARED(g_lo, ntgrid, wdrift) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do ig = -ntgrid, ntgrid-1
          wdrift(ig,:,iglo) = 0.5*(wdrift(ig,:,iglo) + wdrift(ig+1,:,iglo))
       end do
    end do
    !$OMP END PARALLEL DO

    if (debug) write(6,*) 'init_wdrift: driftknob, tpdriftknob=',driftknob,tpdriftknob

  end subroutine init_wdrift

  !> Calculates the value of the magnetic drifts (curvature + grad-B)
  !> corresponding to the theta grid location given by `ig` for the iglo
  !> index of interest.
  pure real function wdrift_func (ig, iglo)
    use dist_fn_arrays, only: vperp2, vpa
    use theta_grid, only: gbdrift, gbdrift0, cvdrift, cvdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use gs2_time, only: code_dt, wunits
    use gs2_layouts, only: g_lo, it_idx, ik_idx
    implicit none
    integer, intent (in) :: ig, iglo
    integer :: it, ik
    real :: vperp2_local, vpa2_local

    vperp2_local = vperp2(ig, iglo)
    vpa2_local = vpa(ig, 1, iglo)**2

    it=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)

    ! note that wunits=aky/2 (for wstar_units=F)
    if (aky(ik) == 0.0) then
       wdrift_func = akx(it)/shat &
                    *(cvdrift0(ig)*vpa2_local &
                      + gbdrift0(ig)*0.5*vperp2_local) &
                     *code_dt/2.0
    else
       wdrift_func = ((cvdrift(ig) + theta0(it,ik)*cvdrift0(ig)) &
                        *vpa2_local &
                      + (gbdrift(ig) + theta0(it,ik)*gbdrift0(ig)) &
                        *0.5*vperp2_local) &
                     *code_dt*wunits(ik)
    end if
  end function wdrift_func

  !> Calculates the value of the Coriolis drift at the location
  !> corresponding to the passed indices.
  pure real function wcoriolis_func (ig, iglo)
    use dist_fn_arrays, only: vpa
    use theta_grid, only: cdrift, cdrift0, shat
    use kt_grids, only: aky, theta0, akx
    use gs2_time, only: code_dt, wunits
    use species, only: spec
    use gs2_layouts, only: g_lo, it_idx, ik_idx, is_idx
    implicit none
    integer, intent (in) :: ig, iglo
    integer :: it, ik, is

    it=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)
    is=is_idx(g_lo,iglo)

    if (aky(ik) == 0.0) then
       wcoriolis_func = mach*vpa(ig,1,iglo) &
            * cdrift0(ig) * code_dt * akx(it)/(2.*shat*spec(is)%stm)
    else
       wcoriolis_func = mach*vpa(ig,1,iglo) &
            * (cdrift(ig) + theta0(it,ik)*cdrift0(ig))*code_dt*wunits(ik)/spec(is)%stm
    end if

  end function wcoriolis_func

  !> FIXME : Add documentation  
  subroutine init_vperp2
    use theta_grid, only: ntgrid, bmag
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    use dist_fn_arrays, only: vperp2
    use le_grids, only: energy, al, forbid
    use array_utils, only: zero_array
    implicit none
    real :: al1, e1
    integer :: iglo, is, il
    if (.not.allocated(vperp2)) then
       allocate (vperp2 (-ntgrid:ntgrid,  g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    call zero_array(vperp2)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, al1, is, e1, il) &
    !$OMP SHARED(g_lo, al, bmag, energy, vperp2, forbid) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       al1 = al(il)
       is = is_idx(g_lo,iglo)
       e1 = energy(ie_idx(g_lo,iglo),is)

       where (.not. forbid(:, il))
          vperp2(:,iglo) = bmag*al1*e1
       end where
    end do
    !$OMP END PARALLEL DO
  end subroutine init_vperp2

  !> Initialised the on grid v|| (vpa) and grid-centred v|| (vpac) arrays
  subroutine init_vpa_vpac
    use dist_fn_arrays, only: vpa, vpac
    use dist_fn_arrays, only: vpa_gf, vperp2_gf
    use species, only: nspec
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy, al, negrid, nlambda, forbid
    use gs2_layouts, only: g_lo, ik_idx, il_idx, ie_idx, is_idx
    use array_utils, only: zero_array
    implicit none
    integer :: iglo, ik, is
!AJ
    integer :: il, ie, ig
    real :: al1, e1
    
    if (.not.allocated(vpa)) then
       allocate (vpa    (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (vpac   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    call zero_array(vpa) ; call zero_array(vpac)

    !AJ if we are doing gf_lo integrate then pre-calcuate some of the
    !data used in the integrate.
    if(gf_lo_integrate) then
       allocate(vpa_gf    (-ntgrid:ntgrid, 2, negrid, nlambda, nspec))
       allocate(vperp2_gf (-ntgrid:ntgrid, negrid, nlambda, nspec))
       
       do il = 1,nlambda
          do ie = 1,negrid
             do is = 1,nspec
                e1 = energy(ie,is)
                al1 = al(il)
                vpa_gf(:, 1, ie, il,is) = sqrt(e1*max(0.0, 1.0 - al1*bmag))
                vpa_gf(:, 2, ie, il,is) = -vpa_gf(:, 1, ie, il,is)
                do ig = -ntgrid,ntgrid
                   if(1.0 - al1*bmag(ig) < 100.0*epsilon(0.0)) then
                      vpa_gf(ig, :, ie, il,is) = 0.0
                   end if
                   ! Should this be moved to init_vperp2?
                   vperp2_gf(ig, ie, il,is) = e1*al1*bmag(ig)
                end do
             end do
          end do
       end do
    end if

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, al1, e1, il, ik, is) &
    !$OMP SHARED(g_lo, al, energy, vpa, bmag, vparknob, ntgrid, vpac, forbid) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo, iglo)
       al1 = al(il)
       e1 = energy(ie_idx(g_lo,iglo), is)

       vpa(:,1,iglo) = sqrt(e1*max(0.0, 1.0 - al1*bmag)) * vparknob
       vpa(:,2,iglo) = - vpa(:,1,iglo)

       where (forbid(:, il))
          vpa(:,1,iglo) = 0.0
          vpa(:,2,iglo) = 0.0
       end where

! Where vpac /= 1, it could be weighted by bakdif for better consistency??
!CMR, 4/8/2011:
!CMR : vpa is parallel velocity at grid points (normalised to v_ts)
!CMR : vpac is grid centered parallel velocity
!CMR : vpar = q_s/sqrt{T_s m_s}*DELT/DTHETA * vpac |\gradpar(theta)| 
!                                     where gradpar(theta) is centered
!  ie  vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT
! 
!   comments on vpac, vpar
!    (i) should some be weighted by bakdif?
!CMR 

       where (forbid(:ntgrid-1, il) .or. forbid(-ntgrid+1:, il))
          vpac(-ntgrid:ntgrid-1,1,iglo) = 0.0
          vpac(-ntgrid:ntgrid-1,2,iglo) = 0.0
       elsewhere
          vpac(-ntgrid:ntgrid-1,1,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,1,iglo) + vpa(-ntgrid+1:ntgrid,1,iglo))
          vpac(-ntgrid:ntgrid-1,2,iglo) = &
              0.5*(vpa(-ntgrid:ntgrid-1,2,iglo) + vpa(-ntgrid+1:ntgrid,2,iglo))
       end where
       vpac(ntgrid,:,iglo) = 0.0
    end do
    !$OMP END PARALLEL DO
  end subroutine init_vpa_vpac

  !> FIXME : Add documentation
  subroutine init_vpar
    use dist_fn_arrays, only: vpac
    use species, only: spec
    use theta_grid, only: ntgrid, delthet, gradpar
    use gs2_time, only: code_dt, tunits
    use gs2_layouts, only: g_lo, ik_idx, il_idx, is_idx
    use array_utils, only: zero_array
    implicit none
    integer :: iglo, ik, is, il

    if (.not.allocated(vpar)) then
       allocate (vpar   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    call zero_array(vpar)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ik, is) &
    !$OMP SHARED(g_lo, ntgrid, vpac, vpar, spec, tunits, code_dt, delthet, gradpar) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo, iglo)

!CMR, 4/8/2011:
!CMR : vpac is grid centered parallel velocity
!CMR : vpar = q_s/sqrt{T_s m_s}*DELT/DTHETA * vpac |\gradpar(theta)|
!                                     where gradpar(theta) is centered
!  ie  vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT
!
!   comments on vpac, vpar
!    (i) should some be weighted by bakdif?
!CMR

       ik = ik_idx(g_lo,iglo)
       vpar(-ntgrid:ntgrid-1,1,iglo) = &
            spec(is)%zstm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,1,iglo)
       vpar(-ntgrid:ntgrid-1,2,iglo) = &
            spec(is)%zstm*tunits(ik)*code_dt &
            *0.5/delthet(-ntgrid:ntgrid-1) &
            *(abs(gradpar(-ntgrid:ntgrid-1)) + abs(gradpar(-ntgrid+1:ntgrid)))&
            *vpac(-ntgrid:ntgrid-1,2,iglo)
       vpar(ntgrid,:,iglo) = 0.0

    end do
    !$OMP END PARALLEL DO
  end subroutine init_vpar

  !> FIXME : Add documentation  
  subroutine init_wstar
    use species, only: nspec, dlnf0drho
    use kt_grids, only: naky
    use le_grids, only: negrid
    use gs2_time, only: code_dt, wunits

    implicit none
    integer :: ik, ie, is

    if(.not.allocated(wstar)) allocate (wstar(naky,negrid,nspec))

    do is = 1, nspec
       do ie = 1, negrid
          do ik = 1, naky
             wstar(ik,ie,is) = - code_dt*wunits(ik) * dlnf0drho(ie,is)
          end do
       end do
    end do
  end subroutine init_wstar

  !> Compute and store the Bessel functions required for future usage,
  !> aj0 and aj1.
  !>
  !> @note j1 returns and aj1 stores J_1(x)/x (NOT J_1(x)),
  subroutine init_bessel
    use dist_fn_arrays, only: aj0, aj1, aj0_gf, aj1_gf
    use kt_grids, only: kperp2
    use species, only: spec, nspec
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: energy,  al,negrid, nlambda, forbid
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx, gf_lo
    use spfunc, only: j0, j1
    implicit none
    integer :: ik, it, il, ie, is, iglo, igf
    real, dimension(-ntgrid:ntgrid) :: arg

    if (bessinit) return
    bessinit = .true.

    allocate (aj0(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (aj1(-ntgrid:ntgrid,g_lo%llim_proc:g_lo%ulim_alloc))

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, il, ie, is, arg) &
    !$OMP SHARED(g_lo, spec, energy, al, bmag, kperp2, aj0, aj1, forbid) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       arg = spec(is)%bess_fac * spec(is)%smz * &
            sqrt(energy(ie,is) * al(il) * kperp2(:,it,ik) / bmag)
       where (forbid(:, il))
          aj0(:, iglo) = 0.0
          aj1(:, iglo) = 0.0
       elsewhere
          aj0(:, iglo) = j0(arg)
          aj1(:, iglo) = j1(arg)
       end where
    end do
    !$OMP END PARALLEL DO

    !AJ For the gf_lo integrate we pre-calculate some of the factors
    !used in the velocity space integration, rather than
    !calculating them on the fly.
    if(gf_lo_integrate) then
       allocate (aj0_gf(-ntgrid:ntgrid,nspec,negrid,nlambda,gf_lo%llim_proc:gf_lo%ulim_alloc))
       allocate (aj1_gf(-ntgrid:ntgrid,nspec,negrid,nlambda,gf_lo%llim_proc:gf_lo%ulim_alloc))

       do igf = gf_lo%llim_proc, gf_lo%ulim_proc
          ik = ik_idx(gf_lo,igf)
          it = it_idx(gf_lo,igf)
          do il = 1,gf_lo%nlambda
             do ie = 1,gf_lo%negrid
                do is = 1,gf_lo%nspec
                   arg = spec(is)%bess_fac * spec(is)%smz * &
                        sqrt(energy(ie,is) * al(il) * kperp2(:,it,ik) / bmag)
                   where (forbid(:, il))
                      aj0_gf(:, is, ie, il, igf) = 0.0
                      aj1_gf(:, is, ie, il, igf) = 0.0
                   else where
                      aj0_gf(:, is, ie, il, igf) = j0(arg)
                      aj1_gf(:, is, ie, il, igf) = j1(arg)
                   end where
                end do
             end do
          end do
       end do
    end if
  end subroutine init_bessel

  !> FIXME : Add documentation  
  subroutine init_invert_rhs
    use species, only: spec
    use theta_grid, only: ntgrid
    use le_grids, only: forbid, grid_has_trapped_particles, is_lower_bounce_point, is_ttp
    use constants, only: zi
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use array_utils, only: zero_array
#ifdef LOWFLOW
    use lowflow, only: wstar_neo
#endif
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is, isgn
    real :: wdttp, vp, bd
#ifdef LOWFLOW
    complex :: wd
#else
    real :: wd
#endif
    if (.not.allocated(a)) then
       allocate (a(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (b(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (ainv(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    endif
    call zero_array(a) ; call zero_array(b)
    call zero_array(r) ; call zero_array(ainv)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, il, ie, is, isgn, ig, wd, wdttp, vp, bd) &
    !$OMP SHARED(g_lo, ntgrid, wdrift, wdriftttp, vpar, bkdiff, ainv, fexp, &
#ifdef LOWFLOW
    !$OMP wstar_neo, &
#endif
    !$OMP spec, forbid, is_ttp, r, a, b) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid-1
             wd = wdrift(ig,isgn,iglo)
             wdttp = wdriftttp(ig, isgn, iglo)

# ifdef LOWFLOW
             wd = wd + wstar_neo(ig,it,ik)*spec(is)%zt
             wdttp = wdttp + wstar_neo(ig,it,ik)*spec(is)%zt
# endif

             ! use positive vpar because we will be flipping sign of d/dz
             ! when doing parallel field solve for -vpar
             ! also, note that vpar is (vpa bhat + v_{magnetic}) . grad theta / delthet
             ! if compile with LOWFLOW defined
             vp = vpar(ig,1,iglo) 
             bd = bkdiff(is)

             ainv(ig,isgn,iglo) &
                  = 1.0/(1.0 + bd &
                  + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0+bd) + 2.0*vp))
             r(ig,isgn,iglo) &
                  = (1.0 - bd &
                  + (1.0-fexp(is))*spec(is)%tz*(zi*wd*(1.0-bd) - 2.0*vp)) &
                  *ainv(ig,isgn,iglo)
             a(ig,isgn,iglo) &
                  = 1.0 + bd &
                  + fexp(is)*spec(is)%tz*(-zi*wd*(1.0+bd) - 2.0*vp)
             b(ig,isgn,iglo) &
                  = 1.0 - bd &
                  + fexp(is)*spec(is)%tz*(-zi*wd*(1.0-bd) + 2.0*vp)
          
             if (grid_has_trapped_particles()) then
                ! zero out forbidden regions
                if (forbid(ig,il) .or. forbid(ig+1,il)) then
                   r(ig,isgn,iglo) = 0.0
                   ainv(ig,isgn,iglo) = 0.0
                end if
             
! CMR, DD, 25/7/2014: 
!  set ainv=1 just left of lower bounce points ONLY for RIGHTWARDS travelling 
!  trapped particles. Part of multiple trapped particle algorithm
!  NB not applicable to ttp or wfb!
                if( isgn.eq.1 )then
                   if (is_lower_bounce_point(ig+1, il)) then
                      ainv(ig,isgn,iglo) = 1.0 + ainv(ig,isgn,iglo)
                   end if
                endif
                ! ???? mysterious mucking around with totally trapped particles
                ! part of multiple trapped particle algorithm
                if (is_ttp(ig, il)) then
                   ainv(ig,isgn,iglo) = 1.0/(1.0 + zi*(1.0-fexp(is))*spec(is)%tz*wdttp)
                   a(ig,isgn,iglo) = 1.0 - zi*fexp(is)*spec(is)%tz*wdttp
                   r(ig,isgn,iglo) = 0.0
                end if
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call init_homogeneous_g(g_h)
  end subroutine init_invert_rhs

  !> Calculates the homogeneous solution `g_h` used as a part of the
  !> trapped particle and linked boundary condition solve in
  !> [[invert_rhs]]. This solution only depends on `r` so can be
  !> calculated during initialisation. Note that `r` depends on
  !> physics terms that contain `dt` so we must recalculate the
  !> solution `g_h` if the timestep changes.
  subroutine init_homogeneous_g(g_h, force_passing)
    use array_utils, only: zero_array
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, ik_idx, il_idx
    use kt_grids, only: aky
    use le_grids, only: mixed_wfb, passing_wfb, trapped_wfb, &
         il_is_wfb, is_ttp, il_is_passing, il_is_trapped, &
         is_upper_bounce_point, is_lower_bounce_point, ng2
    use optionals, only: get_option_with_default
    use redistribute, only: fill
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(out) :: g_h
    logical, intent(in), optional :: force_passing
    complex, dimension(-ntgrid:ntgrid, 2) :: g2
    integer :: ig, iglo
    integer :: ik, il
    integer :: ilmin
    logical :: use_pass_homog_boundary, use_pass_homog
    logical :: force_include_passing, is_wfb, is_passing, is_trapped

    ! Initialise entire array to 0, for this reason we can take g_h as
    ! intent(out) as we don't use any existing values on input.
    call zero_array(g_h)

    ! Decide if we need to use the homogeneous solution due to the
    ! parallel boundary conditions.  This primarily impacts passing
    ! particles -- trapped particles will always require the
    ! homogeneous solution.
    use_pass_homog_boundary = &
         (boundary_option_switch == boundary_option_self_periodic) .or. &
         (boundary_option_switch == boundary_option_linked)

    ! Decide if we want to force the inclusion of the passing contribution
    force_include_passing = get_option_with_default(force_passing, .false.)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, il, g2, use_pass_homog, ig, &
    !$OMP ilmin, is_wfb, is_passing, is_trapped) &
    !$OMP SHARED(g_lo, g_h, use_pass_homog_boundary, is_ttp, mixed_wfb, &
    !$OMP aky, wfb, passing_wfb, trapped_wfb, ntgrid, r, ng2, force_include_passing) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       !/////////////////////////
       !// Constants, flags etc.
       !/////////////////////////

       ik = ik_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)

       use_pass_homog = use_pass_homog_boundary .or. aky(ik) == 0.0 .or. force_include_passing

       if (use_pass_homog) then
          ilmin = 1
       else
          ilmin = ng2 + 1
       end if

       ! If il < il min then g_h has 0 boundary/initial condition and hence should be
       ! exactly zero everywhere, and therefore there is no further work to do so cycle
       ! We might still want to consider calculating the homogeneous solution in this
       ! case "for interest". To do this we would remove this cycle and always set the boundary
       ! values for passing particles.
       if (il < ilmin) cycle

       ! ng2+1 is WFB.
       is_wfb = il_is_wfb(il)

       ! Is this particle passing. Note we exclude the wfb here as we can choose
       ! how we treat this as given by passing_wfb, trapped_wfb etc.
       is_passing = il_is_passing(il)

       ! Is this particle trapped? Note we exclude the wfb here as we can choose
       ! how we treat this as given by passing_wfb, trapped_wfb etc.
       is_trapped = il_is_trapped(il)

       !/////////////////////////
       !// Initialisation
       !/////////////////////////

       !CMR
       ! g2 simply stores trapped homogeneous boundary conditions at bounce points
       !     g2(iub,2) = 1.0        iub is UPPER bounce point, trapped bc for vpar < 0
       !     g2(ilb,1) = g1(ilb,2)  ilb is LOWER bounce point, trapped bc for vpar > 0
       ! otherwise g2 is zero
       !      -- g2 = 0 for all passing particles
       !      -- g2 = 0 for wfb as forbid always false
       !      -- g2 = 0 for ttp too
       ! g2 NOT used for totally trapped particles at il=nlambda
       !
       ! NB with trapped particles lmax=nlambda-1, and ttp is at il=nlambda
       g2 = 0.0

       !/////////////////////////
       !// Boundaries
       !/////////////////////////

       if (use_pass_homog .and. is_passing) then
          g_h(-ntgrid, 1, iglo) = 1.0
          g_h( ntgrid, 2, iglo) = 1.0
       end if

       if (is_wfb .and. mixed_wfb ) then
          ! wfb should be unity here; variable is for testing
          g_h(-ntgrid, 1, iglo) = wfb
          g_h( ntgrid, 2, iglo) = wfb
       else if (is_wfb .and. trapped_wfb) then
          ! MRH Only set the upper bounce point initial condition
          ! Should this be wfb for consistency with mixed treatment above?
          g_h( ntgrid, 2, iglo) = 1.0
       else if (is_wfb .and. passing_wfb) then
          ! MRH give the passing initial condition
          ! Should this be wfb for consistency with mixed treatment above?
          g_h(-ntgrid, 1, iglo) = 1.0
          g_h( ntgrid, 2, iglo) = 1.0
       endif

       ! Set g2=1 at UPPER bounce point for trapped (not ttp or wfb) at vpar<0
       if (is_trapped) then
          do ig=-ntgrid, ntgrid
             ! Skip for ttp
             if (is_ttp(ig, il)) cycle
             if (is_upper_bounce_point(ig, il)) g2(ig,2) = 1.0
          end do
       end if

       !/////////////////////////
       !// Time advance
       !/////////////////////////

       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! time advance vpar < 0   !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! time advance vpar < 0 homogeneous part: g_h
       do ig = ntgrid-1, -ntgrid, -1
          g_h(ig, 2, iglo) = -g_h(ig+1, 2, iglo) * r(ig,2,iglo) + g2(ig,2)
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! time advance vpar > 0   !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! First set BCs for trapped particles at lower bounce point

       if (is_trapped) then
          ! match boundary conditions at lower bounce point
          do ig = -ntgrid, ntgrid-1
             ! Skip for ttp
             if (is_ttp(ig, il)) cycle
             ! Note could change ig loop to full theta grid now we are using
             ! is_lower_bounce_point
             if (is_lower_bounce_point(ig+1, il)) then
                !CMR, 17/4/2012: set g2=(ig+1,1) = g_h(ig+1,2, iglo) where ig+1
                !     is LOWER bounce point (previously g2 was set to 1 at ig
                !      just LEFT of lower bounce point but this was handled
                !      consistently in integration of g1)
                !
                g2(ig+1,1) = g_h(ig+1, 2, iglo)
             end if
          end do
       end if

       ! If trapped wfb enforce the bounce condition at the lower bounce point
       ! Note this treats the wfb as bouncing at the end of the parallel domain.
       ! For nperiod > 1 this means the wfb does not bounce at the interior points
       ! where bmag = bmax
       if (is_wfb .and. trapped_wfb) then
          ig = -ntgrid
          g_h(ig, 1, iglo) = g_h(ig, 2, iglo)
       endif

       ! time advance vpar > 0 homogeneous part
       do ig = -ntgrid, ntgrid-1
          ! Skip for ttp
          if (is_ttp(ig, il)) cycle

          !CMR, 17/4/2012:  use consistent homogeneous trapped solution (g2) at lbp
          g_h(ig+1, 1, iglo) = -g_h(ig, 1, iglo) * r(ig, 1, iglo) + g2(ig+1, 1)
       end do
    end do
    !$OMP END PARALLEL DO

    ! Communicate the homogeneous solution boundary values for use in the
    ! twist-shift boundary conditions if required.
    if (allocated(g_adj) .and. .not. no_connections) then
       call fill (links_h, g_h, g_adj)
       if (mixed_wfb) call fill(wfb_h, g_h, g_adj)
    end if

  end subroutine init_homogeneous_g

  !> FIXME : Add documentation  
  subroutine init_bc
    use kt_grids, only: naky, ntheta0
    implicit none

    if (.not. allocated(l_links)) then
       allocate (l_links(ntheta0, naky))
       allocate (r_links(ntheta0, naky))
       allocate (n_links(2, ntheta0, naky))
    end if
    l_links = 0. ; r_links = 0. ; n_links = 0.

    select case (boundary_option_switch)
    case (boundary_option_linked)
       call init_connected_bc
       if(def_parity)call init_enforce_parity(parity_redist)
    case default
       i_class = 1
       if (.not. allocated(M_class)) then
          allocate (M_class(i_class))
          allocate (N_class(i_class))
       end if
       M_class = naky*ntheta0 ; N_class = 1
    end select

    allocate(supercell_labels(ntheta0, naky))
    allocate(n_supercells(naky))
    allocate(n_cells(ntheta0, naky))
    call calculate_supercell_labels(supercell_labels, n_supercells, n_cells)
  end subroutine init_bc

  !> FIXME : Add documentation  
  subroutine init_source_term
    use dist_fn_arrays, only: vpac, aj0
    use run_parameters, only: has_apar
    use gs2_time, only: wunits, code_dt
    use species, only: spec,nspec,nonmaxw_corr
    use hyper, only: D_res
    use theta_grid, only: ntgrid, itor_over_b
    use le_grids, only: negrid, energy
    use constants, only: zi,pi
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
#ifdef LOWFLOW
    use lowflow, only: vparterm, wdfac, wstarfac
#endif
    implicit none
    integer :: iglo
    integer :: ig, ik, it, il, ie, is, isgn

    !Initialise ufac for use in set_source
    if (.not. allocated(ufac)) then
       allocate (ufac(negrid, nspec))
       do ie = 1, negrid
          do is = 1, nspec
             ufac(ie, is) = (2.0*spec(is)%uprim &
                  + spec(is)%uprim2*energy(ie,is)**(1.5)*sqrt(pi)/4.0)
          end do
       end do
    endif

    if(.not.opt_source) return

    !Setup the source coefficients
    !See comments in get_source_term and set_source
    !for information on these terms.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ie, il, ik, it, is, isgn, ig) &
    !$OMP SHARED(g_lo, ntgrid, source_coeffs, vpar, &
#ifdef LOWFLOW
    !$OMP wdfac, wstarfac, vparterm, &
#endif
    !$OMP vpac, code_dt, wunits, ufac, omprimfac, g_exb, itor_over_B, &
    !$OMP spec, has_apar, aj0, D_res, nonmaxw_corr, wstar, wdrift) &
    !$OMP SCHEDULE(static)
    do iglo=g_lo%llim_proc,g_lo%ulim_proc
       ie=ie_idx(g_lo,iglo)
       il=il_idx(g_lo,iglo)
       ik=ik_idx(g_lo,iglo)
       it=it_idx(g_lo,iglo)
       is=is_idx(g_lo,iglo)
       do isgn=1,2
          do ig=-ntgrid,ntgrid-1
#ifdef LOWFLOW
             !Phi m
             source_coeffs(1,ig,isgn,iglo)=-2*vparterm(ig,isgn,iglo)
             !Phi p
             source_coeffs(2,ig,isgn,iglo)=-zi*wdfac(ig,isgn,iglo) &
              + zi*(wstarfac(ig,isgn,iglo) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm)
             if(has_apar)then
                !Apar m
                source_coeffs(3,ig,isgn,iglo)=-spec(is)%zstm*&
                     vpac(ig,isgn,iglo)*(aj0(ig+1,iglo)+aj0(ig,iglo))*0.5
                !Apar p
                source_coeffs(4,ig,isgn,iglo)=D_res(it,ik)*spec(is)%zstm*&
                     vpac(ig,isgn,iglo)-spec(is)%stm*vpac(ig,isgn,iglo)*&
                     zi*(wstarfac(ig,isgn,iglo) &
                     + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
                     -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) 
             endif
#else

             !Phi m
             source_coeffs(1,ig,isgn,iglo)=-2*vpar(ig,isgn,iglo)*nonmaxw_corr(ie,is)
             !Phi p
             source_coeffs(2,ig,isgn,iglo)=-zi*wdrift(ig,isgn,iglo)*nonmaxw_corr(ie,is) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm)
             if(has_apar)then
                !Apar m
                source_coeffs(3,ig,isgn,iglo)=-spec(is)%zstm*&
                     vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is)*(aj0(ig+1,iglo)+aj0(ig,iglo))*0.5
                !Apar p
                source_coeffs(4,ig,isgn,iglo)=D_res(it,ik)*spec(is)%zstm*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is)&
                     -spec(is)%stm*vpac(ig,isgn,iglo)*&
                     zi*(wstar(ik,ie,is) &
                     + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
                     -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) 
             endif
#endif
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine init_source_term

  !> Extracts the field like data corresponding to the supercell with
  !> the {it, ik} wavenumber pair as a member and constructs the
  !> ballooning space form. By default includes duplicate theta points
  !> but can drop these on request.
  subroutine flux_tube_field_like_to_ballooning_space(it, ik, field_like_ft, &
       theta_bs, field_like_bs, drop_duplicates)
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: theta0, ntheta0
    use sorting, only: quicksort
    use optionals, only: get_option_with_default
    implicit none
    integer, intent(in) :: it, ik
    complex, dimension(-ntgrid:, :, :), intent(in) :: field_like_ft
    real, dimension(:), allocatable, intent(out) :: theta_bs
    complex, dimension(:), allocatable, intent(out) :: field_like_bs
    logical, intent(in), optional :: drop_duplicates
    real, dimension(:), allocatable :: theta0_values
    real, dimension(:), allocatable :: theta0_indices_real
    integer, dimension(:), allocatable :: theta0_indices
    logical, dimension(:), allocatable :: supercell_mask
    real, dimension(:, :), allocatable :: theta_2d
    complex, dimension(:, :), allocatable :: extracted_field
    integer :: isc, ncell, it0, icell
    integer :: nextend, iend, istart, ntheta_minus_1
    logical :: drop_duplicate_theta_points

    ! Get the supercell label related to the point of interest and
    ! form a logical mask selecting members of this supercell
    isc = supercell_labels(it, ik)
    allocate(supercell_mask(ntheta0))
    supercell_mask = supercell_labels(:, ik) == isc

    ! Count how many cells there are.
    ncell = count(supercell_mask)

    theta0_values = pack(theta0(:, ik), supercell_mask)
    theta0_indices_real = pack([(it0*1.0, it0 = 1, ntheta0)], supercell_mask)

    ! Here we're sorting theta0 indices for member of the supercell
    ! with the theta0 values as a key so that theta0_indices_real will give
    ! the it values of the supercell members in order of increasing theta0.
    ! We could skip this sorting and just work out where we transition from the
    ! +ve kx to the -ve kx and slice there, but sorting should work in general
    ! and avoids potentially confusing integer index calculations
    call quicksort(ncell, theta0_values, theta0_indices_real)

    ! Convert theta0_indices_real to integer to allow use for indexing
    ! arrays later.
    allocate(theta0_indices(size(theta0_indices_real)))
    theta0_indices = int(theta0_indices_real)
    theta0_indices = theta0_indices(ncell:1:-1)

    allocate(extracted_field(-ntgrid:ntgrid, ncell))
    allocate(theta_2d(-ntgrid:ntgrid, ncell))

    do icell = 1, ncell
       extracted_field(:, icell) = field_like_ft(:, theta0_indices(icell), ik)
       theta_2d(:, icell) = theta - theta0(theta0_indices(icell), ik)
    end do

    drop_duplicate_theta_points = get_option_with_default(drop_duplicates, .false.)
    if (drop_duplicate_theta_points) then
       ! Get the size of the domain, dropping one theta point for all but the last cell
       nextend = size(theta_2d) - (ncell - 1)
       allocate(theta_bs(nextend))
       allocate(field_like_bs(nextend))

       ntheta_minus_1 = size(theta) - 1
       istart = 1
       iend = istart + ntheta_minus_1 - 1

       do icell = 1, ncell - 1
          theta_bs(istart:iend) = theta_2d(:ntgrid-1, icell)
          field_like_bs(istart:iend) = extracted_field(:ntgrid-1, icell)
          istart = iend + 1
          iend = istart + ntheta_minus_1 - 1
       end do
       theta_bs(istart:) = theta_2d(:, ncell)
       field_like_bs(istart:) = extracted_field(:, ncell)
    else
       theta_bs = reshape(theta_2d, [size(theta_2d)])
       field_like_bs = reshape(extracted_field, [size(extracted_field)])
    end if

  end subroutine flux_tube_field_like_to_ballooning_space

  !> Assigns a label to each {it, ik} point denoting which unique
  !> supercell the point belongs to.
  subroutine calculate_supercell_labels(supercell_labels, n_supercells, n_cells)
    use kt_grids, only: naky, ntheta0
    implicit none
    integer, dimension(:, :), intent(out) :: supercell_labels, n_cells
    integer, dimension(:), intent(out) :: n_supercells
    integer, dimension(:), allocatable :: temporary_it_leftmost
    integer :: it, ik, isc, ncell

    supercell_labels = -1
    n_supercells = 0
    n_cells = 0
    isc = -1

    allocate(temporary_it_leftmost(ntheta0))
    temporary_it_leftmost = -1

    do ik = 1, naky
       do it = 1, ntheta0
          temporary_it_leftmost(it) = get_leftmost_it(it, ik)
       end do

       do it = 1, ntheta0
          if (temporary_it_leftmost(it) >= 0) then
             ! Note ncell isn't currently used but records the number
             ! of cells for the current supercell. We struggle to
             ! store it efficiently at this point as we don't yet know
             ! how many supercells there are. We currently store it in
             ! an array of shape [ntheta0, naky] so that we can look
             ! up the number of cells given {it, ik} rather than
             ! {isc}, so there is some duplication.
             ncell = count(temporary_it_leftmost == temporary_it_leftmost(it))
             isc = isc + 1
             where (temporary_it_leftmost == temporary_it_leftmost(it))
                temporary_it_leftmost = -1
                supercell_labels(:, ik) = isc
                n_cells(:, ik) = ncell
             end where
             n_supercells(ik) = n_supercells(ik) + 1
          end if
       end do
    end do

  end subroutine calculate_supercell_labels

  !> Compute the spacing in kx indices between connected kx domains
  !> (cells) for the lowest non-zonal ky.
  integer function compute_jshift0() result(jshift0)
    use theta_grid, only: ntgrid, theta
    use kt_grids, only: naky, ntheta0, theta0, aky
    use mp, only: mp_abort
    implicit none
    real :: theta_extent, theta0_spacing

    ! jshift0 corresponds to J (not delta j) from Beer thesis (unsure about +0.01) -- MAB
    ! delta j is the number of akxs (i.e. theta0s) 'skipped' when connecting one
    ! parallel segment to the next to satisfy the twist and shift
    ! boundary conditions. delta j = (ik - 1) * jshift0 . EGH

    theta_extent = theta(ntgrid) - theta(-ntgrid)
    if (naky > 1 .and. ntheta0 > 1) then
       theta0_spacing = theta0(2, 2) - theta0(1, 2)
    else if (naky == 1 .and. ntheta0 > 1 .and. aky(1) /= 0.0) then
       ! Not clear we _ever_ get here as box mode usually requires
       ! aky(1) == 0. One could imagine potentially using linked boundary
       ! conditions and kt_grids_range to end up here.
       theta0_spacing = theta0(2, 1) - theta0(1, 1)
    else
       ! Not possible to work out a theta0 spacing here
       ! so just set up such that jshift0 == 1
       theta0_spacing = theta_extent
    end if

    ! Work out how many theta0 "grid points" we must step over
    ! to travel a distance equiavlent from one end of a cell
    ! to the other end.
    jshift0 = nint(theta_extent/theta0_spacing)

    ! Sanity check the setup.
    if (abs(jshift0 * theta0_spacing - theta_extent) > 100*epsilon(0.0)) then
       call mp_abort('Calculation of jshift0 indicates inconsistent grids.', .true.)
    end if
  end function compute_jshift0

  !> For the passed jshift0 value determine the it (kx) indices to
  !> the left and right of each {it, ik} index pair. Used to construct
  !> the connected boundary conditions.
  subroutine compute_itleft_and_itright(jshift0, itleft, itright)
    use kt_grids, only: naky, ntheta0, aky
    implicit none
    integer, intent(in) :: jshift0
    integer, dimension(:, :), intent(out) :: itleft, itright
    integer :: ik, it, it0, itl, itr

    ! Note the below appears to assume ik = 1 is always the ky = 0
    ! mode but this is not consistent with potentially more general
    ! treatement in the init_connected_bc related code. We also revisit
    ! ik = 1 in the below loop so this initialisation is likely
    ! overwritten in all cases.
    itleft(:, 1) = -1 ! No connected segments when ky = 0. EGH
    itright(:, 1) = -1

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ik, it, it0, itl, itr) &
    !$OMP SHARED(naky, ntheta0, aky, jshift0, itleft, itright) &
    !$OMP COLLAPSE(2) &
    !$OMP SCHEDULE(static)
    do ik = 1, naky
       do it = 1, ntheta0
          if (it > (ntheta0+1)/2) then
             ! akx is negative (-1 shift because akx(1) = 0) -- MAB
             it0 = it - ntheta0 - 1
          else
             ! akx is positive (-1 shift because akx(1) = 0) -- MAB
             it0 = it - 1
          end if

          if (ik == 1) then
             if (aky(ik) /= 0.0 .and. naky == 1) then
                ! for this case, jshift0 is delta j from Beer thesis -- MAB
                itl = it0 + jshift0
                itr = it0 - jshift0
             else
                itl = it0
                itr = it0
             end if
          else
             ! ik = 1 corresponds to aky=0, so shift index by 1
             ! (ik-1)*jshift0 is delta j from Beer thesis -- MAB
             itl = it0 + (ik-1)*jshift0
             itr = it0 - (ik-1)*jshift0
          end if

          ! remap to usual GS2 indices -- MAB
          if (itl >= 0 .and. itl < (ntheta0+1)/2) then
             itleft(it, ik) = itl + 1
          else if (itl + ntheta0 + 1 > (ntheta0+1)/2 &
               .and. itl + ntheta0 + 1 <= ntheta0) then
             itleft(it, ik) = itl + ntheta0 + 1
          else
             ! for periodic, need to change below -- MAB/BD
             ! j' = j + delta j not included in simulation, so can't connect -- MAB
             itleft(it, ik) = -1
          end if

          ! same stuff for j' = j - delta j -- MAB
          if (itr >= 0 .and. itr < (ntheta0+1)/2) then
             itright(it, ik) = itr + 1
          else if (itr + ntheta0 + 1 > (ntheta0+1)/2 &
               .and. itr + ntheta0 + 1 <= ntheta0) then
             itright(it, ik) = itr + ntheta0 + 1
          else
             ! for periodic, need to change below -- MAB/BD
             itright(it, ik) = -1
          end if
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine compute_itleft_and_itright

  !> Look up and store the iglo index and responsible processor for connections
  !> to the left and right of each local iglo index. Note this is only interested
  !> in passing particles and the non-trapped wfb. Trapped particles are
  !> considered to have no connections.
  subroutine compute_connections(itleft, itright, connections)
    use gs2_layouts, only: g_lo, il_idx, ik_idx, it_idx, ie_idx, is_idx, proc_id, idx
    use le_grids, only: il_is_trapped, il_is_wfb, trapped_wfb
    implicit none
    integer, dimension(:, :), intent(in) :: itleft, itright
    type(connections_type), dimension(g_lo%llim_proc:), intent(out) :: connections
    integer :: iglo, il, ik, it, ie, is, iglo_connection

    ! Note connections_type have initial values reflecting no connections
    ! so we don't need to explicitly populate every instance.

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ik, it, ie, is, iglo_connection) &
    !$OMP SHARED(g_lo, itleft, itright, connections, trapped_wfb) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ! get processors and indices for j' (kx') modes connecting
       ! to mode j (kx) so we can set up communication -- MAB

       il = il_idx(g_lo,iglo)

       ! if non-wfb trapped particle, no connections
       ! or if wfb is trapped, no connections
       if (il_is_trapped(il) .or. (il_is_wfb(il) .and. trapped_wfb)) cycle

       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)

       ! If no links then cycle
       if (itleft(it, ik) < 0 .and. itright(it, ik) < 0) cycle

       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       if (itleft(it, ik) >= 0) then
          iglo_connection = idx(g_lo, ik, itleft(it, ik), il, ie, is)
          connections(iglo)%iproc_left = proc_id(g_lo, iglo_connection)
          connections(iglo)%iglo_left = iglo_connection
          connections(iglo)%neighbor = .true.
       end if

       if (itright(it, ik) >= 0) then
          iglo_connection = idx(g_lo, ik, itright(it, ik), il, ie, is)
          connections(iglo)%iproc_right = proc_id(g_lo, iglo_connection)
          connections(iglo)%iglo_right = iglo_connection
          connections(iglo)%neighbor = .true.
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine compute_connections

  !> Set the save_h flag based on connections data
  subroutine set_save_h(connections, save_h)
    use gs2_layouts, only: g_lo, il_idx, ik_idx
    use kt_grids, only: aky
    use le_grids, only: il_is_wfb, passing_wfb
    implicit none
    type(connections_type), dimension(g_lo%llim_proc:), intent(in) :: connections
    logical, dimension(:, g_lo%llim_proc:), intent(out) :: save_h
    integer :: iglo, ik, il

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, il) &
    !$OMP SHARED(g_lo, connections, aky, save_h, passing_wfb) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       save_h (1, iglo) = connections(iglo)%iglo_left >= 0 .and. aky(ik) /= 0.0
       save_h (2, iglo) = connections(iglo)%iglo_right >= 0 .and. aky(ik) /= 0.0

       il = il_idx(g_lo,iglo)
       if (il_is_wfb(il)) then
          ! wfb (linked)
          if ( connections(iglo)%neighbor .and. &
               (.not. passing_wfb) .and. aky(ik) /= 0.0) then
             ! MRH below should not occur if we are treating wfb as standard passing
             ! since save_h for passing particles is handled above
             ! If trapped_wfb then neighbor will be false, here we are
             ! dealing with the original "mixed" wfb treatment where
             ! trapped_wfb = passing_wfb = .false.
             ! One might expect neighbor to be false if aky == 0
             ! as the zonal mode is periodic rather than linked.
             save_h (:,iglo) = .true.
          end if
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine set_save_h

  !> Count how many links are to the left/right of each {it,ik} point
  subroutine count_links_one_side(it_connections, links)
    use kt_grids, only: ntheta0, naky
    use mp, only: mp_abort
    implicit none
    integer, dimension(:, :), intent(in) :: it_connections
    !> Links is the number of links to the left or right
    integer, dimension(:, :), intent(out) :: links
    integer :: it, ik, it_star, it_star_next

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ik, it, it_star, it_star_next) &
    !$OMP SHARED(naky, ntheta0, links, it_connections) &
    !$OMP COLLAPSE(2) &
    !$OMP SCHEDULE(static)
    do ik = 1, naky
       do it = 1, ntheta0
          links(it, ik) = 0
          it_star = it
          loop_over_connections: do
             it_star_next = it_connections(it_star, ik)
             ! If connected to self (periodic) exit
             if (it_star == it_star_next) exit loop_over_connections
             ! If no more connections exit
             if (it_star_next < 0) exit loop_over_connections
             links(it, ik) = links(it, ik) + 1
             it_star = it_star_next
             if (links(it, ik) > max_allowed_links) then
                call mp_abort('links error',.true.)
             endif
          end do loop_over_connections
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine count_links_one_side

  !> Count the links to the left and right of each {it, ik} point
  !> and work out how many values are required to set the bc.
  subroutine count_links(itleft, itright, l_links, r_links, n_links)
    implicit none
    integer, dimension(:, :), intent(in) :: itleft, itright
    integer, dimension(:, :), intent(out) :: l_links, r_links
    integer, dimension(:, :, :), intent(out) :: n_links

    ! First count the number of links to the left and right
    call count_links_one_side(itleft, l_links)
    call count_links_one_side(itright, r_links)

    ! Then work out the number of values needed for the bc.
    ! 'n_links' complex numbers are needed to specify bc for (it, ik) region
    ! ignoring wfb
    ! n_links(1,:,:) is for v_par > 0, etc.

    where (l_links == 0)
       n_links(1, :, :) = 0
    elsewhere
       n_links(1, :, :) = 2 * l_links - 1
    end where

    where (r_links == 0)
       n_links(2, :, :) = 0
    elsewhere
       n_links(2, :, :) = 2 * r_links - 1
    end where

  end subroutine count_links

  !> Count the number of unique supercell sizes and setup the fields
  !> implicit arrays M_class and N_class.
  subroutine setup_fields_implicit_class_arrays(l_links, r_links, M_class, N_class, i_class)
    use kt_grids, only: naky, ntheta0
    use mp, only: mp_abort, iproc
    use sorting, only: quicksort
    implicit none
    integer, dimension(:, :), intent(in) :: l_links, r_links
    !> N_class(i) = number of linked cells for i_th class (Size of supercell)
    integer, dimension(:), allocatable, intent(out) :: N_class
    !> M_class(i) = number of members in i_th class (How many supercells of this size)
    integer, dimension(:), allocatable, intent(out) :: M_class
    !> i_class = number of classes (unique sizes of supercell)
    integer, intent(out) :: i_class
    integer, dimension(naky*ntheta0) :: n_k, N_class_full
    integer :: k, it, ik, i, j
    ! Now set up class arrays for the implicit fields

    ! First count number of linked cells for each (kx, ky)
    k = 1
    do ik = 1, naky
       do it = 1, ntheta0
          n_k(k) = 1 + l_links(it, ik) + r_links(it, ik)
          k = k + 1
       end do
    end do

    ! Count how many unique values of n_k there are.  This is the number
    ! of classes.

    ! Sort:
    call quicksort(naky*ntheta0, n_k)

    ! Then count how many unique values there are, storing
    ! each unique value as it is encountered.
    i_class = 1
    N_class_full(i_class) = n_k(1)
    do k = 2, naky*ntheta0
       if (n_k(k) == n_k(k-1)) cycle
       i_class = i_class + 1
       N_class_full(i_class) = n_k(k)
    end do

    ! Extract the unique sizes into N_class
    N_class = N_class_full(1:i_class)

    ! Allocate M
    allocate (M_class(i_class))

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(j) &
    !$OMP SHARED(i_class, M_class, n_k, N_class) &
    !$OMP SCHEDULE(static)
    do j = 1, i_class
       ! Note this assumes N_class(j) is an exact factor of count(n_k == N_class(j))
       ! but we have a consistency check for this later.
       M_class(j) = count(n_k == N_class(j)) / N_class(j)
    end do
    !$OMP END PARALLEL DO

    ! Check for consistency:

    ! j is number of linked cells in class structure
    j = 0
    do i = 1, i_class
       j = j + N_class(i)*M_class(i)
       ! Check that we have N_class which is a factor of the matching n_k
       if ( abs(M_class(i) - count(n_k == N_class(i)) / (1.0 * N_class(i))) &
            > 10 * epsilon(0.0)) then
          write(*,*) 'PE ',iproc,'number of instances of supercell size not ' // &
               'an integer times the size of supercell: M,N,n_N',&
               M_class(i), N_class(i), count(n_k == N_class(i))
          call mp_abort('Problem with connected BC fields implicit arrays.')
       end if
    end do

    if (j /= naky*ntheta0) then
       write(*,*) 'PE ',iproc,'has j= ',j,' k= ',naky*ntheta0,' : Stopping'
       call mp_abort('problem with connnected bc')
    end if
  end subroutine setup_fields_implicit_class_arrays

  !> Populates the redistribute types describing the communication pattern
  !> required to deal with linked boundary conditions. Only deals with
  !> passing particles, including the wfb _if_ passing_wfb is true.
  subroutine setup_connected_bc_redistribute(l_links, r_links, n_links_max, &
       links_p, links_h, no_connections)
    use gs2_layouts, only: g_lo, il_idx, ik_idx, it_idx, proc_id, ie_idx, is_idx
    use le_grids, only: il_is_passing, il_is_wfb, passing_wfb
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use mp, only: nproc, max_allreduce, iproc
    use theta_grid, only: ntgrid
    implicit none
    integer, dimension(:, :), intent(in) :: l_links, r_links
    integer, intent(in) :: n_links_max
    type(redist_type), intent(out) :: links_p, links_h
    logical, intent(out) :: no_connections
    type (index_list_type), dimension(0:nproc-1) :: from_p, from_h, to_p, to_h
    integer, dimension (0:nproc-1) :: nn_from, nn_to, nn_from_h, nn_to_h
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: il, ik, it, ncell, ip, iglo_star, j, n, iglo, n_h, nn_max
    integer :: iglo_right, ipright, iglo_left, ipleft, ie, is

    no_connections = .false.
    !<DD>Note the communications setup here are often equivalent to an all-to-all type
    !communication, i.e. when nproc --> 2 nproc, t_fill --> 4 t_fill
    !See comments in invert_rhs_linked for more details.

    !Note: This setup currently involves several loops over the entire domain
    !and hence does not scale well (effectively serial code). This can come to
    !dominate initialisation at large core count.

    nn_to = 0
    nn_from = 0
    nn_to_h = 0
    nn_from_h = 0

    ! Whilst we loop over the entire domain here, we should not have any work
    ! to do for ik, ie, is and il indices which do not belong to this processor's
    ! local g_lo range. We may be better off writing this as an explicit loop over
    ! the xyles dimensions with yles range being set my the min/max that this
    ! processor sees. For now we can just get the indices and cycle.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ik, ie, is, it, ip, ncell, iglo_star, iglo_right, &
    !$OMP ipright, iglo_left, ipleft, j) &
    !$OMP SHARED(g_lo, iproc, r_links, l_links, passing_wfb) &
    !$OMP REDUCTION(+: nn_to, nn_from, nn_to_h, nn_from_h) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_world, g_lo%ulim_world
       !CMR, 20/10/2013:
       !     Communicate pattern sends passing particles to downwind linked cells
       !      (ntgrid,1,iglo)  => each RH connected cell (j,1,iglo_right)
       !      (-ntgrid,2,iglo) => each LH connected cell (j,2,iglo_left)
       !          where j index in receive buffer = #hops in connection
       !                         (nb j=1 is nearest connection)

       il = il_idx(g_lo,iglo)
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ! Exclude trapped particles (inc wfb if it is not passing)
       if (.not. (il_is_passing(il) .or. (il_is_wfb(il) .and. passing_wfb))) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_star = iglo
       do j = 1, r_links(it, ik)
          call get_right_connection (iglo_star, iglo_right, ipright)
          ! sender
          if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
          ! receiver
          if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

          iglo_star = iglo_right

          !Special counting for links_h
          if(l_links(it, ik)>0) then
             if (ip == iproc) nn_from_h(ipright) = nn_from_h(ipright) + 1
             if (ipright == iproc) nn_to_h(ip) = nn_to_h(ip) + 1
          endif
       end do

       iglo_star = iglo
       do j = 1, l_links(it, ik)
          call get_left_connection (iglo_star, iglo_left, ipleft)
          ! sender
          if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
          ! receiver
          if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1

          iglo_star = iglo_left

          !Special counting for links_h
          if(r_links(it, ik)>0) then
             if (ip == iproc) nn_from_h(ipleft) = nn_from_h(ipleft) + 1
             if (ipleft == iproc) nn_to_h(ip) = nn_to_h(ip) + 1
          endif
       end do
    end do
    !$OMP END PARALLEL DO

    ! Check if we have any communication to to do or not. Note here communication
    ! includes local copying -- this is really checking do we need to fix up the
    ! parallel boundaries anywhere or not.
    nn_max = maxval(nn_to)
    call max_allreduce (nn_max)
    if (nn_max == 0) then
       no_connections = .true.
       ! Note we don't set links_p and links_h in this case
       return
    end if

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from_p(ip)%first(nn_from(ip)))
          allocate (from_p(ip)%second(nn_from(ip)))
          allocate (from_p(ip)%third(nn_from(ip)))
       endif
       if (nn_from_h(ip) > 0) then
          allocate (from_h(ip)%first(nn_from_h(ip)))
          allocate (from_h(ip)%second(nn_from_h(ip)))
          allocate (from_h(ip)%third(nn_from_h(ip)))
       endif
       if (nn_to(ip) > 0) then
          allocate (to_p(ip)%first(nn_to(ip)))
          allocate (to_p(ip)%second(nn_to(ip)))
          allocate (to_p(ip)%third(nn_to(ip)))
       endif
       if (nn_to_h(ip)>0) then
          allocate (to_h(ip)%first(nn_to_h(ip)))
          allocate (to_h(ip)%second(nn_to_h(ip)))
          allocate (to_h(ip)%third(nn_to_h(ip)))
       endif
    end do

    nn_from = 0
    nn_from_h=0
    nn_to = 0
    nn_to_h=0

    do iglo = g_lo%llim_world, g_lo%ulim_world

       il = il_idx(g_lo,iglo)
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ! Exclude trapped particles (inc wfb if it is not passing)
       if (.not. (il_is_passing(il) .or. (il_is_wfb(il) .and. passing_wfb))) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_star = iglo
       do j = 1, r_links(it, ik)
          call get_right_connection (iglo_star, iglo_right, ipright)
          ! sender
          if (ip == iproc) then
             n = nn_from(ipright) + 1
             nn_from(ipright) = n
             from_p(ipright)%first(n) = ntgrid
             from_p(ipright)%second(n) = 1
             from_p(ipright)%third(n) = iglo
             !Special restriction for links_h
             if(l_links(it, ik)>0)then
                n_h=nn_from_h(ipright)+1
                nn_from_h(ipright)=n_h
                from_h(ipright)%first(n_h) = ntgrid
                from_h(ipright)%second(n_h) = 1
                from_h(ipright)%third(n_h) = iglo
             endif
          end if
          ! receiver
          if (ipright == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_p(ip)%first(n) = j
             to_p(ip)%second(n) = 1
             to_p(ip)%third(n) = iglo_right
             !Special restriction for links_h
             if(l_links(it, ik)>0)then
                n_h=nn_to_h(ip)+1
                nn_to_h(ip)=n_h
                to_h(ip)%first(n_h) = 2*l_links(it, ik)+j
                to_h(ip)%second(n_h) = 1
                to_h(ip)%third(n_h) = iglo_right
             endif
          end if
          iglo_star = iglo_right
       end do

       iglo_star = iglo
       do j = 1, l_links(it, ik)
          call get_left_connection (iglo_star, iglo_left, ipleft)
          ! sender
          if (ip == iproc) then
             n = nn_from(ipleft) + 1
             nn_from(ipleft) = n
             from_p(ipleft)%first(n) = -ntgrid
             from_p(ipleft)%second(n) = 2
             from_p(ipleft)%third(n) = iglo
             !Special restriction for links_h
             if(r_links(it, ik)>0)then
                n_h=nn_from_h(ipleft)+1
                nn_from_h(ipleft)=n_h
                from_h(ipleft)%first(n_h) = -ntgrid
                from_h(ipleft)%second(n_h) = 2
                from_h(ipleft)%third(n_h) = iglo
             endif
          end if
          ! receiver
          if (ipleft == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_p(ip)%first(n) = j
             to_p(ip)%second(n) = 2
             to_p(ip)%third(n) = iglo_left
             !Special restriction for links_h
             if(r_links(it, ik)>0)then
                n_h = nn_to_h(ip) + 1
                nn_to_h(ip) = n_h
                to_h(ip)%first(n_h) = 2*r_links(it, ik)+j
                to_h(ip)%second(n_h) = 2
                to_h(ip)%third(n_h) = iglo_left
             endif
          end if
          iglo_star = iglo_left
       end do
    end do

    from_low = [-ntgrid, 1, g_lo%llim_proc]
    from_high = [ntgrid, 2, g_lo%ulim_alloc]

    to_low = [1, 1, g_lo%llim_proc]
    to_high = [n_links_max, 2, g_lo%ulim_alloc]

    call init_fill (links_p, 'c', to_low, to_high, to_p, &
         from_low, from_high, from_p)
    call init_fill (links_h, 'c', to_low, to_high, to_h, &
         from_low, from_high, from_h)

    call delete_list (from_p)
    call delete_list (from_h)
    call delete_list (to_p)
    call delete_list (to_h)

  end subroutine setup_connected_bc_redistribute

  !> Populates the redistribute types describing the communication pattern
  !> required to deal with linked boundary conditions for the mixed wfb.
  !> This is used when passing_wfb = trapped_wfb = false.
  !> This is the default.
  subroutine setup_connected_bc_redistribute_mixed_wfb(l_links, r_links, n_links_max, &
       wfb_p, wfb_h)
    use gs2_layouts, only: g_lo, il_idx, ik_idx, it_idx, proc_id, ie_idx, is_idx
    use le_grids, only: il_is_wfb
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use theta_grid, only: ntgrid
    use mp, only: nproc, iproc
    implicit none
    integer, dimension(:, :), intent(in) :: l_links, r_links
    integer, intent(in) :: n_links_max
    type(redist_type), intent(out) ::wfb_p, wfb_h
    type(index_list_type), dimension(0:nproc-1) :: from, to_p, to_h
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: il, ik, it, ncell, ip, iglo_star, j, n, iglo
    integer :: iglo_right, ipright, iglo_left, ipleft, ie, is

    nn_to = 0
    nn_from = 0

    !NOTE: No special restriction/counting for wfb_h unlike links_h
    ! Whilst we loop over the entire domain here, we should not have any work
    ! to do for ik, ie, is and il indices which do not belong to this processor's
    ! local g_lo range. We may be better off writing this as an explicit loop over
    ! the xyles dimensions with yles range being set my the min/max that this
    ! processor sees. For now we can just get the indices and cycle.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ik, ie, is, it, ncell, ip, iglo_right, iglo_left, &
    !$OMP ipleft, j, ipright, iglo_star) &
    !$OMP SHARED(g_lo, iproc, l_links, r_links) &
    !$OMP REDUCTION(+: nn_to, nn_from) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_world, g_lo%ulim_world
       !CMR, 20/10/2013:
       !     Communicate pattern sends wfb endpoint to ALL linked cells
       !      (ntgrid,1,iglo)  => ALL connected cells (j,1,iglo_conn)
       !            where j index in receive buffer = r_links(iglo)+1
       !      (-ntgrid,2,iglo) => ALL connected cells (j,2,iglo_conn)
       !         where j index in receive buffer = l_links(iglo)+1
       !

       il = il_idx(g_lo,iglo)
       if (.not. il_is_wfb(il)) cycle
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

       ! v_par > 0:
       !CMR: introduced iglo_star to make notation below less confusing
       !
       call find_leftmost_link (iglo, iglo_star, ipright)
       do j = 1, ncell
          ! sender
          if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
          ! receiver
          if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
          call get_right_connection (iglo_star, iglo_right, ipright)
          iglo_star=iglo_right
       end do

       ! v_par < 0:
       call find_rightmost_link (iglo, iglo_star, ipleft)
       do j = 1, ncell
          ! sender
          if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
          ! receiver
          if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
          call get_left_connection (iglo_star, iglo_left, ipleft)
          iglo_star=iglo_left
       end do
    end do
    !$OMP END PARALLEL DO

    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from(ip)%first(nn_from(ip)))
          allocate (from(ip)%second(nn_from(ip)))
          allocate (from(ip)%third(nn_from(ip)))
       endif
       if (nn_to(ip) > 0) then
          allocate (to_p(ip)%first(nn_to(ip)))
          allocate (to_p(ip)%second(nn_to(ip)))
          allocate (to_p(ip)%third(nn_to(ip)))
          allocate (to_h(ip)%first(nn_to(ip)))
          allocate (to_h(ip)%second(nn_to(ip)))
          allocate (to_h(ip)%third(nn_to(ip)))
       endif
    end do

    nn_from = 0
    nn_to = 0

    do iglo = g_lo%llim_world, g_lo%ulim_world

       il = il_idx(g_lo,iglo)
       if (.not. il_is_wfb(il)) cycle
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_right = iglo ; iglo_left = iglo ; ipright = ip ; ipleft = ip

       ! v_par > 0:
       call find_leftmost_link (iglo, iglo_star, ipright)
       do j = 1, ncell
          ! sender
          if (ip == iproc) then
             n = nn_from(ipright) + 1
             nn_from(ipright) = n
             from(ipright)%first(n) = ntgrid
             from(ipright)%second(n) = 1
             from(ipright)%third(n) = iglo
          end if
          ! receiver
          if (ipright == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_p(ip)%first(n) = r_links(it, ik) + 1
             to_p(ip)%second(n) = 1
             to_p(ip)%third(n) = iglo_star
             to_h(ip)%first(n) = 2*ncell-r_links(it, ik)
             to_h(ip)%second(n) = 1
             to_h(ip)%third(n) = iglo_star
          end if
          call get_right_connection (iglo_star, iglo_right, ipright)
          iglo_star=iglo_right
       end do

       ! v_par < 0:
       call find_rightmost_link (iglo, iglo_star, ipleft)
       do j = 1, ncell
          ! sender
          if (ip == iproc) then
             n = nn_from(ipleft) + 1
             nn_from(ipleft) = n
             from(ipleft)%first(n) = -ntgrid
             from(ipleft)%second(n) = 2
             from(ipleft)%third(n) = iglo
          end if
          ! receiver
          if (ipleft == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n
             to_p(ip)%first(n) = l_links(it, ik) + 1
             to_p(ip)%second(n) = 2
             to_p(ip)%third(n) = iglo_star
             to_h(ip)%first(n) = 2*ncell-l_links(it, ik)
             to_h(ip)%second(n) = 2
             to_h(ip)%third(n) = iglo_star
          end if
          call get_left_connection (iglo_star, iglo_left, ipleft)
          iglo_star=iglo_left
       end do
    end do

    from_low = [-ntgrid, 1, g_lo%llim_proc]
    from_high = [ntgrid, 2, g_lo%ulim_alloc]

    to_low = [1, 1, g_lo%llim_proc]
    to_high = [n_links_max, 2, g_lo%ulim_alloc]

    call init_fill (wfb_p, 'c', to_low, to_high, to_p, from_low, from_high, from)
    call init_fill (wfb_h, 'c', to_low, to_high, to_h, from_low, from_high, from)

    call delete_list (from)
    call delete_list (to_p)
    call delete_list (to_h)
  end subroutine setup_connected_bc_redistribute_mixed_wfb

  !> Setup the redistribute associated with passing the incoming
  !> boundary values to the previous connected cell, storing at its
  !> connected boundary. This is required in order to allow for non-zero
  !> incoming boundary conditions on internal cells.
  subroutine setup_pass_incoming_boundary_to_connections(l_links, r_links, &
       incoming_links)
    use gs2_layouts, only: g_lo, il_idx, ik_idx, it_idx, proc_id, ie_idx, is_idx, idx
    use le_grids, only: il_is_wfb, trapped_wfb, il_is_trapped, mixed_wfb
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use theta_grid, only: ntgrid
    use mp, only: nproc, iproc
    implicit none
    integer, dimension(:, :), intent(in) :: l_links, r_links
    type(redist_type), intent(out) :: incoming_links
    type (index_list_type), dimension(0:nproc-1) :: from, to
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension (3) :: to_low, from_low, to_high, from_high
    integer :: ncell, iglo_star, j, n, n_h, nn_max
    integer :: iglo_right, ipright, iglo_left, ipleft
    integer :: iglo, il, ik, it, ie, is, ip, it_leftmost, it_rightmost

    nn_to = 0
    nn_from = 0

    ! Whilst we loop over the entire domain here, we should not have any work
    ! to do for ik, ie, is and il indices which do not belong to this processor's
    ! local g_lo range. We may be better off writing this as an explicit loop over
    ! the xyles dimensions with yles range being set my the min/max that this
    ! processor sees. For now we can just get the indices and cycle.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, ik, ie, is, it, ip, ncell, iglo_star, iglo_right, &
    !$OMP ipright, iglo_left, ipleft, j, it_leftmost, it_rightmost) &
    !$OMP SHARED(g_lo, iproc, r_links, l_links, trapped_wfb, mixed_wfb) &
    !$OMP REDUCTION(+: nn_to, nn_from) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ! Communicate pattern sends incoming boundary for passing particles (and mixed_wfb)
       ! to upwind linked cells
       !      (-ntgrid,1,iglo)  => left connected cell (ntgrid,1,iglo_left)
       !      (ntgrid,2,iglo)  => right connected cell (-ntgrid,2,iglo_right)
       ! As we only need one value per {sigma, iglo} we might want to make
       ! the size of the destination array [1, 2, {glo_local}].

       il = il_idx(g_lo,iglo)
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ! Exclude trapped particles (inc wfb if it is not passing/mixed)
       if (il_is_trapped(il) .or. (il_is_wfb(il) .and. trapped_wfb)) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1

       ! Skip if only one cell -- WFB and ky==0 periodicity handled
       ! as self-periodic rather than linked in this case.
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_star = iglo

       ! Get first connection to the right.
       if (r_links(it, ik) > 0) then
          call get_right_connection (iglo_star, iglo_right, ipright)
          ! sender
          if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
          ! receiver
          if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1

          iglo_star = iglo
       end if

       ! Get first connection to the left.
       if (l_links(it, ik) > 0) then
          call get_left_connection (iglo_star, iglo_left, ipleft)
          ! sender
          if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
          ! receiver
          if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
       end if

       ! Special handling for periodicity of mixed wfb
       if (il_is_wfb(il) .and. mixed_wfb) then
          it_leftmost = get_leftmost_it(it, ik)
          it_rightmost = get_rightmost_it(it, ik)

          ! If we're at the right end we need to pass to the leftmost
          if (it == it_rightmost) then
             iglo_left = idx(g_lo, ik, it_leftmost, il, ie, is)
             ipleft = proc_id(g_lo, iglo_left)
             if (ip == iproc) nn_from(ipleft) = nn_from(ipleft) + 1
             if (ipleft == iproc) nn_to(ip) = nn_to(ip) + 1
          end if

          ! If we're at the left end we need to pass to the rightmost
          if (it == it_leftmost) then
             iglo_right = idx(g_lo, ik, it_rightmost, il, ie, is)
             ipright = proc_id(g_lo, iglo_right)
             if (ip == iproc) nn_from(ipright) = nn_from(ipright) + 1
             if (ipright == iproc) nn_to(ip) = nn_to(ip) + 1
          end if

       end if

    end do
    !$OMP END PARALLEL DO


    do ip = 0, nproc-1
       if (nn_from(ip) > 0) then
          allocate (from(ip)%first(nn_from(ip)))
          allocate (from(ip)%second(nn_from(ip)))
          allocate (from(ip)%third(nn_from(ip)))
       endif
       if (nn_to(ip) > 0) then
          allocate (to(ip)%first(nn_to(ip)))
          allocate (to(ip)%second(nn_to(ip)))
          allocate (to(ip)%third(nn_to(ip)))
       endif
    end do

    nn_from = 0
    nn_to = 0
    do iglo = g_lo%llim_world, g_lo%ulim_world
       ! Communicate pattern sends incoming boundary for passing particles (and mixed_wfb)
       ! to upwind linked cells
       !      (-ntgrid,1,iglo)  => left connected cell (ntgrid,1,iglo_left)
       !      (ntgrid,2,iglo)  => right connected cell (-ntgrid,2,iglo_right)
       ! As we only need one value per {sigma, iglo} we might want to make
       ! the size of the destination array [1, 2, {glo_local}].

       il = il_idx(g_lo,iglo)
       if (il > g_lo%il_max .or. il < g_lo%il_min) cycle

       ! Exclude trapped particles (inc wfb if it is not passing/mixed)
       if (il_is_trapped(il) .or. (il_is_wfb(il) .and. trapped_wfb)) cycle

       ik = ik_idx(g_lo,iglo)
       if (ik > g_lo%ik_max .or. ik < g_lo%ik_min) cycle

       ie = ie_idx(g_lo, iglo)
       if (ie > g_lo%ie_max .or. ie < g_lo%ie_min) cycle

       is = is_idx(g_lo, iglo)
       if (is > g_lo%is_max .or. is < g_lo%is_min) cycle

       it = it_idx(g_lo,iglo)
       ncell = r_links(it, ik) + l_links(it, ik) + 1

       ! Skip if only one cell -- WFB and ky==0 periodicity handled
       ! as self-periodic rather than linked in this case.
       if (ncell == 1) cycle

       ip = proc_id(g_lo,iglo)

       iglo_star = iglo

       ! Get first connection to the right.
       if (r_links(it, ik) > 0) then
          call get_right_connection (iglo_star, iglo_right, ipright)
          ! sender
          if (ip == iproc) then
             n = nn_from(ipright) + 1
             nn_from(ipright) = n

             from(ipright)%first(n) = ntgrid
             from(ipright)%second(n) = 2
             from(ipright)%third(n) = iglo
          end if

          ! receiver
          if (ipright == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n

             to(ip)%first(n) = -ntgrid
             to(ip)%second(n) = 2
             to(ip)%third(n) = iglo_right
          end if
       end if

       iglo_star = iglo

       ! Get first connection to the left.
       if (l_links(it, ik) > 0) then
          call get_left_connection (iglo_star, iglo_left, ipleft)
          ! sender
          if (ip == iproc) then
             n = nn_from(ipleft) + 1
             nn_from(ipleft) = n

             from(ipleft)%first(n) = -ntgrid
             from(ipleft)%second(n) = 1
             from(ipleft)%third(n) = iglo
          end if

          ! receiver
          if (ipleft == iproc) then
             n = nn_to(ip) + 1
             nn_to(ip) = n

             to(ip)%first(n) = ntgrid
             to(ip)%second(n) = 1
             to(ip)%third(n) = iglo_left
          end if
       end if

       ! Special handling for periodicity of mixed wfb
       if (il_is_wfb(il) .and. mixed_wfb) then
          it_leftmost = get_leftmost_it(it, ik)
          it_rightmost = get_rightmost_it(it, ik)

          ! If we're at the right end we need to pass to the leftmost
          if (it == it_rightmost) then
             iglo_left = idx(g_lo, ik, it_leftmost, il, ie, is)
             ipleft = proc_id(g_lo, iglo_left)
             if (ip == iproc) then
                n = nn_from(ipleft) + 1
                nn_from(ipleft) = n

                from(ipleft)%first(n) = ntgrid
                from(ipleft)%second(n) = 2
                from(ipleft)%third(n) = iglo
             end if

             if (ipleft == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n

                to(ip)%first(n) = -ntgrid
                to(ip)%second(n) = 2
                to(ip)%third(n) = iglo_left
             end if
          end if

          ! If we're at the left end we need to pass to the rightmost
          if (it == it_leftmost) then
             iglo_right = idx(g_lo, ik, it_rightmost, il, ie, is)
             ipright = proc_id(g_lo, iglo_right)

             if (ip == iproc) then
                n = nn_from(ipright) + 1
                nn_from(ipright) = n

                from(ipright)%first(n) = -ntgrid
                from(ipright)%second(n) = 1
                from(ipright)%third(n) = iglo
             end if

             ! receiver
             if (ipright == iproc) then
                n = nn_to(ip) + 1
                nn_to(ip) = n

                to(ip)%first(n) = ntgrid
                to(ip)%second(n) = 1
                to(ip)%third(n) = iglo_right
             end if
          end if
       end if
    end do

    from_low = [-ntgrid, 1, g_lo%llim_proc]
    from_high = [ntgrid, 2, g_lo%ulim_alloc]

    to_low = from_low
    to_high = from_high

    call init_fill (incoming_links, 'c', to_low, to_high, to, &
         from_low, from_high, from)

    call delete_list(from)
    call delete_list(to)

  end subroutine setup_pass_incoming_boundary_to_connections

  !> Coordinates the calculation of all data required for later application
  !> of linked boundary conditions. Deals with identifying the communication
  !> pattern and creating the associated redistributes used to perform the
  !> required communication.
  subroutine init_connected_bc
    use kt_grids, only: naky, ntheta0
    use le_grids, only: mixed_wfb
    use gs2_layouts, only: g_lo
    use array_utils, only: zero_array
    use mp, only: proc0
    use file_utils, only: error_unit
    implicit none
    integer :: jshift0, n_links_max

    if (connectinit) return
    connectinit = .true.

    ! Skip setup of linked bc data in cases without linked boundaries.
    ! We might want to remove this check so that we can calculate
    ! the linked bc data for testing etc.
    if (boundary_option_switch /= boundary_option_linked) return

    ! Get the kx index spacing between connections
    jshift0 = compute_jshift0()

    allocate (itleft(ntheta0, naky), itright(ntheta0, naky))
    call compute_itleft_and_itright(jshift0, itleft, itright)

    allocate (connections(g_lo%llim_proc:g_lo%ulim_alloc))
    call compute_connections(itleft, itright, connections)

    allocate (save_h(2, g_lo%llim_proc:g_lo%ulim_alloc))
    call set_save_h(connections, save_h)

    call count_links(itleft, itright, l_links, r_links, n_links)
    n_links_max = maxval(n_links)

    ! wfb -- mixed treatment requires some extra links to be communicated
    if (n_links_max > 0 .and. mixed_wfb) n_links_max = n_links_max + 3

    ! Could we move this call to init_bc? Would have to check that it returns
    ! the expected values when l_links = r_links = 0
    call setup_fields_implicit_class_arrays(l_links, r_links, M_class, N_class, i_class)

    ! now set up communication pattern:
    ! excluding wfb
    call setup_connected_bc_redistribute(l_links, r_links, n_links_max, &
         links_p, links_h, no_connections)

    ! If no connections then don't need to do any more work here. One
    ! might expect us to be able to set no_connections = n_links_max == 0
    ! rather than having to call setup_connected_bc_redistribute to
    ! determine this. This would allow a little bit of work to be skipped.
    ! Note that no_connections is a global quantity (i.e. agreed value
    ! on all processors). We may wish to introduce a local version as well.
    ! At scale each processor may only be responsible for one pitch angle,
    ! for example. If this is a trapped pitch angle then it will not have
    ! any work to do associated with the boundaries yet it will still go
    ! through the allocations etc. associated with setting up linked boundaries
    if (no_connections) then
       if (proc0) write(error_unit(), &
            '("No connections found when setting up linked boundaries.")')
       return
    end if

    ! take care of wfb
    if (mixed_wfb) then
       call setup_connected_bc_redistribute_mixed_wfb(l_links, r_links, n_links_max, &
            wfb_p, wfb_h)
    end if

    ! Setup redistributes for passing the incoming boundary to the previous linked
    ! cell. Used to generalise the linked boundary conditions to allow gnew /= 0 at
    ! the linked location.
    if (start_from_previous_solution) then
       call setup_pass_incoming_boundary_to_connections(l_links, r_links, incoming_links)
    end if

    ! n_links_max is typically 2 * number of cells in largest supercell
    ! Note this is potentially quite wasteful in terms of storage. In flux
    ! tube setup there is likely only one supercell with the maximal number
    ! of links and many more with fewer (or zero) links yet we allocate the
    ! full length for _all_ iglo points. For the longest supercell we expect
    ! n_links ~ 2 * ntheta0/jtwist - 1 and hence n_links_max is ~ 2*( 1 + ntheta/jtwist)
    ! For simulations with high nx and/or low jtwist n_links_max could be >> ntheta
    ! and therefore g_adj may dominate memory consumption in some cases.
    allocate (g_adj (n_links_max, 2, g_lo%llim_proc:g_lo%ulim_alloc))
    call zero_array(g_adj)

  end subroutine init_connected_bc

  !> Copy of CMR's init_pass_right (revision 2085) 
  !! but in customisable direction to aid code reuse and with optional
  !!range in theta (i.e. so can pass more than just boundary points)
  !! An example use is:
  !! call init_pass_ends(pass_right,'r',1,'c',3)
  !! call fill(pass_right,gnew,gnew)
  !! This will pass gnew(ntgrid-2:ntgrid,1,iglo) to 
  !! gnew(-ntgrid-2:-ntgrid,1,iglo_linked)
  !!
  !! @TODO:
  !!      1. May be helpful to be able to pass both sigmas at once (e.g. for explicit scheme)
  !! DD 01/02/2013
  subroutine init_pass_ends(pass_obj,dir,sigma,typestr,ngsend)
    use gs2_layouts, only: g_lo, il_idx, idx, proc_id
    use le_grids, only: il_is_trapped
    use mp, only: iproc, nproc, max_allreduce, proc0
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use theta_grid, only:ntgrid

    implicit none

    type (redist_type), intent(out) :: pass_obj   !< Redist type object to hold communication logic
    character(1),intent(in)::dir                  !< Character string for direction of communication, should be 'l' for left and 'r' for right
    character(1),intent(in)::typestr              !< Character string for type of data to be communicated. Should be 'c','r','i' or 'l'
    integer,intent(in) :: sigma                   !< Which sigma index to send
    integer,intent(in),optional::ngsend           !< How many theta grid points are we sending

    !Internal variables
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension(3) :: from_low, from_high, to_low, to_high
    integer :: il, iglo, ip, iglo_con, ipcon, n, nn_max, j
    logical,parameter :: debug=.false.
    integer :: bound_sign
    integer :: local_ngsend

    !Only applies to linked boundary option so exit if not linked
    if (boundary_option_switch .ne. boundary_option_linked) return

    !Handle the direction sign, basically we're either doing
    !     ntgrid --> -ntgrid (passing to right)
    !or 
    !    -ntgrid --> ntgrid  (passing to left)
    if (dir=='r') then
       bound_sign=1
    elseif (dir=='l') then
       bound_sign=-1
    else
       if (proc0) write(6,*) "Invalid direction string passed to init_pass_ends, defaulting to 'r'"
       bound_sign=1
    endif

    !Set the default number of theta grid points to send, if required
    if (present(ngsend)) then
       local_ngsend=ngsend
    else
       local_ngsend=1
    endif

    if (proc0.and.debug) write (6,*) "Initialising redist_type with settings Direction : ",dir," sigma",sigma,"local_ngsend",local_ngsend
    
    ! Need communications to satisfy || boundary conditions
    ! First find required blocksizes 
    
    !Initialise variables used to count how many entries to send and receive
    !for each processor
    nn_to = 0   ! # nn_to(ip) = communicates from ip TO HERE (iproc)
    nn_from = 0 ! # nn_from(ip) = communicates to ip FROM HERE (iproc)
    
    !Now loop over >all< iglo indices and work out how much data needs to be sent and received by each processor
    !Hence this routine does not scale with number of procs, see updated redist object creation in ccfe_opt_test (~r2173)
    !for examples which do scale
    do iglo = g_lo%llim_world, g_lo%ulim_world
       !Get the lambda index so we can skip trapped particles
       il = il_idx(g_lo,iglo)

       !Exclude disconnected trapped particles
       if (il_is_trapped(il)) cycle

       !Get the processor id for the proc holding the current iglo index
       ip = proc_id(g_lo,iglo)
       
       !What iglo connects to the current one in the direction of interest (and what proc is it on)
       !Note ipcon is <0 if no connection in direction of interest
       if (bound_sign==1) then
          call get_right_connection (iglo, iglo_con, ipcon)
       else
          call get_left_connection (iglo, iglo_con, ipcon)
       endif
       
       !Is the connected tube's proc the current processor?
       if (ipcon .eq. iproc ) then
          !If so add an entry recording an extra piece of information is to be sent
          !to proc ip (the proc holding iglo) from this proc (ipcon)
          !Note: Here we assume theta grid points are all on same proc 
          nn_to(ip)=nn_to(ip)+local_ngsend
       endif
       
       !Is the proc holding iglo this proc and is there a connection in the direction of interest?
       if (ip .eq. iproc .and. ipcon .ge. 0 ) then
          !If so add an entry recording an extra piece of information is to be sent
          !from this proc (ipcon) to ip
          !Note: Here we assume theta grid points are all on same proc 
          nn_from(ipcon)=nn_from(ipcon)+local_ngsend
       endif
    end do
    
    !Find the maximum amount of data to be received by a given processor
    !(first do it locally)
    nn_max = maxval(nn_to)
    !(now do it globally)
    call max_allreduce (nn_max)
    
    !Bit of debug printing
    if (proc0.and.debug) then
       write(6,*) 'init_pass_ends (1) processor, nn_to, nn_from:',iproc,nn_to,nn_from
    endif

    !Now that we've worked out how much data needs to be sent and received, define what specific
    !data needs to be sent to where
    if (nn_max.gt.0) then
       ! 
       ! CMR, 25/1/2013: 
       !      communication required to satisfy linked BC
       !      allocate indirect addresses for sends/receives 
       !     
       !      NB communications use "c_fill_3" as g has 3 indices
       !      but 1 index sufficient as only communicating g(ntgrid,1,*)! 
       !      if "c_fill_1" in redistribute we'd only need allocate: from|to(ip)%first 
       !                             could be more efficient
       !  
       !<DD>, 06/01/2013: This redist object consists of a buffer of length n to hold the 
       !data during transfer and (currently) 3*2 integer arrays each of length n to hold
       !the indices of sent and received data.
       !By using c_fill_1 2*2 of these integer arrays would be removed. Assuming a double complex
       !buffer and long integer indices a 4n long array saving would be equivalent to the buffer
       !size and as such should represent a good memory saving but would not effect the amount
       !of data communicated (obviously).

       !Create to and from list objects for each processor and 
       !create storage to hold information about each specific from/to
       !communication
       do ip = 0, nproc-1
          !If proc ip is sending anything to this processor (iproc)
          if (nn_from(ip) > 0) then
             allocate (from(ip)%first(nn_from(ip)))
             allocate (from(ip)%second(nn_from(ip)))
             allocate (from(ip)%third(nn_from(ip)))
          endif
          !If proc ip is receiving anything from this processor (iproc)
          if (nn_to(ip) > 0) then
             allocate (to(ip)%first(nn_to(ip)))
             allocate (to(ip)%second(nn_to(ip)))
             allocate (to(ip)%third(nn_to(ip)))
          endif
       end do

       !Now fill the indirect addresses...
       
       !Initialise counters used to record how many pieces of data to expect
       nn_from = 0 ; nn_to = 0
       
       !Loop over >all< iglo indices
       do iglo = g_lo%llim_world, g_lo%ulim_world
          !Get the lambda index so we can skip trapped particles
          il = il_idx(g_lo,iglo)
          
          !Exclude disconnected trapped particles
          if (il_is_trapped(il)) cycle

          !What's the processor for the current iglo
          ip = proc_id(g_lo,iglo)
          
          !What iglo connects to the current one in the direction of interest (and what proc is it on)?
          !Note ipcon is <0 if no connection in direction of interest
          if (bound_sign==1) then
             call get_right_connection (iglo, iglo_con, ipcon)
          else
             call get_left_connection (iglo, iglo_con, ipcon)
          endif
          
          !For current proc for current iglo if there's connections in direction
          !then add an entry to the connected procs list of data to expect
          if (ip .eq. iproc .and. ipcon .ge. 0 ) then
             !Loop over theta grid indices
             !Note: Here we assume theta grid points are all on same proc 
             !Note: By looping over theta inside iglo loop we should optimise
             !      memory access compared to looping over theta outside.
             do j=0,local_ngsend-1
                n=nn_from(ipcon)+1 ; nn_from(ipcon)=n
                from(ipcon)%first(n)=bound_sign*(ntgrid-j) !Which theta point to send
                from(ipcon)%second(n)=sigma     !Sigma grid index to send
                from(ipcon)%third(n)=iglo   !iglo index to send
             enddo
          endif
          
          !If target iglo (iglo_con) is on this processor then add an entry recording where
          !we need to put the data when we receive it.
          if (ipcon .eq. iproc ) then 
             !Loop over theta grid indices
             !Note: Here we assume theta grid points are all on same proc 
             do j=0,local_ngsend-1
                n=nn_to(ip)+1 ; nn_to(ip)=n
                to(ip)%first(n)=-bound_sign*(ntgrid-j) !Which theta to store received data
                to(ip)%second(n)=sigma       !Sigma grid index to store received data
                to(ip)%third(n)=iglo_con !iglo index to store received data
             enddo
          endif
       end do
       
       !Bit of debug printing,
       if (debug.and.proc0) then
          write(6,*) 'init_pass_ends (2) processor, nn_to, nn_from:',iproc,nn_to,nn_from
       endif

       !Set data ranges for arrays to be passed, not this just effects how
       !arrays are indexed, not how big the buffer is.
       from_low(1)=-ntgrid ; from_low(2)=1  ; from_low(3)=g_lo%llim_proc       
       from_high(1)=ntgrid ; from_high(2)=2 ; from_high(3)=g_lo%ulim_alloc
       to_low(1)=-ntgrid   ; to_low(2)=1    ; to_low(3)=g_lo%llim_proc       
       to_high(1)=ntgrid   ; to_high(2)=2   ; to_high(3)=g_lo%ulim_alloc
       
       !Initialise fill object
       call init_fill (pass_obj, typestr, to_low, to_high, to, &
            from_low, from_high, from)
       
       !Clean up lists
       call delete_list (from)
       call delete_list (to)
    endif
  end subroutine init_pass_ends

  !> FIXME : Add documentation
  subroutine get_left_connection (iglo, iglo_left, iproc_left)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_left, iproc_left
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itleft(it, ik) < 0) then
       iglo_left = -1
       iproc_left = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_left = idx(g_lo,ik,itleft(it, ik),il,ie,is)
    iproc_left = proc_id(g_lo,iglo_left)
  end subroutine get_left_connection

  !> FIXME : Add documentation  
  subroutine get_right_connection (iglo, iglo_right, iproc_right)
    use gs2_layouts, only: g_lo, proc_id, idx
    use gs2_layouts, only: ik_idx, it_idx, il_idx, ie_idx, is_idx
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_right, iproc_right
    integer :: ik, it, il, ie, is

    ik = ik_idx(g_lo,iglo)
    it = it_idx(g_lo,iglo)
    
    if (itright(it, ik) < 0) then
       iglo_right = -1
       iproc_right = -1
       return
    end if

    il = il_idx(g_lo,iglo)
    ie = ie_idx(g_lo,iglo)
    is = is_idx(g_lo,iglo)

    iglo_right = idx(g_lo,ik,itright(it, ik),il,ie,is)
    iproc_right = proc_id(g_lo,iglo_right)
  end subroutine get_right_connection

  !> Helper function for finding the leftmost it of supercell
  elemental integer function get_leftmost_it(it,ik) result(it_cur)
    implicit none
    integer, intent(in) :: it, ik
    !If not linked then no connections so only one cell in supercell
    if (boundary_option_switch == boundary_option_linked) then
       it_cur = it
       do while (itleft(it_cur, ik) >= 0 .and. itleft(it_cur, ik) /= it_cur)
          !Keep updating it_cur with left connected it until there are no
          !connections to left
          it_cur = itleft(it_cur, ik)
       end do
    else
       it_cur = it
    end if
  end function get_leftmost_it

  !> Helper function for finding the rightmost it of supercell
  elemental integer function get_rightmost_it(it, ik) result(it_cur)
    implicit none
    integer, intent(in) :: it, ik
    !If not linked then no connections so only one cell in supercell
    if (boundary_option_switch == boundary_option_linked) then
       it_cur = it
       do while (itright(it_cur, ik) >= 0 .and. itright(it_cur, ik) /= it_cur)
          !Keep updating it_cur with right connected it until there are no
          !connections to right
          it_cur = itright(it_cur, ik)
       end do
    else
       it_cur = it
    end if
  end function get_rightmost_it
  
  !> FIXME : Add documentation
  subroutine allocate_arrays
    use kt_grids, only: ntheta0, naky,  box
    use array_utils, only: zero_array
    use theta_grid, only: ntgrid, shat
    use dist_fn_arrays, only: g, gnew, g_work
    use dist_fn_arrays, only: kx_shift, theta0_shift   ! MR
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3
    use dist_fn_arrays, only: antot, antota, antotp
    use dist_fn_arrays, only: fieldeq, fieldeqa, fieldeqp
    use gs2_layouts, only: g_lo
    use nonlinear_terms, only: nonlin, split_nonlinear
    use run_parameters, only: has_apar
    use array_utils, only: zero_array
#ifdef SHMEM
    use shm_mpi3, only : shm_alloc
#endif
    implicit none

    if (.not. allocated(g)) then
       allocate (g      (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gnew   (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g_work (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (g_h  (-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       call zero_array(g) ; call zero_array(gnew) ; call zero_array(g_work)
       call zero_array(g_h)

       if(opt_source)then
          if(has_apar)then
             allocate (source_coeffs(4,-ntgrid:ntgrid-1,2,g_lo%llim_proc:g_lo%ulim_alloc))
          else
             allocate (source_coeffs(2,-ntgrid:ntgrid-1,2,g_lo%llim_proc:g_lo%ulim_alloc))
          endif
       endif
#ifdef LOWFLOW
#ifndef SHMEM
          allocate (gexp_1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
#else
          call shm_alloc(gexp_1, (/ -ntgrid, ntgrid, 1, 2, &
               g_lo%llim_proc, g_lo%ulim_alloc/))
#endif
       allocate (gexp_2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       allocate (gexp_3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
       call zero_array(gexp_1) ; call zero_array(gexp_2) ; call zero_array(gexp_3)
#else
       if (nonlin .and. .not. split_nonlinear) then
#ifndef SHMEM
          allocate (gexp_1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
#else
          call shm_alloc(gexp_1, (/ -ntgrid, ntgrid, 1, 2, &
               g_lo%llim_proc, g_lo%ulim_alloc/))
#endif
          allocate (gexp_2(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          allocate (gexp_3(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
          call zero_array(gexp_1) ; call zero_array(gexp_2) ; call zero_array(gexp_3)
       end if
#endif

       if (abs(g_exb*g_exbfac) > epsilon(0.)) then           ! MR
          if (box .or. shat .eq. 0.0) then
             allocate (kx_shift(naky))
             kx_shift = 0.
          else
             allocate (theta0_shift(naky))
             theta0_shift = 0.
          endif
       endif                           ! MR end
    endif

    allocate (antot (-ntgrid:ntgrid,ntheta0,naky))
    allocate (antota(-ntgrid:ntgrid,ntheta0,naky))
    allocate (antotp(-ntgrid:ntgrid,ntheta0,naky))
    call zero_array(antot) ; call zero_array(antota) ; call zero_array(antotp)

    allocate (fieldeq (-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqa(-ntgrid:ntgrid,ntheta0,naky))
    allocate (fieldeqp(-ntgrid:ntgrid,ntheta0,naky))
    call zero_array(fieldeq) ; call zero_array(fieldeqa) ; call zero_array(fieldeqp)
  end subroutine allocate_arrays

  !> This function calculates the distribution function at the next timestep. 
  !! It calculates both the inhomogeneous part, gnew, due to the sources
  !! (principly the drive terms and the nonlinear term)
  !! and the homogeneous part, g1. The actual evolution of the dist func
  !! is done in the subroutine invert_rhs. 
  !!
  !! After solving for the new dist funcs, this routine calls hyper_diff, which
  !! adds hyper diffusion if present, and solfp1, from the collisions module,
  !! which adds collisions if present.
  subroutine timeadv (phi, apar, bpar, phinew, aparnew, bparnew, istep, mode)
    use theta_grid, only: ntgrid
    use le_derivatives, only: vspace_derivatives
    use dist_fn_arrays, only: gnew, g, g_work
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3
    use split_nonlinear_terms, only: advance_nonlinear_term, advance_nonadiabatic_dfn
    use nonlinear_terms, only: add_explicit_terms, split_nonlinear, calculate_current_nl_source_and_cfl_limit, nonlin
    use hyper, only: hyper_diff
    use run_parameters, only: reset, immediate_reset
    use unit_tests, only: debug_message
    use collisions, only: split_collisions
    use gs2_time, only: code_dt
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    integer, optional, intent (in) :: mode
    integer, save :: istep_last = -1
    integer, parameter :: verb = 3

    !Calculate the explicit nonlinear terms
    if (split_nonlinear .and. nonlin) then
       ! If we're starting a new full step then advance just the nonlinear term
       ! by code_dt.
       if (istep > 0 .and. istep /= istep_last) then
          if (advance_nonadiabatic_dfn) then
             call advance_nonlinear_term(g, istep, phi, apar, bpar, code_dt, &
                  calculate_potentials_from_nonadiabatic_dfn, &
                  calculate_current_nl_source_and_cfl_limit)
          else
             call advance_nonlinear_term(g, istep, phi, apar, bpar, code_dt, &
                  get_fields_direct_from_dfn, calculate_current_nl_source_and_cfl_limit)
          end if
       end if
    else
       call add_explicit_terms (gexp_1, gexp_2, gexp_3, &
            phi, apar, bpar, istep, bkdiff(1))
    end if

    call debug_message(verb, &
        'dist_fn::timeadv calculated nonlinear term')
    if(reset .and. immediate_reset) return !Return if resetting
    !Solve for gnew
    call invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    call debug_message(verb, &
        'dist_fn::timeadv calculated rhs')

    !Add hyper terms (damping)
    if ((istep == 0 .and. hyper_in_initialisation) .or. istep > 0) then
       !Note we potentially exclude the hyper terms when we calculate the response
       !data (indicated by istep==0) as these terms are typically not
       !linearly independent.
       call hyper_diff (gnew, phinew)
       call debug_message(verb, &
            'dist_fn::timeadv calculated hypviscosity')
    end if

    !Add collisions
    if( .not. split_collisions) then
      call vspace_derivatives (gnew, g, g_work, phi, bpar, phinew, bparnew, mode)
      call debug_message(verb, &
          'dist_fn::timeadv calculated collisions etc')
    end if

    !Enforce parity if desired (also performed in invert_rhs, but this is required
    !as collisions etc. may break parity?)
    if (def_parity) call enforce_parity(parity_redist)
      
    istep_last=istep

  end subroutine timeadv

  !> Advance collision term one timestep by implicit Euler.
  !!
  !! Advance g under collisions one timestep by implicit Euler, then 
  !! solve Maxwell's equations to find consistent fields.
  !! This is only performed if collisions are not evolved as part of the 
  !! usual time advance in timeadv (i.e. if split_collisions=.true.).
  !! Splitting collisions also allows tricks like:
  !!   + not performing redistributes between collisions and the field solve 
  !!     if both are done in gf. 
  !!   + only applying collisions every nth timestep
  !!
  !!  Author : Joseph Parker, STFC
  subroutine collisions_advance (phi, bpar, phinew, aparnew, bparnew, mode)
    use collisions, only: split_collisions
    use dist_fn_arrays, only: gnew, g, g_work
    use le_derivatives, only: vspace_derivatives
    use theta_grid, only: ntgrid
    use unit_tests, only: debug_message
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    complex, dimension (-ntgrid:,:,:), intent (inout) :: phinew, aparnew, bparnew
    integer, optional, intent (in) :: mode
    integer, parameter :: verb = 3

    if(.not. split_collisions) return

    !Add collisions
    call vspace_derivatives (gnew, g, g_work, phi, bpar, phinew, bparnew, mode)
    call debug_message(verb, &
      'dist_fn::collisions_advance calculated collisions etc')

    !Find fields for new distribution function
    call get_init_field (phinew, aparnew, bparnew)

  end subroutine collisions_advance

! communication initializations for exb_shear should be done once and 
! redistribute routines should be used  BD
!
!  subroutine init_exb_shear
!
!    use redistribute, only: index_list_type, delete_list
!    implicit none
!    type (index_list_type), dimension(0:nproc-1) :: to, from
!    integer, dimension (0:nproc-1) :: nn_from, nn_to
!    integer, dimension (3) :: to_low, from_low, to_high, from_high
!
!  end subroutine init_exb_shear

  !> FIXME : Add documentation
  subroutine exb_shear (g0, phi, apar, bpar, istep, field_local)
! MR, 2007: modified Bill Dorland's version to include grids where kx grid
!           is split over different processors
! MR, March 2009: ExB shear now available on extended theta grid (ballooning)
! CMR, May 2009: 2pishat correction factor on extended theta grid (ballooning)
!                so GEXB is same physical quantity in box and ballooning
! CMR, Oct 2010: multiply timestep by tunits(iky) for runs with wstar_units=.t.
! CMR, Oct 2010: add save statements to prevent potentially large and memory 
!                killing array allocations!
   
    use mp, only: send, receive, mp_abort, broadcast
    use gs2_layouts, only: ik_idx, it_idx, g_lo, idx_local, idx, proc_id
    use theta_grid, only: ntgrid, ntheta, shat
    use file_utils, only: error_unit
    use kt_grids, only: akx, aky, naky, ikx, ntheta0, box, theta0
    use le_grids, only: negrid, nlambda
    use species, only: nspec
    use run_parameters, only: has_phi, has_apar, has_bpar
    use dist_fn_arrays, only: kx_shift, theta0_shift
    use gs2_time, only: code_dt, code_dt_old, code_time, tunits
    use constants, only: twopi
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: g0
    integer, intent(in) :: istep
    logical, intent(in), optional :: field_local
    complex, dimension(:,:,:), allocatable :: temp 
    complex, dimension(:,:), allocatable :: temp2
    integer, dimension(1), save :: itmin
    logical :: should_use_kx
    integer :: ierr, j 
    integer :: ik, it, ie, is, il, to_iglo, from_iglo
    integer:: iib, iit, ileft, iright, i
    integer, save :: istep_last = 0
    logical :: field_local_loc
    real, save :: dkx, dtheta0
    real :: gdt
    complex , dimension(-ntgrid:ntgrid, 2) :: z
    character(130) :: str

    ierr = error_unit()

    ! If in flux-tube (box mode) or have zero shear then we need
    ! to work with kx, otherwise use theta0
    should_use_kx = box .or. (shat .eq. 0.0)

! MR, March 2009: remove ExB restriction to box grids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MR, 2007: Works for ALL layouts (some layouts need no communication)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize kx_shift, jump, idx_indexed
    if (exb_first) then
       exb_first = .false.
       allocate (jump(naky)) 
       jump = 0
       if (should_use_kx) then
          allocate (ikx_indexed(ntheta0))
          itmin = minloc (ikx)
          
          do it=itmin(1), ntheta0
             ikx_indexed (it+1-itmin(1)) = it
          end do
          
          do it=1,itmin(1)-1
             ikx_indexed (ntheta0 - itmin(1) + 1 + it)= it
          end do

          if (ntheta0 .gt. 1) then 
             dkx = akx(2)-akx(1)
          else
             write(ierr,*) "exb_shear: ERROR, need ntheta0>1 for sheared flow"
             !Warning, dkx has not been set => Should really halt run or find
             !a suitable default definition for dkx.
          endif
       else
! MR, March 2009: on extended theta grid theta0_shift tracks ExB shear
! CMR, 25 May 2009: fix 2pi error so that: dtheta0/dt = -GEXB/shat
          if (ntheta0 .gt. 1) then 
             dtheta0 = theta0(2,1)-theta0(1,1) 
          else 
             dtheta0=twopi
          endif
       end if
    end if

    !Check if we want to exit without applying flow shear
    if (istep.eq.istep_last) return !Don't allow flow shear to happen more than once per step    
    if (g_exb_start_timestep > istep) return !Flow shear not active yet, set in timesteps
    if (g_exb_start_time >= 0 .and. code_time < g_exb_start_time) return !Flow shear not active yet, set in time
    
    ! BD: To do: Put the right timestep in here.
    ! For now, approximate Greg's dt == 1/2 (t_(n+1) - t_(n-1))
    ! with code_dt.  
    !
    ! Note: at first time step, there is a difference of a factor of 2.
    !
    ! necessary to get factor of 2 right in first time step and
    ! also to get things right after changing time step size
    ! added May 18, 2009 -- MAB
    gdt = 0.5*(code_dt + code_dt_old)

    !Update istep_last
    istep_last = istep

! kx_shift is a function of time.   Update it here:  
! MR, 2007: kx_shift array gets saved in restart file
! CMR, 5/10/2010: multiply timestep by tunits(ik) for wstar_units=.true. 
    ! Here we calculate the continuous kx/theta0 shift away from the current
    ! grid point we are on. Once this reaches at lesat +1/2 * kx/theta0 spacing
    ! jump will be non-zero so we will need to shuffle data. We then subtract
    ! jump times the grid space from our continuous shift to reset this to
    ! effectively -1/2 * kx/theta0 spacing. Note that as this array starts at
    ! 0, the first jump will happen in half the time of subsequent jumps.
    !
    ! @note These statements should be verified.
    if (should_use_kx) then
       do ik=1, naky
          kx_shift(ik) = kx_shift(ik) - aky(ik)*g_exb*g_exbfac*gdt*tunits(ik)
          jump(ik) = nint(kx_shift(ik)/dkx)
          kx_shift(ik) = kx_shift(ik) - jump(ik)*dkx
       end do
    else
       do ik=1, naky
          theta0_shift(ik) = theta0_shift(ik) - g_exb*g_exbfac*gdt/shat*tunits(ik)
          jump(ik) = nint(theta0_shift(ik)/dtheta0)
          theta0_shift(ik) = theta0_shift(ik) - dtheta0*jump(ik)
       enddo 
    end if

    !If jump is zero for all ky then nothing else to do so leave
    if (all(jump == 0)) return

    !If using field_option='local' and x/it is not entirely local
    !then we need to make sure that all procs know the full field
    !THIS IS A TEMPORARY FIX AND WE SHOULD PROBABLY DO SOMETHING BETTER
    field_local_loc = get_option_with_default(field_local, .false.)
    if(any(jump.ne.0).and.(.not.g_lo%x_local).and.field_local_loc) then
       if(has_phi) call broadcast(phi)
       if(has_apar) call broadcast(apar)
       if(has_bpar) call broadcast(bpar)
    endif

    ! We should note that the following branches contain loops over the entire
    ! global g_lo domain in order to do ordered point-wise communication. We should
    ! probably think about constructing redistribute objects to represent the relevant
    ! communications instead. As the communication pattern is ky dependent (as jump
    ! depends on ky) and we may not wish to communicate in every ky at the same time
    ! we probably want naky redistribute objects. This may be less than optimal as
    ! we cannot then group communications for all ky together as in other redistributes.
    ! The other option is to construct a single redistribute each time one is needed.
    ! This would allow consolidation of communications but with the overhead of
    ! redistribute construction. For cases with should_use_kx true, we would not
    ! expect to need to communicate different ky at the same time (i.e. jump is
    ! unlikely to be non-zero for more than one ky at a time). When should_use_kx
    ! is false we will always need to communicate all ky at the same time as jump
    ! is ky independent (unless wstar_units is true).

    if (.not. should_use_kx) then
! MR, March 2009: impact of ExB shear on extended theta grid computed here
!                 for finite shat
       do ik =1,naky
          j=jump(ik)
          if (j .eq. 0) cycle     
          if (abs(j) .ge. ntheta0) then
              write(str,fmt='("in exb_shear: jump(ik)=",i4," > ntheta0 =",i4," for ik=",i4,". => reduce timestep or increase ntheta0")') j,ik,ntheta0
              write(ierr,*) str
              call mp_abort(str)
          endif 
          allocate(temp2(-ntgrid:ntgrid,abs(j)),temp(-ntgrid:ntgrid,2,abs(j)))
          iit=ntheta0+1-abs(j) ; iib=abs(j)
          ileft = -ntgrid+ntheta ; iright=ntgrid-ntheta

          if (has_phi) then
             if (j < 0) then
                temp2 = phi(:,:iib,ik)
                do i=1,iit-1
                   phi(:,i,ik) = phi(:,i-j,ik)
                enddo
                phi(ileft:,iit:,ik) = temp2(:iright,:)
                phi(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = phi(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   phi(:,i,ik) = phi(:,i-j,ik)
                enddo
                phi(:iright,:iib ,ik) = temp2(ileft:,:)
                phi(iright+1:ntgrid,:iib,:) = 0.0
             endif
          endif
          if (has_apar) then
             if (j < 0) then
                temp2 = apar(:,:iib,ik)
                do i=1,iit-1
                   apar(:,i,ik) = apar(:,i-j,ik)
                enddo
                apar(ileft:,iit:,ik) = temp2(:iright,:)
                apar(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = apar(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   apar(:,i,ik) = apar(:,i-j,ik)
                enddo
                apar(:iright,:iib ,ik) = temp2(ileft:,:)
                apar(iright+1:ntgrid,:iib,:) = 0.0
             endif
          endif
          if (has_bpar) then
             if (j < 0) then
                temp2 = bpar(:,:iib,ik)
                do i=1,iit-1
                   bpar(:,i,ik) = bpar(:,i-j,ik)
                enddo
                bpar(ileft:,iit:,ik) = temp2(:iright,:)
                bpar(:ileft-1,iit:,ik) = 0.0
             else 
                temp2 = bpar(:,iit:,ik)
                do i=ntheta0,iib+1,-1 
                   bpar(:,i,ik) = bpar(:,i-j,ik)
                enddo
                bpar(:iright,:iib ,ik) = temp2(ileft:,:)
                bpar(iright+1:ntgrid,:iib,:) = 0.0
             endif
          end if

! now the distribution functions

          do is=1,nspec
             do ie=1,negrid
                do il=1,nlambda

                   if (j < 0) then
                      do it = 1, iib
                         from_iglo = idx(g_lo, ik, it, il, ie, is)
                         if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                      end do

                      do it = 1, iit-1                        
                           to_iglo = idx(g_lo, ik, it,   il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j, il, ie, is)

                         if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            call send(g0(:, :, from_iglo), proc_id (g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive(g0(:, :, to_iglo), proc_id (g_lo, from_iglo))
                         endif
                      enddo

                      do it = iit, ntheta0                     
                         to_iglo = idx(g_lo, ik, it, il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j-ntheta0, il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(ileft:,:,to_iglo) = temp(:iright,:,it-iit+1)
                            g0(:ileft-1,:,to_iglo) = 0.0
                         else if (idx_local(g_lo, from_iglo)) then
                            call send(temp(:, :, it-iit+1), proc_id (g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive(z, proc_id (g_lo, from_iglo))
                            g0(ileft:, :, to_iglo) = z(:iright, :)
                            g0(:ileft-1, :, to_iglo) = 0.0
                         endif
                      enddo

                   else ! j>0

                      do it = 1, j
                         from_iglo = idx(g_lo, ik, iit+it-1, il, ie, is)
                         if (idx_local (g_lo, from_iglo)) temp(:,:,it) = g0(:,:,from_iglo)
                      end do

                      do it = ntheta0, iib+1, -1
                           to_iglo = idx(g_lo, ik, it,   il, ie, is)
                         from_iglo = idx(g_lo, ik, it-j, il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            call send(g0(:, :, from_iglo), proc_id (g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive(g0(:, :, to_iglo), proc_id (g_lo, from_iglo))
                         endif
                      enddo

                      do it = 1, iib
                           to_iglo = idx(g_lo, ik, it,           il, ie, is)
                         from_iglo = idx(g_lo, ik, iit+it-1, il, ie, is)

                         if (idx_local(g_lo, to_iglo).and. idx_local(g_lo, from_iglo)) then
                            g0(:iright,:,to_iglo) = temp(ileft:,:,it)
                            g0(iright+1:,:,to_iglo) = 0.0
                         else if (idx_local(g_lo, from_iglo)) then
                            call send(temp(:, :, it), proc_id (g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive(z, proc_id (g_lo, from_iglo))
                            g0(:iright, :, to_iglo) = z(ileft:, :)
                            g0(iright+1:, :, to_iglo) = 0.0
                         endif
                      enddo
                   endif
                enddo
             enddo
          enddo
          deallocate (temp,temp2)
       enddo
    else
       ! should_use_kx = .true.
       do ik = naky, 2, -1
          if (jump(ik) < 0) then
             if (has_phi) then
                do it = 1, ntheta0 + jump(ik)
                   phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   phi(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (has_apar) then
                do it = 1, ntheta0 + jump(ik)
                   apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   apar (:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (has_bpar) then
                do it = 1, ntheta0 + jump(ik)
                   bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = ntheta0 + jump(ik) + 1, ntheta0
                   bpar (:,ikx_indexed(it),ik) = 0.
                end do
             end if
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda

                      do it = 1, ntheta0 + jump(ik)                        

                           to_iglo = idx(g_lo, ik, ikx_indexed(it),          il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            call send (g0(:, :, from_iglo), proc_id (g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive (g0(:, :, to_iglo), proc_id (g_lo, from_iglo))
                         endif
                      enddo

                      do it = ntheta0 + jump(ik) + 1, ntheta0                     
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
                      enddo

                   enddo
                enddo
             enddo
          endif

          if (jump(ik) > 0) then 
             if (has_phi) then
                do it = ntheta0, 1+jump(ik), -1
                   phi(:,ikx_indexed(it),ik) = phi(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   phi(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (has_apar) then
                do it = ntheta0, 1+jump(ik), -1
                   apar(:,ikx_indexed(it),ik) = apar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   apar(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             if (has_bpar) then
                do it = ntheta0, 1+jump(ik), -1
                   bpar(:,ikx_indexed(it),ik) = bpar(:,ikx_indexed(it-jump(ik)),ik)
                end do
                do it = jump(ik), 1, -1
                   bpar(:,ikx_indexed(it),ik) = 0.
                end do
             end if
             do is=1,nspec
                do ie=1,negrid
                   do il=1,nlambda

                      do it = ntheta0, 1+jump(ik), -1

                           to_iglo = idx(g_lo, ik, ikx_indexed(it),          il, ie, is)
                         from_iglo = idx(g_lo, ik, ikx_indexed(it-jump(ik)), il, ie, is)

                         if (idx_local(g_lo, to_iglo) .and. idx_local(g_lo, from_iglo)) then
                            g0(:,:,to_iglo) = g0(:,:,from_iglo)
                         else if (idx_local(g_lo, from_iglo)) then
                            call send(g0(:, :, from_iglo), proc_id(g_lo, to_iglo))
                         else if (idx_local(g_lo, to_iglo)) then
                            call receive(g0(:, :, to_iglo), proc_id (g_lo, from_iglo))
                         endif
                      enddo

                      do it = jump(ik), 1, -1
                         to_iglo = idx(g_lo, ik, ikx_indexed(it), il, ie, is)
                         if (idx_local (g_lo, to_iglo)) g0(:,:,to_iglo) = 0.
                      enddo

                   enddo
                enddo
             enddo
          endif
       enddo
    end if
  end subroutine exb_shear

  !> Subroutine to setup a redistribute object to be used in enforcing parity
  subroutine init_enforce_parity(redist_obj,ik_ind)
    use theta_grid, only : ntgrid
    use gs2_layouts, only : g_lo,proc_id, ik_idx
    use redistribute, only: index_list_type, init_fill, delete_list, redist_type
    use mp, only: iproc, nproc, max_allreduce
    implicit none
    type(redist_type), intent(out) :: redist_obj
    integer, intent(in), optional :: ik_ind !If present then resulting redistribute object only applies to ik=ik_ind
    integer :: iglo,iglo_conn,ip_conn,ip
    integer :: ig, ndata, nn_max, n,ik_ind_local
    type (index_list_type), dimension(0:nproc-1) :: to, from
    integer, dimension (0:nproc-1) :: nn_from, nn_to
    integer, dimension(3) :: from_low, from_high, to_low, to_high

    !If enforced parity not requested then exit
    if(.not.(def_parity))return
    
    !If not linked then don't need any setup so exit
    if(boundary_option_switch.ne.boundary_option_linked)return
    
    !Deal with optional input
    if(present(ik_ind)) then
       ik_ind_local=ik_ind
    else
       ik_ind_local=-1
    endif

    !NOTE: If we know x is entirely local (g_lo%x_local=.true.) then we
    !know that we don't need to do any comms so our iglo range can be
    !limited to what is on this proc. At the moment we don't take 
    !advantage of this.

    !Use a shorthand for how much data to send per iglo
    ndata=2*ntgrid+1

    !Count how much data is to be sent/recv to/from different procs
    !/First initialise counters
    nn_to =0
    nn_from=0
    
    !/Now loop over iglo to workout what data is to be sent/received
    do iglo=g_lo%llim_world,g_lo%ulim_world
       !Check if we want to skip this
       if(ik_ind_local.gt.0)then
          if(ik_idx(g_lo,iglo).ne.ik_ind_local) cycle
       endif

       !Get proc id of current iglo
       ip=proc_id(g_lo,iglo)

       !Find connected iglo to which we want to send the data
       call get_parity_conn(iglo,iglo_conn,ip_conn)

       !Do we have the data to send?
       if(ip_conn.eq.iproc) nn_to(ip)=nn_to(ip)+ndata

       !Are we going to receive the data?
       if(ip.eq.iproc) nn_from(ip_conn)=nn_from(ip_conn)+ndata

    enddo

    !Now find the maxmimum amount of data to be sent
    nn_max=MAXVAL(nn_to) !Local max
    call max_allreduce(nn_max) !Global max

    !Now define indices of send/receive data
    if(nn_max.gt.0) then
       !Create to/from list objects
       do ip=0,nproc-1
          if(nn_from(ip)>0)then
             allocate(from(ip)%first(nn_from(ip)))
             allocate(from(ip)%second(nn_from(ip)))
             allocate(from(ip)%third(nn_from(ip)))
          endif
          if(nn_to(ip)>0)then
             allocate(to(ip)%first(nn_to(ip)))
             allocate(to(ip)%second(nn_to(ip)))
             allocate(to(ip)%third(nn_to(ip)))
          endif
       enddo

       !Now fill in the lists
       nn_from=0
       nn_to=0

       !Loop over all iglo again (this doesn't scale)
       do iglo=g_lo%llim_world, g_lo%ulim_world
          !Check if we want to skip this
          if(ik_ind_local.gt.0)then
             if(ik_idx(g_lo,iglo).ne.ik_ind_local) cycle
          endif
          
          !Get proc for this iglo
          ip=proc_id(g_lo,iglo)

          !Get connected point
          call get_parity_conn(iglo,iglo_conn,ip_conn)

          !If we're receiving data where do we put it?
          if(ip.eq.iproc)then
             do ig=-ntgrid,ntgrid !Optimised for data send
                n=nn_from(ip_conn)+1
                nn_from(ip_conn)=n
                from(ip_conn)%first(n)=0-ig
                from(ip_conn)%second(n)=1
                from(ip_conn)%third(n)=iglo
             enddo
          endif

          !If we're sending data where do we get it from?
          if(ip_conn.eq.iproc)then
             do ig=-ntgrid,ntgrid !Optimised for data send
                n=nn_to(ip)+1
                nn_to(ip)=n
                to(ip)%first(n)=ig
                to(ip)%second(n)=2
                to(ip)%third(n)=iglo_conn
             enddo
          endif
       enddo

       !Now setup the redistribute object
       from_low(1)=-ntgrid ; from_low(2)=1  ; from_low(3)=g_lo%llim_proc
       from_high(1)=ntgrid ; from_high(2)=2 ; from_high(3)=g_lo%ulim_proc
       to_low(1)=-ntgrid ; to_low(2)=1  ; to_low(3)=g_lo%llim_proc
       to_high(1)=ntgrid ; to_high(2)=2 ; to_high(3)=g_lo%ulim_proc

       !Initialise the fill object
       call init_fill(redist_obj,'c',to_low,to_high,to,&
            from_low,from_high, from)

       !Delete lists
       call delete_list(from)
       call delete_list(to)
    endif
  end subroutine init_enforce_parity

  !> Return the iglo corresponding to the part of the domain given
  !! by iglo reflected in theta=theta0
  subroutine get_parity_conn(iglo,iglo_conn,iproc_conn)
    use gs2_layouts, only: ik_idx,it_idx,proc_id,g_lo
    implicit none
    integer, intent(in) :: iglo
    integer, intent(out) :: iglo_conn, iproc_conn
    integer :: it, ik, it_conn, link, tmp

    !Get indices
    it=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)

    !Initialise to this cell
    tmp=iglo
    
    !Now check number of links
    if(l_links(it, ik).eq.r_links(it, ik))then
       !If in the middle of the domain then iglo doesn't change
       !Just get the proc id
       iproc_conn=proc_id(g_lo,iglo)
       iglo_conn=tmp
    elseif(l_links(it, ik).gt.r_links(it, ik))then
       !If we're on the right then look left
       do link=1,l_links(it, ik)
          !Get the next connected domain
          call get_left_connection(tmp,iglo_conn,iproc_conn)

          !What is the it index here?
          it_conn=it_idx(g_lo,iglo_conn)

          !Update current iglo
          tmp=iglo_conn

          !If the number of right links now matches the left then we've got the match
          if(r_links(it_conn, ik).eq.l_links(it, ik)) exit
       enddo
    else
       !If we're on the left then look right
       do link=1,r_links(it, ik)
          !Get the next connected domain
          call get_right_connection(tmp,iglo_conn,iproc_conn)

          !What is the it index here?
          it_conn=it_idx(g_lo,iglo_conn)

          !Update current iglo
          tmp=iglo_conn

          !If the number of left links now matches the right then we've got the match
          if(l_links(it_conn, ik).eq.r_links(it, ik)) exit
       enddo
    endif
  end subroutine get_parity_conn

  !> Enforce requested parity
  subroutine enforce_parity(redist_obj,ik_ind)
    use theta_grid, only:ntgrid
    use dist_fn_arrays, only: gnew
    use redistribute, only: scatter,redist_type
    use gs2_layouts, only: g_lo,ik_idx
    implicit none
    type(redist_type),intent(in),optional :: redist_obj
    integer, intent(in),optional :: ik_ind
    type(redist_type) :: redist_local
    integer :: ik_local
    integer :: iglo,mult

    !If enforced parity not requested then exit
    if(.not.(def_parity))return
    
    !Set multiplier
    if(even) then
       mult=1
    else
       mult=-1
    endif
    
    !Behaviour depends upon if we're in flux tube or ballooning space
    if(boundary_option_switch.eq.boundary_option_linked) then !Flux-tube
       !Ensure a redist object is present, if not default to parity_redist
       if(present(redist_obj))then
          redist_local=redist_obj
       else
          redist_local=parity_redist
       endif

       !Redistribute the data
       call scatter(redist_local,gnew,gnew)

       !Multiply by factor if required
       if(mult.ne.1) gnew(:,1,:)=mult*gnew(:,1,:)
    else !Ballooning/extended space
       !Ensure ik_local is specified
       if(present(ik_ind))then
          ik_local=ik_ind
       else
          ik_local=-1
       endif

       !Loop over all local iglo
       do iglo=g_lo%llim_proc,g_lo%ulim_alloc
          !Skip if needed
          if(ik_local.gt.0) then
             if(ik_idx(g_lo,iglo).ne.ik_local) cycle
          endif

          !Apply parity filter
          gnew(-ntgrid:-1,1,iglo)=mult*gnew(ntgrid:1:-1,2,iglo)
          gnew(1:ntgrid,1,iglo)=mult*gnew(-1:-ntgrid:-1,2,iglo)
       enddo
    endif
    
  end subroutine enforce_parity

  !> FIXME : Add documentation
  subroutine get_source_term &
       (phi, apar, bpar, phinew, aparnew, bparnew, include_explicit, &
        ab_coefficients, isgn, iglo,ik,it,il,ie,is, sourcefac, source)
#ifdef LOWFLOW
    use lowflow, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    use dist_fn_arrays, only: aj0, aj1, vperp2, vpac, g
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use le_grids, only: lmax, grid_has_trapped_particles, is_ttp, can_be_ttp
    use species, only: spec, nonmaxw_corr
    use run_parameters, only: fphi, fapar, fbpar
    use gs2_time, only: code_dt, wunits
    use hyper, only: D_res
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    logical, intent (in) :: include_explicit
    real, dimension(:), allocatable, intent(in) :: ab_coefficients
    integer, intent (in) :: isgn, iglo, ik, it, il, ie, is
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:), intent (out) :: source

    integer :: ig
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

!CMR, 4/8/2011
! apargavg and phigavg combine to give the GK EM potential chi. 
!          chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! apargavg = apar J0 
! Both quantities are decentred in time and evaluated on || grid points
!
    phigavg  = (fexp(is)*phi(:,it,ik)   + (1.0-fexp(is))*phinew(:,it,ik)) &
                *aj0(:,iglo)*fphi &
             + (fexp(is)*bpar(:,it,ik) + (1.0-fexp(is))*bparnew(:,it,ik))&
                *aj1(:,iglo)*fbpar*2.0*vperp2(:,iglo)*spec(is)%tz
    apargavg = (fexp(is)*apar(:,it,ik)  + (1.0-fexp(is))*aparnew(:,it,ik)) &
                *aj0(:,iglo)*fapar

! source term in finite difference equations
    select case (source_option_switch)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default choice: solve self-consistent equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case (source_option_full)
       if (il <= lmax) then
          call set_source
       else
          source = 0.0
       end if       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve self-consistent terms + include external i omega_d * F_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(source_option_phiext_full)
       if (il <= lmax) then
          call set_source
          if (aky(ik) < epsilon(0.0)) then
             source(:ntgrid-1) = source(:ntgrid-1) &
                  - zi*wdrift(:ntgrid-1,isgn,iglo)*nonmaxw_corr(ie,is) &
                  *2.0*phi_ext*sourcefac*aj0(:ntgrid-1,iglo)
          end if
       else
          source = 0.0
       end if
    case(source_option_homogeneous)
       source = 0.0
    end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do matrix multiplications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (il <= lmax) then
       if (isgn == 1) then
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + b(ig,1,iglo)*g(ig,1,iglo) + a(ig,1,iglo)*g(ig+1,1,iglo)
          end do
       else
          do ig = -ntgrid, ntgrid-1
             source(ig) = source(ig) &
                  + a(ig,2,iglo)*g(ig,2,iglo) + b(ig,2,iglo)*g(ig+1,2,iglo)
          end do
       end if
    end if

!CMR, 21/7/2014: removed redundant line here:    source(ntgrid)=source(-ntgrid) 
!                as source(ntgrid) should never be used.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special source term for totally trapped particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CMR, 13/10/2014: 
! Upper limit of following loops setting source changed from ntgrid to ntgrid-1
! Source is allocated as: source(-ntgrid:ntgrid-1), so ntgrid is out of bounds.

    if (source_option_switch == source_option_full .or. &
        source_option_switch == source_option_phiext_full) then
       if (grid_has_trapped_particles() .and. isgn == 2 .and. can_be_ttp(il)) then
          do ig = -ntgrid, ntgrid-1
             if (.not. is_ttp(ig, il)) cycle
             source(ig) &
                  ! Note we could replace the explicit "2" appearing below with isgn for consistency.
                  = g(ig,2,iglo)*a(ig,2,iglo) &
#ifdef LOWFLOW
                  - zi*(wdttpfac(ig,it,ik,ie,is,2)*hneoc(ig,2,iglo))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*hneoc(ig,2,iglo)*phigavg(ig)
#else
                  - zi*(wdriftttp(ig, 2, iglo))*nonmaxw_corr(ie,is)*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*phigavg(ig)
#endif             
          end do

          if (source_option_switch == source_option_phiext_full .and. &
               aky(ik) < epsilon(0.0)) then
             do ig = -ntgrid, ntgrid-1
                if (.not. is_ttp(ig, il)) cycle
                source(ig) = source(ig) - zi* &
#ifdef LOWFLOW
                     wdttpfac(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#else
                     wdriftttp(ig, isgn, iglo)*nonmaxw_corr(ie,is)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#endif
             end do
          endif

          if (include_explicit) then
             select case (size(ab_coefficients))
             case (1)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*(&
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo))
                end do
             case (2)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*( &
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                        ab_coefficients(2)*gexp_2(ig,isgn,iglo))
                end do
             case (3)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig) = source(ig) + 0.5*code_dt*( &
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                        ab_coefficients(2)*gexp_2(ig,isgn,iglo) + &
                        ab_coefficients(3)*gexp_3(ig,isgn,iglo))
                end do
             end select
          end if
       end if
    end if

  contains
    !> FIXME : Add documentation
    subroutine set_source

      use species, only: spec, nonmaxw_corr
      use theta_grid, only: itor_over_B
      implicit none
      complex :: apar_p, apar_m, phi_p, phi_m!, bpar_p !GGH added bpar_p
      real :: bd, bdfac_p, bdfac_m

! try fixing bkdiff dependence
      bd = bkdiff(1)

      bdfac_p = 1.+bd*(3.-2.*real(isgn))
      bdfac_m = 1.-bd*(3.-2.*real(isgn))



!CMR, 4/8/2011:
! Some concerns, may be red herrings !
! (1) no bakdif factors in phi_m, apar_p, apar_m, vpar !!! 
!                        (RN also spotted this for apar_p)
! (2) can interpolations of products be improved? 
!
!  Attempt at variable documentation:
! phigavg  = phi J0 + 2 T_s/q_s . vperp^2 bpar/bmag J1/Z
! apargavg = apar J0                        (decentered in t) 
! NB apargavg and phigavg combine to give the GK EM potential chi
! chi = phigavg - apargavg*vpa(:,isgn,iglo)*spec(is)%stm
! phi_p = 2 phigavg                      .... (roughly!)
! phi_m = d/dtheta (phigavg)*DTHETA 
! apar_p = 2 apargavg  
! apar_m = 2 d/dt (apar)*DELT  (gets multiplied later by J0 and vpa when included in source)
! => phi_p - apar_p*vpa(:,isgn,iglo)*spec(is)%stm = 2 chi  .... (roughly!)  
! vparterm = -2.0*vpar (IN ABSENCE OF LOWFLOW TERMS)
! wdfac = wdrift + wcoriolis/spec(is)%stm (IN ABSENCE OF LOWFLOW TERMS)
! wstarfac = wstar  (IN ABSENCE OF LOWFLOW TERMS)
! Line below is not true for alphas. vpar = q_s/Tstar * v_|| * the rest. EGH/GW
! vpar = q_s/sqrt{T_s m_s} (v_||^GS2). \gradpar(theta)/DTHETA . DELT (centred)
! wdrift =    q_s/T_s  v_d.\grad_perp . DELT 
! wcoriolis = q_s/T_s  v_C.\grad_perp . DELT 
!
! Definition of source:= 2*code_dt*RHS of GKE
! source     appears to contain following physical terms
!   -2q_s/T_s v||.grad(J0 phi + 2 vperp^2 bpar/bmag J1/Z T_s/q_s).delt 
!   -2d/dt(q v|| J0 apar / T).delt
!   +hyperviscosity
!   -2 v_d.\grad_perp (q J0 phi/T + 2 vperp^2 bpar/bmag J1/Z).delt 
!   -coriolis terms
!   2{\chi,f_{0s}}.delt  (allowing for sheared flow)
!CMRend
!GJW (2018): 
! Ensuring apar respects bakdif. This fixes a numerical instability that prevents 
! nonlinear electromagnetic runs.

      do ig = -ntgrid, ntgrid-1
         phi_p = bdfac_p*phigavg(ig+1)+bdfac_m*phigavg(ig)
         phi_m = phigavg(ig+1)-phigavg(ig)
         ! RN> bdfac factors seem missing for apar_p
         apar_p = bdfac_p*apargavg(ig+1)+bdfac_m*apargavg(ig)
         apar_m = bdfac_p * (aparnew(ig+1,it,ik) - apar(ig+1,it,ik)) + &
              bdfac_m * (aparnew(ig,it,ik) - apar(ig,it,ik))

!MAB, 6/5/2009:
! added the omprimfac source term arising with equilibrium flow shear  
!CMR, 26/11/2010:
! Useful to understand that all source terms propto aky are specified here 
! using Tref=mref vtref^2. See 
! [note by BD and MK on "Microinstabilities in Axisymmetric Configurations"].
! This is converted to  the standard internal gs2 normalisation, 
! Tref=(1/2) mref vtref^2, by wunits, which contains a crucial factor 1/2.
! (Would be less confusing if always used same Tref!)
!
#ifdef LOWFLOW
         source(ig) = (vparterm(ig,isgn,iglo)*phi_m &
              -spec(is)%zstm*vpac(ig,isgn,iglo) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdfac(ig,isgn,iglo)*phi_p) &
              + zi*(wstarfac(ig,isgn,iglo) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
              *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
#else
         source(ig) = (-2.0*vpar(ig,isgn,iglo)*nonmaxw_corr(ie,is)*phi_m &
              -spec(is)%zstm*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is) &
              *((aj0(ig+1,iglo) + aj0(ig,iglo))*0.5*apar_m  &
              + D_res(it,ik)*apar_p) &
              -zi*wdrift(ig,isgn,iglo)*nonmaxw_corr(ie,is)*phi_p) &
              + zi*(wstar(ik,ie,is) &
              + vpac(ig,isgn,iglo)*code_dt*wunits(ik)*ufac(ie,is) &
              -2.0*omprimfac*vpac(ig,isgn,iglo)*nonmaxw_corr(ie,is)*code_dt*wunits(ik)*g_exb*itor_over_B(ig)/spec(is)%stm) &
              *(phi_p - apar_p*spec(is)%stm*vpac(ig,isgn,iglo))
#endif
      end do

      if (include_explicit) then
         select case (size(ab_coefficients))
         case (1)
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*(&
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo))
            end do
         case (2)
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*( &
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                    ab_coefficients(2)*gexp_2(ig,isgn,iglo))
            end do
         case (3)
            do ig = -ntgrid, ntgrid-1
               source(ig) = source(ig) + 0.5*code_dt*( &
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                    ab_coefficients(2)*gexp_2(ig,isgn,iglo) + &
                    ab_coefficients(3)*gexp_3(ig,isgn,iglo))
            end do
         end select
      end if
    end subroutine set_source

  end subroutine get_source_term

  !> This is a version of [[get_source_term]] which does both sign (sigma) together
  !! and uses precalculated constant terms. Leads to more memory usage than 
  !! original version but can be significantly faster (~50%)
  subroutine get_source_term_opt &
       (phi, apar, bpar, phinew, aparnew, bparnew, include_explicit, &
        ab_coefficients, iglo,ik,it,il,ie,is, sourcefac, source)
#ifdef LOWFLOW
    use lowflow, only: hneoc, vparterm, wdfac, wstarfac, wdttpfac
#endif
    use dist_fn_arrays, only: aj0, aj1, vperp2, g
    use dist_fn_arrays, only: gexp_1, gexp_2, gexp_3
    use theta_grid, only: ntgrid
    use kt_grids, only: aky
    use le_grids, only: lmax, grid_has_trapped_particles, is_ttp, can_be_ttp
    use species, only: spec, nonmaxw_corr
    use run_parameters, only: fphi, fapar, fbpar, has_apar, has_bpar, has_phi
    use gs2_time, only: code_dt
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    logical, intent (in) :: include_explicit
    real, dimension(:), allocatable, intent(in) :: ab_coefficients
    integer, intent (in) :: iglo, ik, it, il, ie, is
    complex, intent (in) :: sourcefac
    complex, dimension (-ntgrid:,:), intent (out) :: source
    integer :: ig, isgn
    real :: fac, fac_now, fac_next
    complex, dimension (-ntgrid:ntgrid) :: phigavg, apargavg

    !Temporally weighted fields
    if (has_phi) then
       fac = fphi
       fac_now = fac * fexp(is) ; fac_next = fac * (1 - fexp(is))
       phigavg = (fac_now * phi(:,it,ik) + fac_next * phinew(:,it,ik)) &
                * aj0(:,iglo)
    else
      phigavg = 0.0
    end if

    if (has_bpar) then
       fac = fbpar * 2 * spec(is)%tz
       fac_now = fac * fexp(is) ; fac_next = fac * (1 - fexp(is))
       phigavg = phigavg + &
            (fac_now * bpar(:,it,ik) + fac_next * bparnew(:,it,ik))&
            * aj1(:,iglo) * vperp2(:,iglo)
    end if

    if (has_apar) then
       fac = fapar
       fac_now = fac * fexp(is) ; fac_next = fac * (1 - fexp(is))
       apargavg = (fac_now * apar(:,it,ik)  + fac_next * aparnew(:,it,ik)) &
            * aj0(:,iglo)
    end if

    !Initialise to zero in case where we don't call set_source_opt
    if(il>lmax)then
       source=0.0
    endif

    !Calculate the part of the source related to EM potentials
    !Do both signs at once to improve memory access
    do isgn=1,2
       ! source term in finite difference equations
       select case (source_option_switch)
       case (source_option_full)
          if (il <= lmax) then
             call set_source_opt
          end if
       case(source_option_phiext_full)
          if (il <= lmax) then
             call set_source_opt
             if (aky(ik) < epsilon(0.0)) then
                source(:ntgrid-1,isgn) = source(:ntgrid-1,isgn) &
                  - zi*wdrift(:ntgrid-1,isgn,iglo)*nonmaxw_corr(ie,is) &
                  *2.0*phi_ext*sourcefac*aj0(:ntgrid-1,iglo)
             end if
          end if
       case(source_option_homogeneous)
          source = 0.0
       end select
    enddo
       
    ! Do matrix multiplications
    if (il <= lmax) then
       do ig = -ntgrid, ntgrid-1
          source(ig,1) = source(ig,1) &
               + b(ig,1,iglo)*g(ig,1,iglo) + a(ig,1,iglo)*g(ig+1,1,iglo)
       end do
!CMR, 21/7/2014: removed redundant line here:    source(ntgrid,1)=source(-ntgrid,1) 
!                as source(ntgrid,1) should never be used.
       do ig = -ntgrid, ntgrid-1
          source(ig,2) = source(ig,2) &
               + a(ig,2,iglo)*g(ig,2,iglo) + b(ig,2,iglo)*g(ig+1,2,iglo)
       end do
!CMR, 21/7/2014: removed redundant line here:    source(ntgrid,2)=source(-ntgrid,2) 
!                as source(ntgrid,2) should never be used.
    end if

    ! special source term for totally trapped particles (isgn=2 only)

!CMR, 13/10/2014: 
! Upper limit of following loops setting source changed from ntgrid to ntgrid-1
! Source allocated as: source(-ntgrid:ntgrid-1,2), so ntgrid is out of bounds.

    isgn=2
    if (source_option_switch == source_option_full .or. &
         source_option_switch == source_option_phiext_full) then
       if (grid_has_trapped_particles() .and. can_be_ttp(il)) then
          do ig = -ntgrid, ntgrid-1
             if (.not. is_ttp(ig, il)) cycle
             source(ig,isgn) &
                  = g(ig,2,iglo)*a(ig,2,iglo) &
#ifdef LOWFLOW
                  - zi*(wdttpfac(ig,it,ik,ie,is,2)*hneoc(ig,2,iglo))*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*hneoc(ig,2,iglo)*phigavg(ig)
#else
                  - zi*(wdriftttp(ig, 2, iglo))*nonmaxw_corr(ie,is)*phigavg(ig) &
                  + zi*wstar(ik,ie,is)*phigavg(ig)
#endif             
          end do
             
          if (source_option_switch == source_option_phiext_full .and. &
               aky(ik) < epsilon(0.0)) then
             do ig = -ntgrid, ntgrid-1
                if (.not. is_ttp(ig, il)) cycle
                source(ig,isgn) = source(ig,isgn) - zi* &
#ifdef LOWFLOW
                     wdttpfac(ig,it,ik,ie,is,isgn)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#else
                     wdriftttp(ig, isgn, iglo)*2.0*phi_ext*sourcefac*aj0(ig,iglo)
#endif
             end do
          endif

          ! add in explicit terms
          if (include_explicit) then
             select case (size(ab_coefficients))
             case (1)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*(&
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo))
                end do
             case (2)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*( &
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                        ab_coefficients(2)*gexp_2(ig,isgn,iglo))
                end do
             case (3)
                do ig = -ntgrid, ntgrid-1
                   if (.not. is_ttp(ig, il)) cycle
                   source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*( &
                        ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                        ab_coefficients(2)*gexp_2(ig,isgn,iglo) + &
                        ab_coefficients(3)*gexp_3(ig,isgn,iglo))
                end do
             end select
          end if
       end if
    end if

  contains
    !> FIXME : Add documentation
    subroutine set_source_opt
      implicit none
      complex :: apar_p, apar_m, phi_p, phi_m
      real :: bd, bdfac_p, bdfac_m
      integer :: sgn
! try fixing bkdiff dependence
      bd = bkdiff(1)

      sgn = 3 - 2 * isgn
      bdfac_p = 1 + bd * sgn
      bdfac_m = 1 - bd * sgn

      if (has_apar) then
         do ig = -ntgrid, ntgrid-1
            phi_m = -phigavg(ig) + phigavg(ig+1)
            phi_p = phigavg(ig) * bdfac_m + phigavg(ig+1) * bdfac_p
            apar_m = -apar(ig,it,ik) - apar(ig+1,it,ik) + &
                 aparnew(ig,it,ik) + aparnew(ig+1, it, ik)
            apar_p = apargavg(ig) + apargavg(ig+1)

            source(ig, isgn) = phi_m * source_coeffs(1, ig, isgn, iglo) + &
                 phi_p * source_coeffs(2, ig, isgn, iglo) + &
                 apar_m * source_coeffs(3, ig, isgn, iglo) + &
                 apar_p * source_coeffs(4, ig, isgn, iglo)
         end do
      else
         do ig = -ntgrid, ntgrid-1
            phi_m = -phigavg(ig) + phigavg(ig+1)
            phi_p = phigavg(ig) * bdfac_m + phigavg(ig+1) * bdfac_p
            
            source(ig, isgn) = phi_m * source_coeffs(1, ig, isgn, iglo) + &
                 phi_p * source_coeffs(2, ig, isgn, iglo)
         end do
      end if

      ! add in explicit terms
      if (include_explicit) then
         select case (size(ab_coefficients))
         case (1)
            do ig = -ntgrid, ntgrid-1
               source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*(&
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo))
            end do
         case (2)
            do ig = -ntgrid, ntgrid-1
               source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*( &
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                    ab_coefficients(2)*gexp_2(ig,isgn,iglo))
            end do
         case (3)
            do ig = -ntgrid, ntgrid-1
               source(ig, isgn) = source(ig, isgn) + 0.5*code_dt*( &
                    ab_coefficients(1)*gexp_1(ig,isgn,iglo) + &
                    ab_coefficients(2)*gexp_2(ig,isgn,iglo) + &
                    ab_coefficients(3)*gexp_3(ig,isgn,iglo))
            end do
         end select
      end if
    end subroutine set_source_opt

  end subroutine get_source_term_opt

  !> Calculates the regular source term on the entire local
  !> domain and stores in the output `source`.
  !>
  !> In the future we may wish to move the loop over the domain into the routine
  !> which calculates the source in order to improve memory locality etc.
  subroutine get_source_term_on_local_domain(source, phi, apar, bpar, &
       phinew, aparnew, bparnew, istep, sourcefac)
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use theta_grid, only: ntgrid
    use gs2_time, only: get_adams_bashforth_coefficients
    use nonlinear_terms, only: nonlin, split_nonlinear
    use kt_grids, only: kwork_filter
    implicit none
    complex, dimension (-ntgrid:, :, g_lo%llim_proc:), intent(in out) :: source
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep
    complex, intent (in) :: sourcefac
    logical :: include_explicit
    integer :: iglo, ik, it, il, ie, is, isgn
    real, dimension(:), allocatable :: ab_coefficients

    include_explicit = istep > 0

#ifndef LOWFLOW
    include_explicit = include_explicit .and. nonlin .and. (.not. split_nonlinear)
#endif

    if (include_explicit) allocate(ab_coefficients, &
         source = get_adams_bashforth_coefficients())

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, il, ie, is, isgn) &
    !$OMP SHARED(g_lo, kwork_filter, phi, apar, bpar, phinew, aparnew, bparnew, &
    !$OMP include_explicit, ab_coefficients, sourcefac, source, opt_source) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)

       !Skip work if we're not interested in this ik and it
       if(kwork_filter(it,ik)) cycle

       il = il_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       is = is_idx(g_lo, iglo)

       if(opt_source)then
          call get_source_term_opt (phi, apar, bpar, phinew, aparnew, bparnew, &
               include_explicit, ab_coefficients, iglo, ik, it, il, ie, is, &
               sourcefac, source(:, :, iglo))
       else
          do isgn = 1, 2
             call get_source_term (phi, apar, bpar, phinew, aparnew, bparnew, &
                  include_explicit, ab_coefficients, isgn, iglo, ik, it, il, ie, is, &
                  sourcefac, source(:, isgn, iglo))
          end do
       endif
    end do
    !$OMP END PARALLEL DO
  end subroutine get_source_term_on_local_domain

  !> Calculates the "new" distribution function for the passed source
  !> calculated using the passed fields and current distribution in g.
  !>
  !> Also calculates the homogeneous solution in some cases.
  !>
  !> Note that further work may be required to produce the final distribution function.
  !> In particular with link boundary conditions we may need to add in the homogeneous
  !> solution to ensure continuity across linked boundaries.
  !>
  !> The routine is roughly organised as:
  !>  1. Calculate constants
  !>  2. Loop over domain
  !>     - Calculate flags
  !>     - Set initial conditions
  !>     - Set parallel boundary conditions
  !>     - Solve for homogeneous and inhomogeneous solution
  !>     - Ensure continuity (e.g. periodicity in theta and in sigma at bounce points)
  !>  3. Enforce parity if required.
  subroutine invert_rhs_1 (phinew, bparnew, source)
    use dist_fn_arrays, only: gnew, vperp2, aj1, aj0, set_h_zero, g
    use theta_grid, only: ntgrid
    use le_grids, only: ng2, forbid,&
         is_passing_hybrid_electron, passing_wfb, trapped_wfb, &
         il_is_wfb, is_ttp, il_is_passing, il_is_trapped, &
         is_upper_bounce_point, is_lower_bounce_point, mixed_wfb
    use kt_grids, only: aky
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx
    use run_parameters, only: fbpar, fphi, ieqzip
    use kt_grids, only: kwork_filter
    use species, only: spec, nonmaxw_corr
    use mp, only: mp_abort
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, bparnew
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in out) :: source
    integer :: iglo
    integer :: ig, ik, it, il, ie, is
    integer :: ilmin
    complex :: beta1
    complex :: adjleft, adjright
    logical :: use_pass_homog, speriod_flag
    logical :: is_wfb, is_passing, is_trapped
    logical :: is_self_periodic, is_passing_like

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, il, ie, is, use_pass_homog, speriod_flag, &
    !$OMP ilmin, is_wfb, is_passing, is_trapped, is_self_periodic, &
    !$OMP adjright, adjleft, ig, beta1, is_passing_like) &
    !$OMP SHARED(g_lo, kwork_filter, ieqzip, is_ttp, start_from_previous_solution, &
    !$OMP gnew, phinew, bparnew, fphi, fbpar, boundary_option_switch, g, &
    !$OMP aky, ng2, forbid, g_h, zero_forbid, mixed_wfb, &
    !$OMP connections, trapped_wfb, passing_wfb, wfb, ntgrid, nonad_zero, vperp2, &
    !$OMP aj1, aj0, spec, nonmaxw_corr, l_links, r_links, r, ainv, source) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       !/////////////////////////
       !// Constants, flags etc.
       !/////////////////////////

       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)

       !Skip work if we're not interested in this ik and it
       if(kwork_filter(it,ik)) cycle

       if(ieqzip(it,ik)) cycle

       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)

       !/////////////////////////
       !// Handle hybrid electrons when passing
       !/////////////////////////
       
       if (is_passing_hybrid_electron(is, ik, il))  then
          call set_h_zero (gnew, phinew, bparnew, iglo)
          cycle
       end if

       !CMR, 17/4/2012:
       !  use_pass_homog = T  iff one of following applies:
       !                 boundary_option_self_periodic
       !     OR          boundary_option_linked
       !     OR          aky=0
       !  if use_pass_homog = T, compute homogeneous solution (g1) for passing.
       !        otherwise ONLY compute inhomogenous solution for passing particles.
       !
       !  speriod_flag = T  iff boundary_option_linked AND aky=0
       !  if speriod_flag = T, apply self-periodic bcs for passing and wfb.

       select case (boundary_option_switch)
       case (boundary_option_self_periodic)
          use_pass_homog = .true.
       case (boundary_option_linked)
          use_pass_homog = .true.
          speriod_flag = aky(ik) == 0.0
       case default
          use_pass_homog = .false.
       end select

       use_pass_homog = use_pass_homog .or. aky(ik) == 0.0

       if (use_pass_homog) then
          ilmin = 1
       else
          ilmin = ng2 + 1
       end if

       ! ng2+1 is WFB.
       is_wfb = il_is_wfb(il)

       ! Is this particle passing. Note we exclude the wfb here as we can choose
       ! how we treat this as given by passing_wfb, trapped_wfb etc.
       is_passing = il_is_passing(il)

       ! Should this particle be treated as a passing particle
       is_passing_like = is_passing .or. (is_wfb .and. passing_wfb)

       ! Is this particle trapped? Note we exclude the wfb here as we can choose
       ! how we treat this as given by passing_wfb, trapped_wfb etc.
       is_trapped = il_is_trapped(il)

       ! Decide if this point is self periodic
       !CMR, 17/4/2012:
       ! self_periodic bc is applied as follows:
       ! if boundary_option_linked = .T.
       !      isolated wfb (ie no linked domains)
       !      passing + wfb if aky=0
       ! else (ie boundary_option_linked = .F.)
       !      wfb
       !      passing and wfb particles if boundary_option_self_periodic = T
       !
       !CMR, 6/8/2014:
       ! Parallel BC for wfb is as follows:
       !    ballooning space
       !           wfb ALWAYS self-periodic (-ntgrid:ntgrid)
       !    linked fluxtube
       !           wfb self-periodic (-ntgrid:ntgrid) ONLY if ISOLATED,
       !                                               OR if iky=0
       !           OTHERWISE wfb treated as passing for now in (-ntgrid:ntgrid)
       !                     store homogeneous and inhomog solutions
       !                     AND set amplitude of homog later in invert_rhs_linked
       !                     to satisfy self-periodic bc in fully linked SUPER-CELL
       if (boundary_option_switch == boundary_option_linked) then
          is_self_periodic = (speriod_flag .and. is_passing) & ! Passing particle with ky=0
               .or. (speriod_flag .and. is_wfb .and. .not. trapped_wfb) & ! Non-trapped wfb treatment of wfb with ky = 0
               .or. (is_wfb .and. .not. connections(iglo)%neighbor & ! Isolated (i.e. ky/=0 cell with no connections) wfb
               .and. mixed_wfb )     ! treated as neither passing or trapped particle (i.e. unique)
       else
          is_self_periodic = (use_pass_homog .and. is_passing_like) & ! ky=0 or with self-periodic boundaries and passing like particles
               .or. (is_wfb .and. mixed_wfb) ! Mixed wfb is always self periodic when not linked
       end if

       !/////////////////////////
       !// Initialisation
       !/////////////////////////
       if (start_from_previous_solution) then
          ! Initialise gnew to g. Note we probably only need to do this when
          ! istep == istep_last (i.e. on the second half of a timestep). We
          ! may also be able to get a slight improvement by copying over the
          ! non-adiabatic part of g and then forming gnew using the new fields.
          ! For example we could do the equivalent of
          ! call g_adjust(gnew, phi, bpar, from_g_gs2)
          ! call g_adjust(gnew, phinew, bparnew, to_g_gs2)
          ! after copying g into gnew.
          gnew(:, :, iglo) = g(:, :, iglo)
       else
          gnew(:, :, iglo) = 0.0
       end if

       !/////////////////////////
       !// Boundaries
       !/////////////////////////

       !CMR, 18/4/2012:
       ! What follows is a selectable improved parallel bc for passing particles.
       !                                            (prompted by Greg Hammett)
       ! Original bc is: g_gs2 = gnew = 0 at ends of domain:
       !   ONLY implies zero incoming particles in nonadiabatic delta f if phi=bpar=0
       ! Here ensure ZERO incoming particles in nonadiabatic delta f at domain ends
       !  (exploits code used in subroutine g_adjust to transform g_wesson to g_gs2)
       if ( nonad_zero .and. is_passing_like ) then
          !Only apply the new boundary condition to the leftmost
          !cell for sign going from left to right
          if (l_links(it, ik) .eq. 0) then
             adjleft = 2.0*vperp2(-ntgrid,iglo)*aj1(-ntgrid,iglo) &
                  *bparnew(-ntgrid,it,ik)*fbpar &
                  + spec(is)%zt*phinew(-ntgrid,it,ik)*aj0(-ntgrid,iglo) &
                  *fphi*nonmaxw_corr(ie,is)
             gnew(-ntgrid,1,iglo) = -adjleft
          end if
          !Only apply the new boundary condition to the rightmost
          !cell for sign going from right to left
          if (r_links(it, ik) .eq. 0) then
             adjright = 2.0*vperp2(ntgrid,iglo)*aj1(ntgrid,iglo) &
                  *bparnew(ntgrid,it,ik)*fbpar &
                  + spec(is)%zt*phinew(ntgrid,it,ik)*aj0(ntgrid,iglo) &
                  *fphi*nonmaxw_corr(ie,is)
             gnew(ntgrid,2,iglo) = -adjright
          end if
       end if

       !/////////////////////////
       !// Time advance
       !/////////////////////////

       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! time advance vpar < 0   !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! inhomogeneous part: gnew
       ! r=ainv=0 if forbid(ig,il) or forbid(ig+1,il), so gnew=0 in forbidden
       ! region and at upper bounce point
       do ig = ntgrid-1, -ntgrid, -1
          gnew(ig,2,iglo) = -gnew(ig+1,2,iglo)*r(ig,2,iglo) + ainv(ig,2,iglo)*source(ig,2,iglo)
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! time advance vpar > 0   !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! First set BCs for trapped particles at lower bounce point
       if (is_trapped) then
          ! match boundary conditions at lower bounce point
          ! Note we can almost change this range to the full theta domain
          ! now we use is_lower_bounce_point. We would require a slight
          ! change in how we handle copying the bounce point value into the source
          ! as currently we set it at one point before the bounce point. We might
          ! be able to set gnew at the bounce point directly with some small
          ! changes to ainv etc.
          do ig = -ntgrid, ntgrid - 1
             ! Skip for ttp
             if (is_ttp(ig, il)) cycle
             if (is_lower_bounce_point(ig+1, il)) then
                !CMR: init_invert_rhs  sets ainv=1 at lower bounce point for trapped
                !     source below ensures gnew(lb,1,iglo)=gnew(lb,2,iglo)
                !     where lb is the lower bounce point.
                source(ig,1,iglo) = gnew(ig+1,2,iglo)
             end if
          end do
       end if

       ! If trapped wfb enforce the bounce condition at the lower bounce point
       ! Note this treats the wfb as bouncing at the end of the parallel domain.
       ! For nperiod > 1 this means the wfb does not bounce at the interior points
       ! where bmag = bmax
       if (is_wfb .and. trapped_wfb) then
          ig = -ntgrid
          gnew(ig,1,iglo) = gnew(ig,2,iglo)
       endif

       ! time advance vpar > 0 inhomogeneous part
       do ig = -ntgrid, ntgrid - 1
          if (.not. is_ttp(ig, il)) then
             ! Use two statements here to *prevent* erroneous vectorization of this loop by gfortran 7.2
             gnew(ig + 1, 1, iglo) = ainv(ig, 1, iglo) * source(ig, 1, iglo)
             gnew(ig + 1, 1, iglo) = gnew(ig + 1, 1, iglo) - gnew(ig, 1, iglo) * r(ig, 1, iglo)
          else
             ! If we're totally trapped then just enforce continuity in sigma.
             ! Note we don't allow for totally trapped particles at ig = +/- ntgrid
             ! so it doesn't matter that the outer loop does not cover the full range.
             if (.not. forbid(ig, il)) gnew(ig, 1, iglo) = gnew(ig, 2, iglo)
          end if
       end do

       !/////////////////////////
       !// Periodicity at bounce point and parallel boundary
       !/////////////////////////

       if (is_self_periodic) call self_periodic(gnew, g_h, iglo, is, ie, it, ik, is_wfb, phinew, bparnew)

       ! add correct amount of homogeneous solution for trapped particles to satisfy boundary conditions
       !CMR, 24/7/2014:
       !Revert to looping from il>= ng2+2, i.e. exclude wfb as:
       !          (1) wfb bc is already handled above
       !          (2) forbid never true for wfb, so including ng2+1 in loop is pointless.
       if (is_trapped) then
          beta1 = 0.0
          do ig = ntgrid-1, -ntgrid, -1
             ! If totally trapped or forbidden then cycle
             if (forbid(ig, il) .or. is_ttp(ig, il)) then
                beta1 = 0.0
                cycle !CMR: to avoid pointless arithmetic later in loop
             else if (forbid(ig+1,il)) then !This could be is_upper_bounce_point(ig, il)?
                beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g_h(ig, 1, iglo))
             end if
             gnew(ig,1,iglo) = gnew(ig,1,iglo) + beta1*g_h(ig, 1, iglo)
             gnew(ig,2,iglo) = gnew(ig,2,iglo) + beta1*g_h(ig, 2, iglo)
          end do
       end if

       ! If wfb is trapped, enforce the bounce condition at the upper bounce point
       ! by combining homogeneous and inhomogeneous solutions
       ! Note this treats the wfb as bouncing at the end of the parallel grid.
       ! When nperiod  > 1 this means we do not consider the wfb to bounce at the
       ! interior points with bmag=bmax. This treatment works out slightly similar
       ! to the self_periodic treatment but mixes the different sigma.
       if (is_wfb .and. trapped_wfb) then
          ig = ntgrid
          ! Note here we find the difference between gnew at different sigma
          ! but the same theta, whilst in self_periodic we find the difference
          ! between different theta at the same sigma.
          ! Note if trapped_wfb then g1(ntgrid, 2) == 1.0 so the denominator
          ! is equivalent to g1(ig,2) - g1(ig,1) here.
          beta1 = (gnew(ig,1,iglo) - gnew(ig,2,iglo))/(1.0 - g_h(ig, 1, iglo))
          gnew(:, :, iglo) = gnew(:, :, iglo) + beta1*g_h(:, :, iglo)
       end if

       !CMR,DD, 25/7/2014:
       ! Not keen on following kludge zeroing forbidden region
       ! Introduced new flag: zero_forbid defaults to .false.
       !  Tested and default is fine linearly, expect should work nonlinearly,
       !  (Can easily restore old behaviour by setting: zero_forbid=.true.)
       ! zero out spurious gnew outside trapped boundary
       if(zero_forbid)then
          where (forbid(:,il))
             gnew(:,1,iglo) = 0.0
             gnew(:,2,iglo) = 0.0
          end where
       endif

    end do
    !$OMP END PARALLEL DO

    if (def_parity) call enforce_parity(parity_redist)
  end subroutine invert_rhs_1

  !> Sets sum of homogeneous and inhomogeneous solutions to give a result
  !>     gnew(ntgrid,2) = gnew(-ntgrid,2)
  !>     gnew(ntgrid,1) = gnew(-ntgrid,1)
  !> ie periodic bcs at the ends of the flux tube.
  !>
  !> CMR, 25/7/2014:
  !> self-periodic applied to zonal modes (makes sense)
  !>                      and wfb        (seems strange)
  !> adding adjr, adjl to cope with application of self-periodic to WFB
  !>   dadj=adjl-adjr will be ZERO for ky=0 modes, but NOT for WFB.
  !> This change sets g_wesson (or h) to be self-periodic for wfb, not g !!!
  !> NB this code change will implement this only in ballooning space,
  !> and not in a linked fluxtube.
  subroutine self_periodic(gnew, g_h, iglo, is, ie, it, ik, is_wfb, phinew, bparnew)
    use species, only: spec, nonmaxw_corr
    use dist_fn_arrays, only: vperp2, aj1, aj0
    use run_parameters, only: fphi, fbpar
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use le_grids, only: passing_wfb
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in out) :: gnew
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_h
    integer, intent(in) :: iglo, is, ie, it, ik
    logical, intent(in) :: is_wfb
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, bparnew
    complex :: adjl, adjr, dadj, beta1

    if ((is_wfb .and. nonad_zero) .and. .not. passing_wfb) then
       adjl = 2.0*vperp2(-ntgrid,iglo)*aj1(-ntgrid,iglo) &
            *bparnew(-ntgrid,it,ik)*fbpar &
            + spec(is)%zt*phinew(-ntgrid,it,ik)*aj0(-ntgrid,iglo) &
            *nonmaxw_corr(ie,is)*fphi
       adjr = 2.0*vperp2(ntgrid,iglo)*aj1(ntgrid,iglo) &
            *bparnew(ntgrid,it,ik)*fbpar &
            + spec(is)%zt*phinew(ntgrid,it,ik)*aj0(ntgrid,iglo) &
            *nonmaxw_corr(ie,is)*fphi
       dadj = adjl-adjr
    else
       dadj = cmplx(0.0,0.0)
    endif

    if (g_h(ntgrid, 1, iglo) /= 1.) then
       beta1 = (gnew(ntgrid,1,iglo) - gnew(-ntgrid,1,iglo) - dadj)/(1.0 - g_h(ntgrid, 1, iglo))
       gnew(:,1,iglo) = gnew(:,1,iglo) + beta1*g_h(:, 1, iglo)
    end if

    if (g_h(-ntgrid, 2, iglo) /= 1.) then
       beta1 = (gnew(-ntgrid,2,iglo) - gnew(ntgrid,2,iglo) + dadj)/(1.0 - g_h(-ntgrid, 2, iglo))
       gnew(:,2,iglo) = gnew(:,2,iglo) + beta1*g_h(:, 2, iglo)
    end if

    ! Should we abort if either of the /= 1. checks return false? In this case we end up doing nothing
    ! so we might want to at least check that gnew _is_ self-periodic already and if not abort as we've
    ! asked to make this point self-periodic and we've not been able to.
  end subroutine self_periodic

  !> A shorthand for [[apply_linked_boundary_conditions]] passing the regular
  !> module level instances of each of the required arguments. Primarily
  !> intended for testing to avoid having to expose additional data arrays.
  subroutine enforce_linked_boundary_conditions
    use fields_arrays, only: phinew, bparnew
    use dist_fn_arrays, only: gnew
    implicit none
    call apply_linked_boundary_conditions(gnew, g_h, g_adj, phinew, bparnew)
  end subroutine enforce_linked_boundary_conditions

  !> Applies linked (twist-shift) boundary conditions to passed distribution function
  subroutine apply_linked_boundary_conditions (g_in, g_h, g_adj, phinew, bparnew)
    use dist_fn_arrays, only: g_work
    use theta_grid, only: ntgrid, bmag
    use le_grids, only: is_passing_hybrid_electron, il_is_wfb, mixed_wfb, energy, al
    use gs2_layouts, only: g_lo, ik_idx, it_idx, il_idx, ie_idx, is_idx, idx
    use redistribute, only: fill
    use run_parameters, only: fbpar, fphi
    use species, only: spec, nonmaxw_corr
    use array_utils, only: zero_array
    use kt_grids, only: kwork_filter, kperp2
    use spfunc, only: j0, j1
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in out) :: g_in, g_h
    complex, dimension(:, :, g_lo%llim_proc:), intent(in out) :: g_adj
    complex, dimension (-ntgrid:, :, :), intent(in) :: phinew, bparnew
    complex :: b0, fac, facd
    real :: vperp2left, vperp2right, argl, argr, aj0l, aj0r, aj1l, aj1r
    integer :: il, ik, it, n, i
    integer :: iglo, ncell
!
! adding adjr, adjl to cope with application of self-periodic to WFB
!   dadj=adjl-adjr will be ZERO for ky=0 modes, but NOT for WFB.
! This change sets g_wesson (or h) to be self-periodic for wfb, not g !!!
! NB this code change implement this in a linked fluxtube.
!
    integer :: ie, is, itl, itr
    complex :: adjl, adjr, dadj

    ! No more work to do if no connections so return
    if (no_connections) return

    ! Some extra work to do if we allow non-zero incoming boundary conditions
    if (start_from_previous_solution) then
       ! Pass the incoming boundary values to the previous connected cell, storing
       ! at its connected boundary. This requires two redistribute calls _and_ a full
       ! dist_fn sized array even though there should be at most two non-zero values
       ! per iglo point. Future work should be able to make this much more efficient
       call zero_array(g_work)
       call fill(incoming_links, g_in, g_work)

       ! Construct g_work = g_in - g_work where g_work initially contains
       ! the incoming boundary values from the upstream connected
       ! cell. This adapts the linked boundary conditions to allow for
       ! non-zero incoming boundary conditions (gnew /= 0 ) for internal
       ! (i.e. connected cells)
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo) &
       !$OMP SHARED(g_lo, g_work, g_in, ntgrid) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(ntgrid, 1, iglo) = g_in(ntgrid, 1, iglo) - g_work(ntgrid, 1, iglo)
          g_work(-ntgrid, 2, iglo) = g_in(-ntgrid, 2, iglo) - g_work(-ntgrid, 2, iglo)
       end do
       !$OMP END PARALLEL DO

       !<DD>Note these fill routines are often equivalent to an all-to-all type
       !communication, i.e. when nproc --> 2 nproc, t_fill --> 4 t_fill
       !By only communicating with our direct neighbours we would significantly
       !reduce the amount of data to communicate and we should improve the communication
       !scaling. However, if we do this then we lose the ability to perform the linked
       !update (i.e. what we do below) in a parallel manner, so the below code would
       !become partially serial.
       call fill (links_p, g_work, g_adj)
       if( mixed_wfb ) call fill (wfb_p, g_work, g_adj)
    else
       call fill (links_p, g_in, g_adj)
       if( mixed_wfb ) call fill (wfb_p, g_in, g_adj)
    end if

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, ik, it, ncell, il, is, ie, itl, itr, &
    !$OMP vperp2left, vperp2right, argl, argr, aj0l, aj0r, aj1l, aj1r, &
    !$OMP adjl, adjr, dadj, facd, b0, fac, i, n) &
    !$OMP SHARED(g_lo, r_links, l_links, mixed_wfb, nonad_zero, kperp2, &
    !$OMP ntgrid, spec, bparnew, fbpar, nonmaxw_corr, kwork_filter, al, &
    !$OMP fphi, phinew, save_h, g_adj, g_in, g_h, n_links, bmag, energy) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       if (.not. any(save_h(:, iglo))) cycle
       ik = ik_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       if (kwork_filter(it, ik)) cycle
       ncell = r_links(it, ik) + l_links(it, ik) + 1
       if (ncell == 1) cycle
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       if (is_passing_hybrid_electron(is, ik, il))  cycle
       ! MRH above line causes us to skip the loops below which
       ! build the solution along the connected cells
       ! MRH we skip this in the case that we have taken the non-zonal
       ! passing electron (il<=ng2+1 including wfb) solution to be h=0

       ! wfb
       if ( il_is_wfb(il) .and. mixed_wfb ) then
          if (nonad_zero) then
             ie = ie_idx(g_lo,iglo)

             itl=get_leftmost_it(it,ik)
             itr=get_rightmost_it(it,ik)

             vperp2left = bmag(-ntgrid)*al(il)*energy(ie,is)
             vperp2right = bmag(ntgrid)*al(il)*energy(ie,is)

             argl = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(-ntgrid)*kperp2(-ntgrid,itl,ik))
             argr = spec(is)%bess_fac*spec(is)%smz*sqrt(energy(ie,is)*al(il)/bmag(ntgrid)*kperp2(ntgrid,itr,ik))

             aj0l = j0(argl)
             aj0r = j0(argr)
             aj1l = j1(argl)
             aj1r = j1(argr)

             !CMR, 29/8/2014:
             !(1) compute adjl, adjr: the corrections in mapping from g_gs2 to g_wesson
             !                     at the extreme left and right of the supercell
             !(2) dadj=adjl-adjr then used to apply the self-periodic bc to g_wesson
             !    (was previously applied to g)
             !
             adjl = 2.0*vperp2left*aj1l &
                  *bparnew(-ntgrid,itl,ik)*fbpar &
                  + spec(is)%zt*phinew(-ntgrid,itl,ik)*aj0l &
                  *nonmaxw_corr(ie,is)*fphi
             adjr = 2.0*vperp2right*aj1r &
                  *bparnew(ntgrid,itr,ik)*fbpar &
                  + spec(is)%zt*phinew(ntgrid,itr,ik)*aj0r &
                  *nonmaxw_corr(ie,is)*fphi
             dadj = adjl-adjr
          else
             dadj = 0.0
          end if

          if (save_h(1, iglo)) then
             !CMR, 7/8/2014:
             ! This code implements "self-periodic" parallel BC for wfb
             ! g_adj contains exit point inhomog and homog sol'ns from ALL wfb cells
             ! g_adj(j,1,iglo):    bcs for rightwards travelling particles (RP)
             !  j=1:ncell : inhom sol'n at ntgrid in cell j from RIGHT for RP
             !  j=ncell+1:2ncell : hom sol'n at ntgrid in cell (2ncell+1-j) from RIGHT for RP
             ! g_adj(j,2,iglo):    bcs for leftwards travelling particles (LP)
             !  j=1:ncell : inhom sol'n at -ntgrid in cell j from LEFT for LP
             !  j=ncell+1:2ncell : hom sol'n at -ntgrid in cell (2ncell+1-j) from LEFT for LP
             !
             !  facd= 1/(1-\Prod_{j=1,ncell} h^r_j)  (see CMR Parallel BC note)
             !    where h^r_j is the homogeneous exit solution from cell j for RP
             !
             facd = 1./(1. - product(g_adj(ncell + 1 : 2 * ncell, 1, iglo)))

             b0 = 0.

             ! i labels cell counting from LEFTmost cell.
             fac = 1
             do i = ncell - 1, 1, -1
                ! fac is product of all homog sol's from cells to RIGHT on cell i
                ! g_adj(ncell+j,1,iglo) accesses these homog solutions
                fac = fac * g_adj(ncell + i + 1, 1, iglo)
                ! g_adj(ncell+1-i,1,iglo) accesses inhom solution from cell i
                b0 = b0 + fac * g_adj(ncell+1-i,1,iglo)
             end do

             ! b0 computed next line is homog amplitude in leftmost cell  (see CMR note)
             b0 = (b0 + g_adj(1,1,iglo)-dadj)*facd

             ! BUT we NEED homog amplitude in THIS cell.
             ! Solve matrix BC equation by cascading homog amplitude, b0, rightwards
             !        from leftmost cell to THIS cell.

             do i = 1, l_links(it, ik)
                !  Loop implements cascade to right, using CMR note equ'n:
                !           a^r_{j+1} = a^r_j h^r_j + i^r_j

                b0 = b0 * g_adj(ncell+i,1,iglo) + g_adj(ncell+1-i,1,iglo)
             end do

             ! Resultant b0 is homog amplitude for THIS cell.
             ! Finally add correct homog amplitude to gnew to match the parallel BC.
             g_in(:,1,iglo) = g_in(:,1,iglo) + b0*g_h(:,1,iglo)
          endif

          ! Following code repeats same calculation for LEFTWARD travelling wfb.
          !CMRend
          if (save_h(2, iglo)) then
             facd = 1./(1. - product(g_adj(ncell + 1 : 2 * ncell, 2, iglo)))

             b0 = 0.
             fac = 1.0
             do i = ncell - 1, 1, -1
                fac = fac * g_adj(ncell + i + 1, 2, iglo)
                b0 = b0 + fac * g_adj(ncell+1-i,2,iglo)
             end do
             b0 = (b0 + g_adj(1,2,iglo)+dadj)*facd

             do i = 1, r_links(it, ik)
                b0 = b0 * g_adj(ncell+i,2,iglo) + g_adj(ncell+1-i,2,iglo)
             end do

             g_in(:,2,iglo) = g_in(:,2,iglo) + b0*g_h(:,2,iglo)
          end if
       else
          !
          ! n_links is the number of complex numbers required to fix the boundary
          ! conditions in each cell that is a member of a supercell with at least two
          ! cells, and for which the bounce point is not at theta=pi.
          !
          if (save_h(1, iglo) .and. l_links(it, ik) > 0) then
             ! The l_links test might already be part of save_h?
             n = n_links(1, it, ik)
             b0 = g_adj(1, 1, iglo)
             fac = 1.0
             do i = 2, l_links(it, ik)
                fac = fac * g_adj(n + 2 - i, 1, iglo)
                b0 = b0 + g_adj(i, 1, iglo) * fac
             end do
             g_in(:, 1, iglo) = g_in(:, 1, iglo) + b0 * g_h(:, 1, iglo)
          end if

          if (save_h(2, iglo) .and. r_links(it, ik) > 0) then
             ! The r_links test might already be part of save_h?
             n = n_links(2, it, ik)
             b0 = g_adj(1, 2, iglo)
             fac = 1.0
             do i = 2, r_links(it, ik)
                fac = fac * g_adj(n + 2 - i, 2, iglo)
                b0 = b0 +  g_adj(i, 2, iglo) * fac
             end do
             g_in(:, 2, iglo) = g_in(:, 2, iglo) + b0 * g_h(:, 2, iglo)
          end if
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine apply_linked_boundary_conditions

  !> FIXME : Add documentation
  subroutine invert_rhs (phi, apar, bpar, phinew, aparnew, bparnew, istep)
    use theta_grid, only: ntgrid
    use gs2_time, only: code_time
    use constants, only: zi, pi
    use dist_fn_arrays, only: g_work, gnew
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    integer, intent (in) :: istep

    real :: time
    complex :: sourcefac

    time = code_time
    !Sourcefac ends up being passed all the way through to get_source_term
    !where it is multiplied by phi_ext (default 0.0) and added to ky<epsilon(0.0)
    !modes source term. Should probably just be calculated in get_source_term and
    !only if min(ky)<epsilon(0.0) & phi_ext/=0 & source_option_switch==source_option_phiext_full
    !
    if (time > t0) then
       sourcefac = source0*exp(-zi*omega0*time+gamma0*time)
    else
       sourcefac = (0.5 - 0.5*cos(pi*time/t0))*exp(-zi*omega0*time+gamma0*time)
    end if

    ! Is this needed?
    call zero_array(g_work)

    call get_source_term_on_local_domain(g_work, phi, apar, bpar, phinew, &
         aparnew, bparnew, istep, sourcefac)

    call invert_rhs_1 (phinew, bparnew, source = g_work)

    select case (boundary_option_switch)
    case (boundary_option_linked)
       ! Note if no_connections default was true we could also unconditionally call apply_linked_boundary_conditions
       ! or just call the below if no_connections is false.
       call apply_linked_boundary_conditions(gnew, g_h, g_adj, phinew, bparnew)
    end select

  end subroutine invert_rhs

  !> Ensure that linked boundary values of passed complex field are
  !> single valued (e.g. kperp2(ntgrid,ikx,iky) is
  !> equal to kperp2(-ntgrid,ikx_link,iky) as these correspond to the same location).
  !>
  !> @note itright(it, ik) == it for ky == 0 modes in the flux tube setup so this
  !> routine will also copy one periodic boundary value into the other one.
  subroutine ensure_single_val_fields_pass(arr, force_all)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use optionals, only: get_option_with_default
    use kt_grids, only: naky, ntheta0
    implicit none
    integer :: it, ik, link_it, ik_min, ik_max, it_min, it_max
    complex, dimension (-ntgrid:,:,:), intent (in out) :: arr
    logical, intent(in), optional :: force_all
    logical :: force
    if (boundary_option_switch /= boundary_option_linked) return
    force = get_option_with_default(force_all, .false.)
    if (force) then
       ik_min = 1; ik_max = naky ; it_min = 1 ; it_max = ntheta0
    else
       ik_min = g_lo%ik_min; ik_max = g_lo%ik_max
       it_min = g_lo%it_min; it_max = g_lo%it_max
    end if

    do ik = ik_min, ik_max
       do it = it_min, it_max
          link_it = itright(it, ik)
          if (link_it < 0) cycle
          arr(-ntgrid, link_it, ik) = arr(ntgrid, it, ik)
       end do
    end do
  end subroutine ensure_single_val_fields_pass

  !> Ensure that linked boundary values of passed real field are
  !> single valued (e.g. kperp2(ntgrid,ikx,iky) is
  !> equal to kperp2(-ntgrid,ikx_link,iky) as these correspond to the same location).
  !>
  !> @note itright(it, ik) == it for ky == 0 modes in the flux tube setup so this
  !> routine will also copy one periodic boundary value into the other one.
  subroutine ensure_single_val_fields_pass_r(arr, force_all)
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo
    use optionals, only: get_option_with_default
    use kt_grids, only: naky, ntheta0
    implicit none
    integer :: it, ik, link_it, ik_min, ik_max, it_min, it_max
    real, dimension (-ntgrid:,:,:), intent (in out) :: arr
    logical, intent(in), optional :: force_all
    logical :: force
    if (boundary_option_switch /= boundary_option_linked) return
    force = get_option_with_default(force_all, .false.)
    if (force) then
       ik_min = 1; ik_max = naky ; it_min = 1 ; it_max = ntheta0
    else
       ik_min = g_lo%ik_min; ik_max = g_lo%ik_max
       it_min = g_lo%it_min; it_max = g_lo%it_max
    end if

    do ik = ik_min, ik_max
       do it = it_min, it_max
          link_it = itright(it, ik)
          if (link_it < 0) cycle
          arr(-ntgrid, link_it, ik) = arr(ntgrid, it, ik)
       end do
    end do
  end subroutine ensure_single_val_fields_pass_r

  !> Compute velocity space integrals over \(h_s\):
  !> $$\mathrm{antot} = n_s q_s \int d^3 v f_0 J_0(Z_s) h_s$$
  !> $$\mathrm{antota} = 2 \beta n_s q_s \sqrt{\tfrac{T_s}{m_s}} \int d^3 v f_0 J_0(Z_s) h_s$$
  !> $$\mathrm{antotp} = n_s T_s \int d^3 v f_0 \frac{v_\perp^2 J_1(Z_s)}{Z_s} h_s$$
  !>
  !> See Colin's field equation notes for more details
  !>
  !> Takes `g_in` as \(h_s\)
  !>
  !> FIXME: Get these notes on the website
  subroutine getan_from_dfn (g_in, antot, antota, antotp, local_only)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, g_work
    use kt_grids, only: kperp2, naky, ntheta0
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species
    use run_parameters, only: beta, has_phi, has_apar, has_bpar
    use gs2_layouts, only: g_lo
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in out) :: antot, antota, antotp
    logical, optional :: local_only
    real, dimension (nspec) :: wgt
    integer :: isgn, iglo, ig
    logical :: local_local_only

    local_local_only = get_option_with_default(local_only, .false.)

    !Don't do this as integrate_species will fill in all values
    !    antot=0. ; antota=0. ; antotp=0.
    !Need to set individual arrays to zero if not using integrate_species for
    !that field. (NOTE this probably isn't actually needed as we shouldn't
    !use the various antots if we're not calculating/using the related field).

    if (has_phi) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, g_work, aj0, g_in) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g_work(ig,isgn,iglo) = aj0(ig,iglo)*g_in(ig,isgn,iglo)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       wgt = spec%z*spec%dens
       call integrate_species (g_work, wgt, antot)

       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)
       if (esv) call ensure_single_val_fields_pass(antot, force_all = .true.)
    end if

    if (has_apar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, g_work, aj0, vpa, g_in) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g_work(ig,isgn,iglo) = aj0(ig,iglo)*vpa(ig,isgn,iglo)*g_in(ig,isgn,iglo)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       
       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g_work, wgt, antota)
       if (esv) call ensure_single_val_fields_pass(antota, force_all = .true.)
    end if

    if (has_bpar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(iglo, isgn, ig) &
       !$OMP SHARED(g_lo, ntgrid, g_work, aj1, vperp2, g_in) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g_work(ig,isgn,iglo) = aj1(ig,iglo)*vperp2(ig,iglo)*g_in(ig,isgn,iglo)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       wgt = spec%temp*spec%dens
       call integrate_species (g_work, wgt, antotp)
       if (esv) call ensure_single_val_fields_pass(antotp, force_all = .true.)
    end if

  end subroutine getan_from_dfn

  !> Compute velocity space integrals over \(h_s\):
  !> $$\mathrm{antot} = n_s q_s \int d^3 v f_0 J_0(Z_s) h_s$$
  !> $$\mathrm{antota} = 2 \beta n_s q_s \sqrt{\tfrac{T_s}{m_s}} \int d^3 v f_0 J_0(Z_s) h_s$$
  !> $$\mathrm{antotp} = n_s T_s \int d^3 v f_0 \frac{v_\perp^2 J_1(Z_s)}{Z_s} h_s$$
  !>
  !> See Colin's field equation notes for more details
  !>
  !> Always uses `gnew` as \(h_s\)
  !>
  !> FIXME: Get these notes on the website
  subroutine getan (antot, antota, antotp, local_only)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: antot, antota, antotp
    logical, intent(in), optional :: local_only
    call getan_from_dfn(gnew, antot, antota, antotp, local_only)
  end subroutine getan


  !> Check if the passed distribution function shaped array has the same
  !> value at the duplicate theta points arising from the linked boundary conditions.
  !>
  !> This is primarily intended as a utility to aid debugging.
  subroutine check_linked_boundaries_are_satisfied(g_array, pass_right, pass_left, &
       identifier, tolerance)
    use mp, only: mp_abort, max_allreduce
    use theta_grid, only: ntgrid
    use gs2_layouts, only: g_lo, il_idx, it_idx, ik_idx
    use le_grids, only: il_is_trapped, il_is_wfb, trapped_wfb
    use optionals, only: get_option_with_default
    use redistribute, only: redist_type, fill
    use array_utils, only: zero_array
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_array
    type(redist_type), intent(in out) :: pass_right, pass_left
    character(len=*), intent(in), optional :: identifier
    real, intent(in), optional :: tolerance
    complex, dimension(:, :, :), allocatable :: g_pass
    character(len=:), allocatable :: id_message
    integer :: iglo, il, it, ik
    real :: violation_flag, local_tolerance, err

    ! Exit if no connections as nothing to check
    if (no_connections) return
    if (boundary_option_switch /= boundary_option_linked) return

    local_tolerance = get_option_with_default(tolerance, epsilon(0.0) * 1.0e2)
    id_message = get_option_with_default(identifier, '')
    allocate(g_pass(-ntgrid:ntgrid, 2, g_lo%llim_proc:g_lo%ulim_alloc))

    err = 0.0

    call zero_array(g_pass)
    call fill(pass_right, g_array, g_pass)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, it, ik) &
    !$OMP SHARED(g_lo, trapped_wfb, g_array, g_pass, l_links, ntgrid) &
    !$OMP REDUCTION(max: err) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ! Don't check for trapped particles
       if (il_is_trapped(il)) cycle
       if (il_is_wfb(il) .and. trapped_wfb) cycle

       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)

       if (l_links(it, ik) > 0) then
          err = max(err, abs(g_array(-ntgrid, 1, iglo) - g_pass(-ntgrid, 1, iglo)))
       end if
    end do
    !$OMP END PARALLEL DO

    call zero_array(g_pass)
    call fill(pass_left, g_array, g_pass)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, il, it, ik) &
    !$OMP SHARED(g_lo, trapped_wfb, g_array, g_pass, r_links, ntgrid) &
    !$OMP REDUCTION(max: err) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       il = il_idx(g_lo, iglo)
       ! Don't check for trapped particles
       if (il_is_trapped(il)) cycle
       if (il_is_wfb(il) .and. trapped_wfb) cycle

       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)

       if (r_links(it, ik) > 0) then
          err = max(err, abs(g_array(ntgrid, 2, iglo) - g_pass(ntgrid, 2, iglo)))
       end if
    end do
    !$OMP END PARALLEL DO

    ! Use a real to represent the outcome so we can use existing mp routines
    violation_flag = err
    call max_allreduce(violation_flag)
    if (violation_flag > local_tolerance) then
       call mp_abort('The passed array has inconsistent linked boundaries. '//id_message, .true.)
    end if
  end subroutine check_linked_boundaries_are_satisfied

  !> Calculate various moments of the distribution function
  subroutine getmoms (phinew, bparnew, ntot, density, upar, tpar, tperp, qparflux, pperpj1, qpperpj1)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_adjust, to_g_gs2, from_g_gs2, g_work
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, energy
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    !> Total perturbed density
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: ntot
    !> Non-adiabatic part of the perturbed density. Normalised to \(n_s n_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: density
    !> Parallel velocity. Normalised to \(v_{t,s} = v_{t,ref}\sqrt{T_s / m_s}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: upar
    !> Parallel temperature. Normalised to \(T_s T_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tpar
    !> Perpendicular temperature. Normalised to \(T_s T_{ref}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp
    !> Parallel heat flux. Normalised to \(n_s n_{ref} T_s T_{ref} v_{t,s}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: qparflux
    !> Pressure part of the particle flux. Normalised to \(\frac{n_s T_s}{q_s} \frac{n_{ref} T_{ref}}{q_{ref}}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: pperpj1
    !> Pressure part of the heat flux. Normalised to \(\frac{n_s T_s^2}{q_s} \frac{n_{ref} T_{ref}^2}{q_{ref}}\)
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: qpperpj1
    !> Electrostatic potential, \(\phi\)
    complex, dimension (-ntgrid:,:,:), intent(in) :: phinew
    !> Parallel magnetic field, \(B_\parallel\)
    complex, dimension (-ntgrid:,:,:), intent(in) :: bparnew

    integer :: ik, it, isgn, ie, is, iglo

! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
!                   from gnew, phinew and bparnew.
!           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
!           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
!
! use g_work as dist_fn dimensioned working space for all moments
! (avoid making more copies of gnew to save memory!)
!
! set gnew = g_wesson, but return gnew to entry state prior to exit 
    call g_adjust(gnew,phinew,bparnew, direction = from_g_gs2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate moments over the nonadiabatic part of <delta_f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g_work= J0 g_wesson = nonadiabatic piece of <delta_f>/F_m 
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

! CMR: density is the nonadiabatic piece of perturbed density
! NB normalised wrt equ'm density for species s: n_s n_ref  
!    ie multiply by (n_s n_ref) to get abs density pert'n
    call integrate_moment (g_work, density,  moment_to_allprocs, full_arr)

! DJA/CMR: upar and tpar moments 
! (nb adiabatic part of <delta f> does not contribute to upar, tpar or tperp)
! NB UPAR is normalised to vt_s = sqrt(T_s/m_s) vt_ref
!    ie multiply by spec(is)%stm * vt_ref to get abs upar
    g_work = vpa*g_work
    call integrate_moment (g_work, upar,  moment_to_allprocs, full_arr)

    g_work = 2.*vpa*g_work
    call integrate_moment (g_work, tpar,  moment_to_allprocs, full_arr)
! tpar transiently stores ppar, nonadiabatic perturbed par pressure 
!      vpa normalised to: sqrt(2 T_s T_ref/m_s m_ref)
!  including factor 2 in g_work product ensures 
!     ppar normalised to: n_s T_s n_ref T_ref 
!                         ppar = tpar + density, and so:
    tpar = tpar - density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = vperp2(:,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)
       end do
    end do
    call integrate_moment (g_work, tperp,  moment_to_allprocs, full_arr)
! tperp transiently stores pperp, nonadiabatic perturbed perp pressure
!                          pperp = tperp + density, and so:
    tperp = tperp - density
! NB TPAR, and TPERP are normalised by T_s T_ref
!    ie multiply by T_s T_ref to get abs TPAR, TPERP

! Now compute QPARFLUX
! NB QPARFLUX is normalised to n_s n_ref T_s T_ref v_ts
!    ie multiply by (n_s T_s spec(is)%stm) n_ref T_ref vt_ref to get abs qparflux
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = vpa(:,isgn,iglo)*gnew(:,isgn,iglo)*aj0(:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
    end do 
    call integrate_moment (g_work, qparflux,  moment_to_allprocs, full_arr)
   
! Now compute PPERPJ1, a modified p_perp which gives particle flux from Bpar
! NB PPERPJ1 is normalised to (n_s T_s/q_s)  n_ref T_ref/q_ref 
!    ie multiply by (n_s spec(is)%tz) n_ref T_ref/q_ref to get abs PPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) &
               = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
       end do
    end do
    call integrate_moment (g_work, pperpj1,  moment_to_allprocs, full_arr)

! Now compute QPPERPJ1, a modified p_perp*energy which gives heat flux from Bpar
! NB QPPERPJ1 is normalised to (n_s T_s^2/q_s)  n_ref  T_ref^2/q_ref
!    ie multiply by (n_s T_s spec(is)%tz) n_ref T_ref^2/q_ref to get abs QPPERPJ1
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
    end do
    call integrate_moment (g_work, qpperpj1, moment_to_allprocs, full_arr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now include the adiabatic part of <delta f>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set g_work = <delta_f>/F_m, including the adiabatic term
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - phinew(:,it,ik)*spec(is)%zt*nonmaxw_corr(ie,is)
       end do
    end do

! total perturbed density
    call integrate_moment (g_work, ntot,  moment_to_allprocs, full_arr)

!CMR, now multiply by species dependent normalisations to leave only reference normalisations
    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       density(:,:,:,is)=density(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
       qparflux(:,:,:,is)=qparflux(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%stm
       pperpj1(:,:,:,is)=pperpj1(:,:,:,is)*spec(is)%dens*spec(is)%tz
       qpperpj1(:,:,:,is)=qpperpj1(:,:,:,is)*spec(is)%dens*spec(is)%temp*spec(is)%tz
    end do

! return gnew to its initial state, the variable evolved in GS2
    call g_adjust(gnew,phinew,bparnew, direction = to_g_gs2)
  end subroutine getmoms

  !> FIXME : Add documentation
  subroutine getmoms_gryfx_dist (density_gryfx, upar_gryfx, tpar_gryfx, tperp_gryfx, qpar_gryfx, qperp_gryfx, phi_gryfx)
    use dist_fn_arrays, only: vpa, vperp2, gnew, g_work
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment
    use fields_arrays, only: phinew
    use kt_grids, only: ntheta0, naky
    use mp, only: proc0
    use iso_fortran_env, only: real32
    implicit none
    integer :: ik, it, isgn, is, iglo, ig, iz, index_gryfx
    complex(real32), dimension(:), intent(out) :: density_gryfx, upar_gryfx, &
         tpar_gryfx, tperp_gryfx, qpar_gryfx, qperp_gryfx
    complex(real32), dimension(:), intent(out) :: phi_gryfx
    complex, dimension(:,:,:,:), allocatable :: total

    allocate(total(-ntgrid:ntgrid,ntheta0,naky,nspec))
   
    ! gryfx uses vt=sqrt(T/m) normalization, 
    ! so vt_gs2 = sqrt(2T/m) = sqrt(2) vt_gryfx.
    ! since vpa and vperp2 are normalized to vt_gs2, below we multiply by
    ! factors of sqrt(2) to get vpa and vperp2 in gryfx units before taking
    ! integrals, i.e.
    ! vpa_gryfx = sqrt(2) vpa
    ! vperp2_gryfx = 2 vperp2

    ! dens = < f >
    g_work = gnew 
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            density_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! upar = < (vpar/vt_gryfx) f > = < ( sqrt(2) vpar/vt_gs2 ) f > = upar1/vt_gryfx
    g_work = sqrt(2.)*vpa*gnew  
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            upar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! tpar = < (vpar^2/vt_gryfx^2 - 1) f > = tpar1/T_0
    g_work = (2.*vpa*vpa - 1.)*gnew
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            tpar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! qpar = < (vpar^3/vt_gryfx^3 - 3*vpar/vt_gryfx) f > = qpar1/( n0 m vt_gryfx^3 )
    g_work = (2.*sqrt(2.)*vpa*vpa*vpa - 3.*sqrt(2.)*vpa)*gnew
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            qpar_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! tprp = < ( 1/2 vprp/vt_gryfx^2 - 1 ) f > = tprp1/T_0
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = (vperp2(:,iglo)-1.)*gnew(:,isgn,iglo) 
       end do
    end do
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            tperp_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! qprp = < vpar/vt_gryfx*( 1/2 vprp/vt_gryfx^2 - 1 ) f > = qprp1/( n0 m vt_gryfx^3 )
    g_work = sqrt(2.)*vpa*g_work
    call integrate_moment (g_work, total)
    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1,g_lo%ntheta0
        do ik = 1, g_lo%naky
          do is = 1,g_lo%nspec
            index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + &
                            g_lo%naky*g_lo%ntheta0*(iz-1) + &
                            2*ntgrid*g_lo%naky*g_lo%ntheta0*(is-1)
            qperp_gryfx(index_gryfx) = total(ig, it, ik, is)
          end do
        end do
      end do
    end do
    end if

    ! we also must account for the overall a/rho_i factor on all quantities. 
    ! so we multiply all moments and phi by an additional factor of sqrt(2)
    if(proc0 ) then 
      ! normalize to gryfx units
      density_gryfx = sqrt(2.)*density_gryfx
      upar_gryfx = sqrt(2.)*upar_gryfx
      tpar_gryfx = sqrt(2.)*tpar_gryfx
      tperp_gryfx = sqrt(2.)*tperp_gryfx
      qpar_gryfx = sqrt(2.)*qpar_gryfx
      qperp_gryfx = sqrt(2.)*qperp_gryfx
    endif
      

    if(proc0) then
    do ig = -ntgrid, ntgrid-1
      iz = ig + ntgrid + 1
      do it = 1, g_lo%ntheta0
        do ik = 1, g_lo%naky
          index_gryfx = 1 + (ik-1) + g_lo%naky*((it-1)) + g_lo%naky*g_lo%ntheta0*(iz-1)
          phi_gryfx(index_gryfx) = phinew(ig, it, ik)*sqrt(2.)
        end do
      end do
    end do
    end if

    deallocate(total)


  end subroutine getmoms_gryfx_dist

  !> Calculates hermite-laguerre moments of g
  !! uses gryfx normalization convention (see getmoms_gryfx above)
  subroutine get_mom_glm (glm, l, m)
    use dist_fn_arrays, only: vpa, vperp2, gnew, g_work
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment, hermite_prob
    use gauss_quad, only: laguerre_norm

    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: glm
    integer, intent(in) :: l, m
    integer :: isgn, iglo, ig

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
      do isgn = 1, 2
        do ig = -ntgrid, ntgrid
          g_work(ig,isgn,iglo) = hermite_prob(l, sqrt(2.0)*vpa(ig,isgn,iglo)) * &
                        laguerre_norm(m, vperp2(ig,iglo)) * gnew(ig,isgn,iglo) 
        end do
      end do
    end do
    call integrate_moment (g_work, glm,  moment_to_allprocs, full_arr)
  end subroutine get_mom_glm

  !> Returns electron density and Tperp moment integrals to PE 0
  subroutine getemoms (phinew, bparnew, ntot, tperp)
    use dist_fn_arrays, only: vperp2, aj0, gnew, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_layouts, only: is_idx, ie_idx, g_lo, ik_idx, it_idx
    use species, only: nspec, spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_moment

    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, dimension (-ntgrid:,:,:,:), intent (out) :: tperp, ntot
    complex, dimension (-ntgrid:,:,:), intent(in) :: phinew, bparnew

    integer :: ik, it, isgn, ie, is, iglo

!
! What are <delta_f> and g_wesson in the note below?
! g_wesson = <delta_f> + q phi/T    [ignore F_m factors for simplicity]
!
! Electrostatically (for simplicity), g_adjust produces:
!
! h = g_gs2 + q <phi> / T
! 
! then in the subsequent code they calculate for ntot:
!
! ntot = integral[   J0 h - q phi / T  ]
!
! so g_wesson == h.  What is odd in our notation is the LHS.  
! We typically indicate the perturbed distribution at fixed spatial position 
! by delta_f.  In the note below, they must mean delta_f = delta_f (R), so that 
! they are gyro-averaging to get the distribution at fixed spatial position.
!
! In summary, DJA and CMR are calculating the moments at fixed spatial position
! rather than at fixed guiding centers, and using different notation than appears
! in most of our papers.
!
! DJA+CMR: 17/1/06, use g_adjust routine to extract g_wesson
!                   from gnew, phinew and bparnew.
!           nb  <delta_f> = g_wesson J0 - q phi/T F_m  where <> = gyroaverage
!           ie  <delta_f>/F_m = g_wesson J0 - q phi/T
!
! use g_work as dist_fn dimensioned working space for all moments
! (avoid making more copies of gnew to save memory!)
!
! set gnew = g_wesson, but return gnew to entry state prior to exit 
    call g_adjust(gnew, phinew, bparnew, direction = from_g_gs2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo) - phinew(:,it,ik)*spec(is)%zt*nonmaxw_corr(ie,is)
       end do
    end do

! total perturbed density
    call integrate_moment (g_work, ntot,  moment_to_allprocs, full_arr)

! vperp**2 moment:
    do iglo = g_lo%llim_proc, g_lo%ulim_proc       
       ie = ie_idx(g_lo,iglo) ; is = is_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo) ; it = it_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)*vperp2(:,iglo) &
             - phinew(:,it,ik)*spec(is)%zt*nonmaxw_corr(ie,is)*vperp2(:,iglo)
       end do
    end do

! total perturbed perp pressure
    call integrate_moment (g_work, tperp,  moment_to_allprocs, full_arr)

! tperp transiently stores pperp, 
!       pperp = tperp + density, and so:
    tperp = tperp - ntot

    do is=1,nspec
       ntot(:,:,:,is)=ntot(:,:,:,is)*spec(is)%dens
       tperp(:,:,:,is)=tperp(:,:,:,is)*spec(is)%temp
    end do

! return gnew to its initial state, the variable evolved in GS2
    call g_adjust(gnew,phinew,bparnew, direction = to_g_gs2)

  end subroutine getemoms
  
  !> Calculates moments at not guiding center coordinate
  subroutine getmoms_notgc (phinew, bparnew, dens, upar, tpar, tper, ntot, jpar, reset_coefficients)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, gnew, g_work
    use gs2_layouts, only: g_lo, is_idx, ik_idx, it_idx, ie_idx
    use species, only: nspec, spec, nonmaxw_corr
    use theta_grid, only: ntgrid
    use kt_grids, only: nakx => ntheta0, naky
    use le_grids, only: integrate_moment
    use optionals, only: get_option_with_default
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    complex, intent (out) :: &
         & dens(-ntgrid:,:,:,:), upar(-ntgrid:,:,:,:), &
         & tpar(-ntgrid:,:,:,:), tper(-ntgrid:,:,:,:)
    complex, intent (out), optional :: ntot(-ntgrid:,:,:,:)
    complex, intent (out), optional :: jpar(-ntgrid:,:,:)
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew, bparnew
    logical, intent(in), optional :: reset_coefficients
    real, dimension(:, :, :, :), allocatable, save :: mom_coeff
    real, dimension(:, :, :), allocatable, save :: mom_coeff_npara, mom_coeff_nperp
    real, dimension(:, :, :), allocatable, save :: mom_coeff_tpara, mom_coeff_tperp
    real, dimension(:, :, :), allocatable, save :: mom_shift_para, mom_shift_perp

    integer :: isgn, iglo, is, ie

    real :: a, b, tpar2, tper2
    integer :: it, ik, ig

    ! returns moment integrals to PE 0

    if (get_option_with_default(reset_coefficients, .false.)) then
       if (allocated(mom_coeff)) deallocate(mom_coeff)
       if (allocated(mom_coeff_npara)) deallocate(mom_coeff_npara)
       if (allocated(mom_coeff_nperp)) deallocate(mom_coeff_nperp)
       if (allocated(mom_coeff_tpara)) deallocate(mom_coeff_tpara)
       if (allocated(mom_coeff_tperp)) deallocate(mom_coeff_tperp)
       if (allocated(mom_shift_para)) deallocate(mom_shift_para)
       if (allocated(mom_shift_perp)) deallocate(mom_shift_perp)
    end if

    ! Note we mark mom_coeff* as 'save' so that the result of this call is cached
    ! and not recalculated each time. This can help save overhead and should be
    ! fine. The result of this call should only depend on the velocity grids, the
    ! wavenumbers and the species parameters (and the derived Bessel functions).
    ! These should not change during a single standard run, but if we wish to enable
    ! more exotic behaviour we might want to remove the save from above declarations and
    ! recalculate these coefficients each time. For now just add an optional flag to
    ! request recalculation.
    call init_mom_coeff(mom_coeff, mom_coeff_npara, mom_coeff_nperp, mom_coeff_tpara, &
       mom_coeff_tperp, mom_shift_para, mom_shift_perp)

    ! not guiding center n_total
    if(present(ntot)) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          it = it_idx(g_lo,iglo)
          ie = ie_idx(g_lo,iglo)

          do isgn = 1, 2
             g_work(:,isgn,iglo) = aj0(:,iglo) * gnew(:,isgn,iglo)
          end do
          do isgn = 1, 2
             g_work(:,isgn,iglo) = g_work(:,isgn,iglo) + phinew(:,it,ik) &
                  & *(aj0(:,iglo)**2-1.0) * spec(is)%zt * nonmaxw_corr(ie,is)
          end do
          do isgn = 1, 2
             g_work(:,isgn,iglo) = g_work(:,isgn,iglo) &
                  & + 2.*vperp2(:,iglo)*aj1(:,iglo)*aj0(:,iglo) &
                  & * bparnew(:,it,ik) * nonmaxw_corr(ie,is)
          end do
       end do
       call integrate_moment (g_work, ntot,  moment_to_allprocs,full_arr)
    endif

! not guiding center density
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, dens,  moment_to_allprocs,full_arr)

! not guiding center upar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, upar,  moment_to_allprocs,full_arr)

! not guiding center tpar
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = 2.*aj0(:,iglo)*vpa(:,isgn,iglo)**2*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, tpar,  moment_to_allprocs,full_arr)
    
! not guiding center tperp
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = 2.*vperp2(:,iglo)*aj1(:,iglo)*gnew(:,isgn,iglo)
       end do
    end do

    call integrate_moment (g_work, tper,  moment_to_allprocs,full_arr)

    do ig=-ntgrid,ntgrid
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                tpar2=tpar(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_npara(it,ik,is)
                tper2=tper(ig,it,ik,is) &
                     & - dens(ig,it,ik,is)*mom_coeff_nperp(it,ik,is)

                a=mom_coeff_tperp(it,ik,is)
                b=mom_coeff_tpara(it,ik,is)

                tpar(ig,it,ik,is)=(   tpar2-a*tper2)/(1.-a*b)
                tper(ig,it,ik,is)=(-b*tpar2+  tper2)/(1.-a*b)
             end do
          end do
       end do
    end do

    do is=1,nspec
       dens(:,:,:,is)=dens(:,:,:,is)*spec(is)%dens
       upar(:,:,:,is)=upar(:,:,:,is)*spec(is)%stm
       tpar(:,:,:,is)=tpar(:,:,:,is)*spec(is)%temp
       tper(:,:,:,is)=tper(:,:,:,is)*spec(is)%temp
    end do

    if(present(jpar)) then
       jpar(:,:,:)=cmplx(0.,0.)
       do is=1,nspec
          jpar(:,:,:)=jpar(:,:,:)+spec(is)%z*spec(is)%dens*upar(:,:,:,is)
       end do
    endif
  end subroutine getmoms_notgc

  !> FIXME : Add documentation  
  subroutine init_fieldeq
    use dist_fn_arrays, only: aj0, aj1, vperp2, g_work
    use species, only: nspec, spec, has_electron_species, has_ion_species, nonmaxw_corr
    use theta_grid, only: ntgrid
    use kt_grids, only: naky, ntheta0, aky, kperp2
    use le_grids, only: integrate_species
    use gs2_layouts, only: g_lo, ie_idx, is_idx
    use run_parameters, only: tite
    implicit none
    integer :: iglo, isgn
    integer :: ik, it, ie, is
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: tot
    real, dimension (nspec) :: wgt

    if (feqinit) return
    feqinit = .true.

    allocate (gridfac1(-ntgrid:ntgrid,ntheta0,naky))
    gridfac1 = 1.0
    select case (boundary_option_switch)
    case (boundary_option_self_periodic)
       ! nothing
    case (boundary_option_linked)
       do it = 1, ntheta0
          do ik = 1, naky
             if (aky(ik) == 0.0) cycle
             if (itleft(it, ik) < 0) gridfac1(-ntgrid,it,ik) = gridfac
             if (itright(it, ik) < 0) gridfac1(ntgrid,it,ik) = gridfac
          end do
       end do
    case default
       do ik = 1, naky
          if (aky(ik) == 0.0) cycle
          gridfac1(-ntgrid,:,ik) = gridfac
          gridfac1(ntgrid,:,ik) = gridfac
       end do
    end select

    allocate (gamtot(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot1(-ntgrid:ntgrid,ntheta0,naky))
    allocate (gamtot2(-ntgrid:ntgrid,ntheta0,naky))
    if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
       allocate (gamtot3(-ntgrid:ntgrid,ntheta0,naky))
    endif
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = (1.0 - aj0(:,iglo)**2)*nonmaxw_corr(ie,is)
       end do
    end do
    wgt = spec%z*spec%z*spec%dens/spec%temp
    call integrate_species (g_work, wgt, tot)
    gamtot = real(tot) + kperp2*poisfac
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*aj1(:,iglo) &
               *2.0*vperp2(:,iglo) * nonmaxw_corr(ie,is)
       end do
    end do
    wgt = spec%z*spec%dens
    call integrate_species (g_work, wgt, tot)
    gamtot1 = real(tot)
    
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       ie = ie_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj1(:,iglo)**2*2.0*vperp2(:,iglo)**2*nonmaxw_corr(ie,is)
       end do
    end do
    wgt = spec%temp*spec%dens
    call integrate_species (g_work, wgt, tot)
    gamtot2 = real(tot)
    
    !<DD>Make sure gamtots are single valued
    if(esv)then
       call ensure_single_val_fields_pass_r(gamtot)
       call ensure_single_val_fields_pass_r(gamtot1)
       call ensure_single_val_fields_pass_r(gamtot2)
    endif

! adiabatic electrons 
    if (.not. has_electron_species(spec) .or. .not. has_ion_species(spec) ) then
       if (adiabatic_option_switch == adiabatic_option_yavg) then
          do ik = 1, naky
             if (aky(ik) > epsilon(0.0)) gamtot(:,:,ik) = gamtot(:,:,ik) + tite
          end do
       elseif (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          gamtot  = gamtot + tite
          gamtot3 = (gamtot-tite) / gamtot
          where (gamtot3 < 2.*epsilon(0.0)) gamtot3 = 1.0
       else
          gamtot = gamtot + tite 
       endif
    endif

    
  end subroutine init_fieldeq

  !> FIXME : Add documentation    
  subroutine getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp, local_only)
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0, kperp2
    use run_parameters, only: has_phi, has_apar, has_bpar, beta
    use species, only: spec, has_electron_species
    use optionals, only: get_option_with_default
    use gs2_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (in out) ::fieldeq,fieldeqa,fieldeqp
    logical, intent(in), optional :: local_only
    logical :: local_local_only
    integer :: ik, it
    integer :: it_llim, it_ulim, ik_llim, ik_ulim

    local_local_only = get_option_with_default(local_only, .false.)

    if (local_local_only) then
       it_llim = g_lo%it_min ; it_ulim = g_lo%it_max
       ik_llim = g_lo%ik_min ; ik_ulim = g_lo%ik_max
    else
       it_llim = 1 ; it_ulim = ntheta0
       ik_llim = 1 ; ik_ulim = naky
    end if

    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          call calculate_flux_surface_average(fl_avg,antot)
       end if
    end if

    if (has_phi) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ik, it) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, fieldeq, antot, &
       !$OMP bpar, gamtot1, gamtot, gridfac1, phi, fl_avg) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             fieldeq(:,it,ik) = antot(:,it,ik) + &
                  bpar(:,it,ik)*gamtot1(:,it,ik) - &
                  gamtot(:,it,ik)*gridfac1(:,it,ik)*phi(:,it,ik) + &
                  fl_avg(it,ik)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (has_apar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ik, it) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, fieldeqa, antota, &
       !$OMP kperp2, gridfac1, apar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             fieldeqa(:, it, ik) = antota(:, it, ik) - &
                  kperp2(:, it, ik)*gridfac1(:, it, ik)*apar(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    ! bpar == delta B_parallel / B_0(theta) b/c of the factor of 1/bmag(theta)**2
    ! in the following
    if (has_bpar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ik, it) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, fieldeqp, antotp, bpar, &
       !$OMP gamtot2, phi, gamtot1, beta, apfac, bmag, gridfac1)  &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             fieldeqp(:,it,ik) = (antotp(:, it, ik) + &
                  bpar(:, it, ik)*gamtot2(:, it, ik) + &
                  0.5*phi(:, it, ik)*gamtot1(:, it, ik))*(beta*apfac/bmag**2) + &
                  bpar(:, it, ik)*gridfac1(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine getfieldeq1

  !> FIXME : Add documentation    
  subroutine getfieldeq (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp, local_only)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: antot, antota, antotp
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) ::fieldeq,fieldeqa,fieldeqp
    logical, intent(in), optional :: local_only

    call getan (antot, antota, antotp, local_only)
    call getfieldeq1 (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp, local_only)
  end subroutine getfieldeq


  !///////////////////////////////////////
  !// SOME NO GATHER TEST ROUTINES
  !///////////////////////////////////////
  
  !> FIXME : Add documentation    
  subroutine getfieldeq_nogath (phi, apar, bpar, fieldeq, fieldeqa, fieldeqp, local_only)
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: antot, antota, antotp
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in out) ::fieldeq,fieldeqa,fieldeqp
    logical, intent(in), optional :: local_only

    call getan_nogath (antot, antota, antotp, local_only)

    call getfieldeq1_nogath (phi, apar, bpar, antot, antota, antotp, &
         fieldeq, fieldeqa, fieldeqp, local_only)
  end subroutine getfieldeq_nogath

  !> FIXME : Add documentation
  !>
  !> This routine is near identical to getfieldeq1
  subroutine getfieldeq1_nogath (phi, apar, bpar, antot, antota, antotp, &
       fieldeq, fieldeqa, fieldeqp, local_only)
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0, kperp2
    use run_parameters, only: has_phi, has_apar, has_bpar, beta
    use kt_grids, only: kwork_filter
    use species, only: spec, has_electron_species
    use optionals, only: get_option_with_default
    use gs2_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: antot, antota, antotp
    complex, dimension (-ntgrid:,:,:), intent (in out) ::fieldeq,fieldeqa,fieldeqp
    logical, intent(in), optional :: local_only
    logical :: local_local_only
    integer :: ik, it
    integer :: it_llim, it_ulim, ik_llim, ik_ulim

    local_local_only = get_option_with_default(local_only, .false.)

    if (local_local_only) then
       it_llim = g_lo%it_min ; it_ulim = g_lo%it_max
       ik_llim = g_lo%ik_min ; ik_ulim = g_lo%ik_max
    else
       it_llim = 1 ; it_ulim = ntheta0
       ik_llim = 1 ; ik_ulim = naky
    end if

    if (.not. allocated(fl_avg)) allocate (fl_avg(ntheta0, naky))
    fl_avg = 0.

    if (.not. has_electron_species(spec)) then
       if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
          call calculate_flux_surface_average(fl_avg,antot)
       end if
    end if

    if (has_phi) then
       if (has_bpar) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ik, it) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, kwork_filter, fieldeq, &
          !$OMP antot, bpar, gamtot1, gamtot, gridfac1, phi, fl_avg) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                if(kwork_filter(it,ik)) cycle
                fieldeq(:,it,ik) = antot(:,it,ik) + &
                     bpar(:,it,ik)*gamtot1(:,it,ik) - &
                     gamtot(:,it,ik)*gridfac1(:,it,ik)*phi(:,it,ik) + &
                     fl_avg(it,ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ik, it) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, kwork_filter, fieldeq, &
          !$OMP antot, gamtot, gridfac1, phi, fl_avg) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                if(kwork_filter(it,ik)) cycle
                fieldeq(:,it,ik) = antot(:,it,ik) - &
                     gamtot(:,it,ik)*gridfac1(:,it,ik)*phi(:,it,ik) + &
                     fl_avg(it,ik)
             end do
          end do
          !$OMP END PARALLEL DO
       end if
    end if

    if (has_apar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(ik, it) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, kwork_filter, fieldeqa, &
       !$OMP antota, kperp2, gridfac1, apar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             if(kwork_filter(it,ik)) cycle
             fieldeqa(:, it, ik) = antota(:, it, ik) - &
                  kperp2(:, it, ik)*gridfac1(:, it, ik)*apar(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    ! bpar == delta B_parallel / B_0(theta) b/c of the factor of 1/bmag(theta)**2
    ! in the following
    if (has_bpar) then
       if (has_phi) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ik, it) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, kwork_filter, gridfac1, &
          !$OMP fieldeqp, antotp, bpar, gamtot2, phi, gamtot1, beta, apfac, bmag) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                if(kwork_filter(it,ik)) cycle
                fieldeqp(:,it,ik) = (antotp(:, it, ik) + &
                     bpar(:, it, ik)*gamtot2(:, it, ik) + &
                     0.5*phi(:, it, ik)*gamtot1(:, it, ik))*(beta*apfac/bmag**2) + &
                     bpar(:, it, ik)*gridfac1(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(ik, it) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, kwork_filter, &
          !$OMP fieldeqp, antotp, bpar, gamtot2, beta, apfac, bmag, gridfac1) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                if(kwork_filter(it,ik)) cycle
                fieldeqp(:,it,ik) = (antotp(:, it, ik) + &
                     bpar(:, it, ik)*gamtot2(:, it, ik))*(beta*apfac/bmag**2) + &
                     bpar(:, it, ik)*gridfac1(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       end if
    end if
  end subroutine getfieldeq1_nogath

  !> FIXME : Add documentation    
  subroutine check_getan (antot, antota, antotp, tempantot, tempantota, tempantotp)
    use mp, only : iproc
    use theta_grid, only: ntgrid
    use run_parameters, only: has_phi, has_apar, has_bpar
    use kt_grids, only: naky, ntheta0, kwork_filter
    use gs2_layouts, only: g_lo, it_idx, ik_idx, gf_lo, idx, proc_id

    implicit none

    complex, dimension (-ntgrid:ntgrid,ntheta0,naky), intent (inout) :: antot, antota, antotp
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky), intent (inout) :: tempantot, tempantota, tempantotp
    complex :: tol = (1e-14,1e-14)
    integer :: it, ik, igf, ig, iglo

    if (has_phi) then

       if(any(kwork_filter))then
          loop1: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if(kwork_filter(it,ik)) cycle
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antot(ig,it,ik)-tempantot(ig,it,ik))).gt.aimag(tol) .or. abs(real(antot(ig,it,ik)-tempantot(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antot gf_lo integration',igf,it,ik,antot(ig,it,ik),tempantot(ig,it,ik)
                      exit loop1
                   end if
                end do
             end if
          end do loop1
       else
          loop2: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antot(ig,it,ik)-antot(ig,it,ik))).gt.aimag(tol) .or. abs(real(antot(ig,it,ik)-antot(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antot gf_lo integration',igf,it,ik,antot(ig,it,ik),tempantot(ig,it,ik)
                      exit loop2
                   end if
                end do
             end if
          end do loop2
       end if

    end if

    if (has_apar) then

       if(any(kwork_filter))then
          loop3: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if(kwork_filter(it,ik)) cycle
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antota(ig,it,ik)-tempantota(ig,it,ik))).gt.aimag(tol) .or. abs(real(antota(ig,it,ik)-tempantota(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antota gf_lo integration',igf,it,ik,antota(ig,it,ik),tempantota(ig,it,ik)
                      exit loop3
                   end if
                end do
             end if
          end do loop3
       else
          loop4: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antota(ig,it,ik)-antota(ig,it,ik))).gt.aimag(tol) .or. abs(real(antota(ig,it,ik)-antota(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antota gf_lo integration',igf,it,ik,antota(ig,it,ik),tempantota(ig,it,ik)
                      exit loop4
                   end if
                end do
             end if
          end do loop4
       end if


    end if


    if (has_bpar) then

       if(any(kwork_filter))then
          loop5: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             if(kwork_filter(it,ik)) cycle
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antotp(ig,it,ik)-tempantotp(ig,it,ik))).gt.aimag(tol) .or. abs(real(antotp(ig,it,ik)-tempantotp(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antotp gf_lo integration',igf,it,ik,antotp(ig,it,ik),tempantotp(ig,it,ik)
                      exit loop5
                   end if
                end do
             end if
          end do loop5
       else
          loop6: do iglo = g_lo%llim_proc, g_lo%ulim_proc
             ik = ik_idx(g_lo,iglo)
             it = it_idx(g_lo,iglo)
             igf = idx(gf_lo,ik,it)
             if(proc_id(gf_lo,igf) .eq. iproc) then
                do ig=-ntgrid,ntgrid
                   if(abs(aimag(antotp(ig,it,ik)-antotp(ig,it,ik))).gt.aimag(tol) .or. abs(real(antotp(ig,it,ik)-antotp(ig,it,ik))).gt.real(tol)) then
                      write(*,*) iproc,'problem with antotp gf_lo integration',igf,it,ik,antotp(ig,it,ik),tempantotp(ig,it,ik)
                      exit loop6
                   end if
                end do
             end if
          end do loop6
       end if

    end if

  end subroutine check_getan

  !> Getan_nogath has been substantially changed to enable gf_lo field calculations, it can now operate in
  !> 2 modes; using the standard integrate_species (integrate_species_sub in le_grids) when called from
  !> fields_local and gf_lo_integrate is false (which is teh default, or doing the integrate in place when
  !> called from fields_gf_local or from fields_local when gf_lo_integrate is true.  Note, if this is called from
  !> fields_local and gf_lo_integrate is true then the calculation will be done locally but there is a function
  !> called after this has finished in fields_local that sends the data back (from gf_lo to g_lo layouts).
  !> When call with gf_lo_integrate = .true. this routine does a gather that converts the gnew array from g_lo
  !> data distribution to gf_lo data distribution and stores the result in gfarray.  With gf_lo_integrate = .false.
  !> g_in is used instead and this is in g_lo data distribution.
  !>
  !> AJ
  subroutine getan_nogath_from_dfn (g_in, antot, antota, antotp, local_only)
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1, g_work
    use dist_fn_arrays, only: vpa_gf, vperp2_gf, aj0_gf, aj1_gf
    use species, only: nspec, spec
    use theta_grid, only: ntgrid
    use le_grids, only: integrate_species, g2gf, negrid, nlambda, w, wl
    use run_parameters, only: beta, has_phi, has_apar, has_bpar
    use gs2_layouts, only: g_lo, it_idx,ik_idx, gf_lo, proc_id, idx, ie_idx, is_idx, il_idx
    use kt_grids, only: kwork_filter, kperp2, naky, ntheta0
    use redistribute, only: gather
    use mp, only: barrier
    use optionals, only: get_option_with_default
    implicit none
    complex, dimension (-ntgrid:, :, g_lo%llim_proc:), intent(in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in out) :: antot, antota, antotp
    logical, optional :: local_only
    real, dimension (nspec) :: wgt
    integer :: isgn, iglo, it, ik, igf
    complex, dimension (:, :, :, :, :, :), allocatable :: gfarray
    integer :: ie, il, is
    real, dimension(-ntgrid:ntgrid) :: temparg3
    logical :: local_local_only
    integer :: it_llim, it_ulim, ik_llim, ik_ulim

    local_local_only = get_option_with_default(local_only, .false.)
    if (local_local_only) then
       it_llim = g_lo%it_min ; it_ulim = g_lo%it_max
       ik_llim = g_lo%ik_min ; ik_ulim = g_lo%ik_max
    else
       it_llim = 1 ; it_ulim = ntheta0
       ik_llim = 1 ; ik_ulim = naky
    end if

    ! Initialise the outputs if with gf_lo_integrate. We note that
    ! the arrays are initialised on allocation so we do not need
    ! to zero in the gf_lo_integrate = .false. case even though
    ! integrate species may not populate all indices.
    if (gf_lo_integrate) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, antot) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             antot(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, antota) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             antota(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, antotp) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             antotp(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO

       !AJ if we are performing the integrate in gf_lo space (i.e. the result will be returned in gf_lo)
       !AJ then do the initial redistribution of the data from g_lo (where g_in is current) to gf_lo
       !AJ where it will end up in gfarray)
       allocate( gfarray(-ntgrid:ntgrid,2,nspec,negrid,nlambda,gf_lo%llim_proc:gf_lo%ulim_alloc))
       call gather(g2gf, g_in, gfarray, ntgrid)
    end if

    if (has_phi) then
       !<DD>NOTE: It's possible to rewrite this loop as simply
       !g_work=g_in*spread(aj0,2,2)
       !but this seems to be slower than an explicit loop and
       !the ability to skip certain it/ik values is lost.

       wgt = spec%z*spec%dens

       if(.not. gf_lo_integrate) then
          if(any(kwork_filter))then
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn, it, ik) &
             !$OMP SHARED(g_lo, kwork_filter, g_work, aj0, g_in) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                it=it_idx(g_lo,iglo)
                ik=ik_idx(g_lo,iglo)
                if(kwork_filter(it,ik))cycle
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj0(:,iglo)*g_in(:,isgn,iglo)
                end do
             end do
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn) &
             !$OMP SHARED(g_lo, g_work, aj0, g_in) &
             !$OMP COLLAPSE(2) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj0(:,iglo)*g_in(:,isgn,iglo)
                end do
             end do
            !$OMP END PARALLEL DO
          endif
          call integrate_species (g_work, wgt, antot, nogath=.true.)
          if(esv) call ensure_single_val_fields_pass(antot)
       else
          !AJ The code below undertakes the functionality of integrate_species (or specifically integrate_species_gf_nogather from le_grids).
          !AJ It is possible to call that routine instead of the code below, but the code was brought into here to enable the combinations of the
          !AJ pre-calculation loop (see the call above where aj0 is multipled by g_in) with the integrate species itself.  This should be more
          !AJ efficient as we then are only traversing the gf_array once for the integration of this field.  However, the opposite could also have been
          !AJ done (i.e. aj0_gf taken into integrate_species) and if that would be cleaner then that refactoring can be done.
          if(any(kwork_filter))then
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antot(:,it,ik) =0.
                if(kwork_filter(it,ik)) cycle
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         antot(:,it,ik) = antot(:,it,ik) + wgt(is)*w(ie,is)*wl(:,il)*((aj0_gf(:,is,ie,il,igf)*gfarray(:,1,is,ie,il,igf))+(aj0_gf(:,is,ie,il,igf)*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          else
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antot(:,it,ik) =0.
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         antot(:,it,ik) = antot(:,it,ik) + wgt(is)*w(ie,is)*wl(:,il)*((aj0_gf(:,is,ie,il,igf)*gfarray(:,1,is,ie,il,igf))+(aj0_gf(:,is,ie,il,igf)*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          endif
       end if

       if (afilter > epsilon(0.0)) antot = antot * exp(-afilter**4*kperp2**2/4.)
       !NOTE: We don't do ensure_single_val_fields here as we're not certain we
       !have the full data
    end if

    if (has_apar) then

       wgt = 2.0*beta*spec%z*spec%dens*spec%stm

       if(.not. gf_lo_integrate) then
          if(any(kwork_filter))then
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn, it, ik) &
             !$OMP SHARED(g_lo, kwork_filter, g_work, aj0, vpa, g_in) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                it=it_idx(g_lo,iglo)
                ik=ik_idx(g_lo,iglo)
                if(kwork_filter(it,ik))cycle
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*g_in(:,isgn,iglo)
                end do
             end do
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn) &
             !$OMP SHARED(g_lo, g_work, aj0, vpa, g_in) &
             !$OMP COLLAPSE(2) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj0(:,iglo)*vpa(:,isgn,iglo)*g_in(:,isgn,iglo)
                end do
             end do
             !$OMP END PARALLEL DO
          end if
          call integrate_species (g_work, wgt, antota, nogath=.true.)
          if(esv) call ensure_single_val_fields_pass(antota)
       else
          !AJ The code below undertakes the functionality of integrate_species (or specifically integrate_species_gf_nogather from le_grids).
          !AJ It is possible to call that routine instead of the code below, but the code was brought into here to enable the combinations of the
          !AJ pre-calculation loop (see the call above where aj0 and vpa are multipled by g_in) with the integrate species itself.  This should be more
          !AJ efficient as we then are only traversing the gf_array once for the integration of this field.  However, the opposite could also have been
          !AJ done (i.e. aj0_gf and vpa_gf taken into integrate_species) and if that would be cleaner then that refactoring can be done.
          if(any(kwork_filter))then
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antota(:,it,ik) =0.
                if(kwork_filter(it,ik)) cycle
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         antota(:,it,ik) = antota(:,it,ik) + &
                              wgt(is)*w(ie,is)*wl(:,il)*((aj0_gf(:,is,ie,il,igf)*vpa_gf(:,1,ie,il,is)*gfarray(:,1,is,ie,il,igf))+(aj0_gf(:,is,ie,il,igf)*vpa_gf(:,2,ie,il,is)*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          else
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antota(:,it,ik) =0.
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         antota(:,it,ik) = antota(:,it,ik) + &
                              wgt(is)*w(ie,is)*wl(:,il)*((aj0_gf(:,is,ie,il,igf)*vpa_gf(:,1,ie,il,is)*gfarray(:,1,is,ie,il,igf))+(aj0_gf(:,is,ie,il,igf)*vpa_gf(:,2,ie,il,is)*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          end if
       end if



       !NOTE: We don't do ensure_single_val_fields here as we're not certain we
       !have the full data
    end if

    if (has_bpar) then

       wgt = spec%temp*spec%dens

       if(.not. gf_lo_integrate) then

          if(any(kwork_filter))then
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn, it, ik) &
             !$OMP SHARED(g_lo, kwork_filter, g_work, aj1, vperp2, g_in) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                it=it_idx(g_lo,iglo)
                ik=ik_idx(g_lo,iglo)
                if(kwork_filter(it,ik))cycle
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj1(:,iglo)*vperp2(:,iglo)*g_in(:,isgn,iglo)
                end do
             end do
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO DEFAULT(none) &
             !$OMP PRIVATE(iglo, isgn) &
             !$OMP SHARED(g_lo, g_work, aj1, vperp2, g_in) &
             !$OMP COLLAPSE(2) &
             !$OMP SCHEDULE(static)
             do iglo = g_lo%llim_proc, g_lo%ulim_proc
                do isgn = 1, 2
                   g_work(:,isgn,iglo) = aj1(:,iglo)*vperp2(:,iglo)*g_in(:,isgn,iglo)
                end do
             end do
             !$OMP END PARALLEL DO
          end if
          call integrate_species (g_work, wgt, antotp, nogath=.true.)
          if(esv) call ensure_single_val_fields_pass(antotp)
       else
          !AJ The code below undertakes the functionality of integrate_species (or specifically integrate_species_gf_nogather from le_grids).
          !AJ It is possible to call that routine instead of the code below, but the code was brought into here to enable the combinations of the
          !AJ pre-calculation loop (see the call above where aj1 and vperp2 are multipled by g_in) with the integrate species itself.  This should be more
          !AJ efficient as we then are only traversing the gf_array once for the integration of this field.  However, the opposite could also have been
          !AJ done (i.e. aj1_gf and vperp2_gf taken into integrate_species) and if that would be cleaner then that refactoring can be done.
          if(any(kwork_filter))then
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antotp(:,it,ik) =0.
                if(kwork_filter(it,ik)) cycle
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         !AJ temparg3 is simply an attempt at optimising out a common factor for the calculation below it
                         temparg3 = aj1_gf(:,is,ie,il,igf)*vperp2_gf(:,ie,il,is)
                         antotp(:,it,ik) = antotp(:,it,ik) + wgt(is)*w(ie,is)*wl(:,il)*((temparg3*gfarray(:,1,is,ie,il,igf))+(temparg3*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          else
             do igf = gf_lo%llim_proc,gf_lo%ulim_proc
                it = it_idx(gf_lo,igf)
                ik = ik_idx(gf_lo,igf)
                antotp(:,it,ik) =0.
                do il = 1,gf_lo%nlambda
                   do ie = 1,gf_lo%negrid
                      do is = 1,gf_lo%nspec
                         temparg3 = aj1_gf(:,is,ie,il,igf)*vperp2_gf(:,ie,il,is)
                         antotp(:,it,ik) = antotp(:,it,ik) + wgt(is)*w(ie,is)*wl(:,il)*((temparg3*gfarray(:,1,is,ie,il,igf))+(temparg3*gfarray(:,2,is,ie,il,igf)))
                      end do
                   end do
                end do
             end do

          endif
       end if

    end if

    if (allocated(gfarray)) deallocate(gfarray)
  end subroutine getan_nogath_from_dfn

  !> Calls getan_nogath_from_dfn to calculate antot arrays for
  !> current gnew value
  subroutine getan_nogath (antot, antota, antotp, local_only)
    use dist_fn_arrays, only: gnew
    use theta_grid, only: ntgrid
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in out) :: antot, antota, antotp
    logical, intent(in), optional :: local_only
    call getan_nogath_from_dfn(gnew, antot, antota, antotp, local_only)
  end subroutine getan_nogath

!///////////////////////////////////////
!///////////////////////////////////////
!///////////////////////////////////////


  !> Calculates the potentials consistent with the passed nonadiabatic
  !> distribution function, `f_in`. Note this is not what GS2 usually
  !> evolves as `g`/`gnew` but the adjusted version.
  !>
  !> This is closely related to [[get_init_field]] which does the same job
  !> but works with the `g`/`gnew` modified distribution function instead.
  !>
  !> Currently we force potentials to zero if they are not included in the
  !> simulation. We should provide a method to calculate these fields even
  !> if they are not included in the evolution (e.g. for diagnostics,
  !> collisions etc.). Perhaps an optional flag here would be enough.
  !>
  !> @note We could precalculate the denominators used here saving a little
  !> computation and removing the need for a branch in the Apar case.
  subroutine calculate_potentials_from_nonadiabatic_dfn(f_in, phi, apar, bpar, &
       gf_lo, local_only)
    use gs2_layouts, only: g_lo
    use theta_grid, only: ntgrid, bmag
    use run_parameters, only: has_phi, has_apar, has_bpar, beta, tite
    use kt_grids, only: naky, ntheta0, inv_kperp2
    use species, only: spec, has_electron_species, has_ion_species
    use optionals, only: get_option_with_default
    use dist_fn_arrays, only: antot, antota, antotp
    implicit none
    complex, dimension(-ntgrid:, :, g_lo%llim_proc:), intent(in) :: f_in
    complex, dimension(-ntgrid:, :, :), intent(out) :: phi, apar, bpar
    logical, optional :: gf_lo, local_only
    real :: phi_weight
    real, dimension(-ntgrid:ntgrid) :: bpar_weight
    logical :: local_gf_lo, local_local_only
    integer :: it, ik, ig
    integer :: it_llim, it_ulim, ik_llim, ik_ulim

    local_gf_lo = get_option_with_default(gf_lo, .false.)
    local_local_only = get_option_with_default(local_only, .false.)

    ! Get the weighted velocity space integrals of the passed
    ! distribution function.
    if(local_gf_lo) then
       call getan_nogath_from_dfn (f_in, antot, antota, antotp, local_only)
    else
       call getan_from_dfn (f_in, antot, antota, antotp, local_only)
    end if

    if (local_local_only) then
       it_llim = g_lo%it_min ; it_ulim = g_lo%it_max
       ik_llim = g_lo%ik_min ; ik_ulim = g_lo%ik_max
    else
       it_llim = 1 ; it_ulim = ntheta0
       ik_llim = 1 ; ik_ulim = naky
    end if

    if (has_phi) then
       ! We could add in kperp2*poisfac to the weight to retain this feature
       if (.not. has_electron_species(spec) .or. .not. has_ion_species(spec) ) then
          phi_weight = 1.0 / (sum(spec%dens*spec%z*spec%z / spec%temp) + tite)
       else
          phi_weight = 1.0 / sum(spec%dens*spec%z*spec%z / spec%temp)
       end if

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, phi, antot, phi_weight) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             phi(:, it, ik) = antot(:, it, ik) * phi_weight
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, phi) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             phi(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (has_apar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik, ig) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, ntgrid, inv_kperp2, apar, antota) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             do ig = -ntgrid, ntgrid
                apar(ig, it, ik) = antota(ig, it, ik) * inv_kperp2(ig, it, ik)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, apar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             apar(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (has_bpar) then
       bpar_weight = -beta / bmag**2
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, bpar, bpar_weight, antotp) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             bpar(:, it, ik) = antotp(:, it, ik) * bpar_weight
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, bpar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             bpar(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine calculate_potentials_from_nonadiabatic_dfn

  !> Precomputes constants used in the inversion of the field equations:
  !>   gamtot * phi - gamtot1 * bpar = antot
  !>   kperp2 * apar = antota
  !>   beta/2 * gamtot1 * phi + (beta * gamtot2 + 1) * bpar = - beta * antotp
  !>
  !> @note As apar_denom is exactly inv_kperp2 we do not calculate this here
  !> and just use inv_kperp2 in get_fields_direct_from_dfn
  subroutine setup_get_fields_direct_from_dfn_denominators(phi_denom, bpar_denom)
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: naky, ntheta0
    use run_parameters, only: beta, has_phi, has_bpar
    use array_utils, only: copy
    implicit none
    real, dimension(:, :, :), allocatable, intent(in out) :: phi_denom, bpar_denom
    integer :: it, ik, ig

    ! get phi denominator
    if (has_phi) then
       ! Note we allocate and set as the full size to support calls to
       ! this routine with either choice for local_only.
       if (.not. allocated(phi_denom)) allocate(phi_denom(-ntgrid:ntgrid, ntheta0, naky))

       if (has_bpar) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ntheta0, naky, beta, gamtot2, &
          !$OMP bmag, gamtot1, phi_denom, gamtot) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = 1, naky
             do it = 1, ntheta0
                phi_denom(:, it, ik) = (beta * gamtot2(:, it, ik) + bmag**2) &
                     * gamtot(:, it, ik) &
                     + (beta/2.0) * gamtot1(:, it, ik) * gamtot1(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          ! If fbpar = 0 then the terms in earlier branch proportional to beta
          ! should be zero as these terms arrive from the bpar parts of quasineutrality
          ! as expressed in terms of gnew. In this case we get the below.
          call copy(gamtot, phi_denom)
       end if

       ! Now invert the denominator
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik, ig) &
       !$OMP SHARED(ntheta0, naky, ntgrid, phi_denom) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                if (abs(phi_denom(ig, it, ik)) > epsilon(0.0)) then
                   phi_denom(ig, it, ik) = &
                        1.0 / phi_denom(ig, it, ik)
                else
                   phi_denom(ig, it, ik) = 0.0
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

    ! get bpar
    if (has_bpar) then
       ! Note we allocate and set as the full size to support calls to
       ! this routine with either choice for local_only.
       if (.not. allocated(bpar_denom)) allocate(bpar_denom(-ntgrid:ntgrid, ntheta0, naky))

       ! As in our calculation of phi, we see that the equations for
       ! bpar and phi are coupled when written in terms of gnew. This
       ! means we get contributions arising from phi in the expression
       ! for B|| here. In particular, they are the terms involving
       ! gamtot1. However, if we have disabled phi then they should
       ! not appear in our calculation, so handle that here.
       if (has_phi) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ntheta0, naky, beta, gamtot, &
          !$OMP gamtot1, bpar_denom, gamtot2, bmag) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = 1, naky
             do it = 1, ntheta0
                bpar_denom(:, it, ik) = gamtot(:, it, ik) * &
                     (beta * gamtot2(:, it, ik) + bmag**2) + &
                     (beta/2.0) * gamtot1(:, it, ik) * gamtot1(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ntheta0, naky, beta, &
          !$OMP bpar_denom, gamtot2, bmag) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = 1, naky
             do it = 1, ntheta0
                bpar_denom(:, it, ik) = (beta * gamtot2(:, it, ik) + bmag**2)
             end do
          end do
          !$OMP END PARALLEL DO
       end if

       ! Now invert the denominator
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik, ig) &
       !$OMP SHARED(ntheta0, naky, ntgrid, bpar_denom) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = 1, naky
          do it = 1, ntheta0
             do ig = -ntgrid, ntgrid
                if (abs(bpar_denom(ig, it, ik)) > epsilon(0.0)) then
                   bpar_denom(ig, it, ik) = &
                        1.0 / bpar_denom(ig, it, ik)
                else
                   bpar_denom(ig, it, ik) = 0.0
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine setup_get_fields_direct_from_dfn_denominators

  !> Inverts the field equations:
  !>   gamtot * phi - gamtot1 * bpar = antot
  !>   kperp2 * apar = antota
  !>   beta/2 * gamtot1 * phi + (beta * gamtot2 + 1) * bpar = - beta * antotp
  !>
  !> @note I haven't made any check for use_Bpar=T case.
  !>
  !> TT> Given initial distribution function this obtains consistent fields
  !> MAB> ported from agk
  !> CMR, 1/8/2011> corrections below for inhomogeneous bmag
  subroutine get_fields_direct_from_dfn(g_in, phi, apar, bpar, gf_lo, local_only)
    use run_parameters, only: beta, has_phi, has_apar, has_bpar
    use species, only: spec, has_electron_species
    use theta_grid, only: ntgrid, bmag
    use kt_grids, only: ntheta0, naky, inv_kperp2
    use optionals, only: get_option_with_default
    use gs2_layouts, only: g_lo
    use dist_fn_arrays, only: antot, antota, antotp
    implicit none
    complex, dimension (:,:,:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (out) :: phi, apar, bpar
    logical, optional :: gf_lo, local_only
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky) :: numerator
    complex, dimension (ntheta0,naky) :: fl_avg
    logical :: local_gf_lo, local_local_only
    integer :: it, ik
    integer :: it_llim, it_ulim, ik_llim, ik_ulim

    local_gf_lo = get_option_with_default(gf_lo, .false.)

    local_local_only = get_option_with_default(local_only, .false.)

    if (local_local_only) then
       it_llim = g_lo%it_min ; it_ulim = g_lo%it_max
       ik_llim = g_lo%ik_min ; ik_ulim = g_lo%ik_max
    else
       it_llim = 1 ; it_ulim = ntheta0
       ik_llim = 1 ; ik_ulim = naky
    end if

    ! Ensure the denominator arrays are setup.
    if (.not. allocated(inv_phi_denominator_g)) then
       call setup_get_fields_direct_from_dfn_denominators(inv_phi_denominator_g, &
            inv_bpar_denominator_g)
    end if

    if(local_gf_lo) then
       call getan_nogath_from_dfn (g_in, antot, antota, antotp, local_only)
    else
       call getan_from_dfn (g_in, antot, antota, antotp, local_only)
    end if

    fl_avg = 0.

    ! Todo:
    ! This logical is a constant so should perhaps be calculated once during
    ! initialisiation and given a helpful name.
    if( .not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg ) then
      call calculate_flux_surface_average(fl_avg,antot)
    end if

    ! get phi
    if (has_phi) then

       !CMR, 1/8/2011:  bmag corrections here:
       !BSP 20/04/2021: Check if fbpar= 0.0 as beta terms should not be present in this case.
       if (has_bpar) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, numerator, beta, gamtot2, &
          !$OMP bmag, antot, fl_avg, gamtot1, antotp) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                numerator(:, it, ik) = (beta * gamtot2(:, it, ik) + bmag**2) * (antot(:, it, ik)+fl_avg(it, ik)) - (beta * gamtot1(:, it, ik)) * antotp(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          ! If fbpar = 0 then the terms in earlier branch proportional to beta
          ! should be zero as these terms arrive from the bpar parts of quasineutrality
          ! as expressed in terms of gnew. In this case we get the below.
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, numerator, antot, fl_avg) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                numerator(:, it, ik) = antot(:, it, ik)+fl_avg(it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       end if

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, inv_phi_denominator_g, &
       !$OMP phi, numerator) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             phi(:, it, ik) = numerator(:, it, ik) * inv_phi_denominator_g(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, phi) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             phi(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    ! get apar
    if (has_apar) then
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, antota, apar, inv_kperp2) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             apar(:, it, ik) = antota(:, it, ik) * inv_kperp2(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, apar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             apar(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    ! get bpar
    if (has_bpar) then
       ! As in our calculation of phi, we see that the equations for
       ! bpar and phi are coupled when written in terms of gnew. This
       ! means we get contributions arising from phi in the expression
       ! for B|| here. In particular, they are the terms involving
       ! gamtot1. However, if we have disabled phi then they should
       ! not appear in our calculation, so handle that here.
       if (has_phi) then
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, numerator, beta, gamtot, &
          !$OMP antotp, gamtot1, antot, bmag) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                numerator(:, it, ik) = - (beta * gamtot(:, it, ik)) * antotp(:, it, ik) - (beta/2.0) * gamtot1(:, it, ik) * antot(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO DEFAULT(none) &
          !$OMP PRIVATE(it, ik) &
          !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, numerator, beta, antotp) &
          !$OMP COLLAPSE(2) &
          !$OMP SCHEDULE(static)
          do ik = ik_llim, ik_ulim
             do it = it_llim, it_ulim
                numerator(:, it, ik) = - beta * antotp(:, it, ik)
             end do
          end do
          !$OMP END PARALLEL DO
       end if

       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, inv_bpar_denominator_g, &
       !$OMP numerator, bpar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             bpar(:, it, ik) = numerator(:, it, ik) * inv_bpar_denominator_g(:, it, ik)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP PRIVATE(it, ik) &
       !$OMP SHARED(ik_llim, ik_ulim, it_llim, it_ulim, bpar) &
       !$OMP COLLAPSE(2) &
       !$OMP SCHEDULE(static)
       do ik = ik_llim, ik_ulim
          do it = it_llim, it_ulim
             bpar(:, it, ik) = 0.0
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine get_fields_direct_from_dfn

  !> Calls get_fields_direct_from_dfn to solve for fields consistent with current gnew
  subroutine get_init_field (phi, apar, bpar, gf_lo)
    use dist_fn_arrays, only: gnew
    implicit none
    complex, dimension (:,:,:), intent (out) :: phi, apar, bpar
    logical, optional :: gf_lo
    call get_fields_direct_from_dfn(gnew, phi, apar, bpar, gf_lo)
  end subroutine get_init_field

  !> Calculates the potentials directly from the passed distribution
  !> function and compares the result to the passed
  !> potentials. Returns the maximum relative error for each
  !> potential.
  subroutine get_field_inconsistency(g_in, phi_in, apar_in, bpar_in, phi_err, apar_err, bpar_err)
    use run_parameters, only: has_phi, has_apar, has_bpar
    use array_utils, only: zero_array
    implicit none
    complex, dimension(:,:,:), intent(in) :: g_in, phi_in, apar_in, bpar_in
    real, intent(out) :: phi_err, apar_err, bpar_err
    complex, dimension(:,:,:), allocatable :: phi_new, apar_new, bpar_new
    real, dimension(:,:,:), allocatable :: err

    allocate(err, mold = real(phi_in))
    allocate(phi_new, mold = phi_in); call zero_array(phi_new)
    allocate(apar_new, mold = phi_in); call zero_array(apar_new)
    allocate(bpar_new, mold = phi_in); call zero_array(bpar_new)

    phi_err = 0. ; apar_err = 0. ; bpar_err = 0.

    call get_fields_direct_from_dfn(g_in, phi_new, apar_new, bpar_new)

    if (has_phi) phi_err = get_field_relative_error(phi_in, phi_new)
    if (has_apar) apar_err = get_field_relative_error(apar_in, apar_new)
    if (has_bpar) bpar_err = get_field_relative_error(bpar_in, bpar_new)

  contains

    real function get_field_relative_error(field_in, field_new) result(rel_err)
      use mp, only: max_allreduce
      implicit none
      complex, dimension(:, :, :), intent(in) :: field_in, field_new
      err = 0.
      where(abs(field_in) > 0)
         err = abs(field_in - field_new)/abs(field_in)
      end where
      rel_err = maxval(err)
      call max_allreduce(rel_err)
    end function get_field_relative_error

  end subroutine get_field_inconsistency

  !> Routine to dump the current source term used in [[invert_rhs]].
  !>
  !> This might be relatively expensive so care should be taken in
  !> calling this.  It makes use of the existing code for saving
  !> restart files and as such produces either one file per processor
  !> or a single file depending on if the build uses parallel i/o or
  !> not. These files are automatically overwritten so calling it
  !> repeatedly within a single run might not be useful.
  subroutine dump_current_source_term(istep, phi, apar, bpar, phinew, &
       aparnew, bparnew, sourcefac)
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: g_work
    use theta_grid, only: ntgrid
    use run_parameters, only: has_phi, has_apar, has_bpar
    use gs2_time, only: user_time, code_dt, code_dt_prev1, code_dt_prev2, code_dt_max
    use collisions, only: vnmult
    use array_utils, only: zero_array
    implicit none
    integer, intent(in) :: istep
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi,    apar,    bpar
    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, aparnew, bparnew
    complex, intent(in) :: sourcefac
    character(len = 20), parameter :: extension = '.source'
    integer :: istatus

    ! Is this needed?
    call zero_array(g_work)

    call get_source_term_on_local_domain(g_work, phi, apar, bpar, phinew, &
         aparnew, bparnew, istep, sourcefac)

    call gs2_save_for_restart (g_work, user_time, vnmult, istatus, &
         has_phi, has_apar, has_bpar, &
         code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, fileopt = extension)
  end subroutine dump_current_source_term

  !> Routine to dump the homogeneous solution
  subroutine dump_homogeneous_solution(force_passing)
    use gs2_save, only: gs2_save_for_restart
    use dist_fn_arrays, only: g_work
    use run_parameters, only: has_phi, has_apar, has_bpar
    use gs2_time, only: user_time, code_dt, code_dt_prev1, code_dt_prev2, code_dt_max
    use collisions, only: vnmult
    implicit none
    logical, intent(in), optional :: force_passing
    character(len = 20), parameter :: extension = '.homogeneous'
    integer :: istatus
    call init_homogeneous_g(g_work, force_passing)
    call gs2_save_for_restart (g_work, user_time, vnmult, istatus, &
         has_phi, has_apar, has_bpar, &
         code_dt, code_dt_prev1, code_dt_prev2, code_dt_max, fileopt = extension)
  end subroutine dump_homogeneous_solution

  !> Identical to get_flux except don't integrate over poloidal angle ! JRB
  subroutine get_flux_dist (g_in, fld, flx_dist)
    use theta_grid, only: ntgrid, grho, theta, bmag, gradpar
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use integration, only: trapezoidal_integration
    use gs2_time, only: woutunits
    use gs2_layouts, only: g_lo
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    !> Input weighted distribution    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (-ntgrid:,:,:,:), intent (in out) :: flx_dist
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    ! EGH added final parameter 'all'
    ! to the call below for new parallel output
    ! This is temporary until distributed fields
    ! are implemented 1/2014
    call integrate_moment (g_in, total, moment_to_allprocs, full_arr)

    !Note that jacob = 1/(bmag*gradpar*drhodpsi). This is similar to the factor
    !we have below without the drhodpsi.
    wgt = 1.0/trapezoidal_integration(theta, grho / (bmag * gradpar))

    ! EGH for new parallel I/O everyone
    ! calculates fluxes
    ! Note: This doesn't include the field line jacobian (jacob) and as such
    ! it's not clear that the above wgt is correct (as it uses 1/(bmag*gradpar)
    ! which is jacob * drhodpsi, rather than jacob as one might expect).
    !if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                flx_dist(:,it,ik,is) = aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *aky(ik) * wgt / woutunits(ik)
             end do
          end do
       end do
       flx_dist = flx_dist*0.5
    ! end if

    deallocate (total)

  end subroutine get_flux_dist

#ifdef LOWFLOW
  !> Calculate toroidal angular momentum flux due to field
  !>
  !> FIXME: confirm & expand
  subroutine get_flux_tormom (g_in, fld, flx)
    use theta_grid, only: ntgrid
    use theta_grid, only : theta, grho, bmag, jacob
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    use kt_grids, only: aky
    use lowflow, only: mach_lab
    use integration, only trapezoidal_integration
    use gs2_layouts, only: g_lo
    implicit none
    !> Input weighted distribution    
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    !> Input field
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    !> Output flux
    real, dimension (:,:,:), intent (out) :: flx
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is
    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g_in, total)

    ! This is essentially 2 pi / surfarea
    wgt = 1.0/trapezoidal_integration(theta, grho * jacob)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                flx(it,ik,is) = trapezoidal_integration(theta, aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     *Rplot(:)**2 * jacob)*aky(ik)*wgt*mach_lab(is)
             end do
          end do
       end do
       flx = flx*0.5
    end if

    deallocate (total)
  end subroutine get_flux_tormom !JPL
#endif

  !> Calculate the momentum flux as a function of \((v_\parallel, \theta, t)\)
  subroutine flux_vs_theta_vs_vpa (phinew,vflx)
    use constants, only: zi
    use dist_fn_arrays, only: gnew, vperp2, aj1, aj0, vpac, g_work
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use theta_grid, only: ntgrid, bmag, gds21, gds22, qval, shat
    use theta_grid, only: Rplot, Bpol, rhoc
    use kt_grids, only: aky, theta0
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only: spec, nspec

    implicit none
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew
    real, dimension (-ntgrid:,:,:), intent (out) :: vflx
    
    integer :: all = 1
    integer :: iglo, isgn, ig, it, ik, is

    real, dimension (:,:,:), allocatable :: g0r
    real, dimension (:,:,:,:,:), allocatable :: gavg

    allocate (g0r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gavg(-ntgrid:ntgrid,nlambda,negrid,2,nspec))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             ! Not clear that this should use the grid-centre v|| (vpac)
             ! as this is not defined at ig = ntgrid and all other quantities
             ! are evaluated on the grid points.
             g_work(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)*vpac(ig,isgn,iglo) &
                  *Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2) &
                  -zi*aky(ik)*gnew(ig,isgn,iglo)*aj1(ig,iglo) &
                  *rhoc*(gds21(ig)+theta0(it,ik)*gds22(ig))*vperp2(ig,iglo)*spec(is)%smz/(qval*shat*bmag(ig)**2)
             g0r(ig,isgn,iglo) = aimag(g_work(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)
          end do
       end do
    end do
    
    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)

  end subroutine flux_vs_theta_vs_vpa

  !> Diagnose particle flux contribution to toroidal momentum flux
  !> in the lab frame in terms of vpar and theta.
  subroutine pflux_vs_theta_vs_vpa (vflx)
#ifdef LOWFLOW   
    use dist_fn_arrays, only: gnew, aj0, g_work
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use theta_grid, only: Rplot !JPL
    use fields_arrays, only: phinew
    use kt_grids, only: aky
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only:  nspec
    use lowflow, only: mach_lab    
#endif
    use theta_grid, only: ntgrid
    implicit none

    real, dimension (-ntgrid:,:,:), intent (out) :: vflx
#ifdef LOWFLOW  
    integer :: all = 1
    integer :: iglo, isgn, ig, it, ik, is

    real, dimension (:,:,:), allocatable :: g0r
    real, dimension (:,:,:,:,:), allocatable :: gavg
    allocate (g0r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gavg(-ntgrid:ntgrid,nlambda,negrid,2,nspec))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g_work(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo) 

             g0r(ig,isgn,iglo) = aimag(g_work(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)*Rplot(ig)**2*mach_lab(is)

          end do
       end do
    end do
    
    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)
#else
    vflx=0.
#endif
  end subroutine pflux_vs_theta_vs_vpa

  !> Diagnose the poloidal distribution of the particle, angular momentum, 
  !! and energy fluxes
  !!
  !! JRB
  subroutine flux_dist (phi, pflux_dist, vflux_par_dist, vflux_perp_dist, qflux_dist)

!CMR, 15/1/08: 
!  Implemented Clemente Angioni's fix for fluxes by replacing g with gnew 
!  so fields and distribution function are evaluated self-consistently in time.
!  This fixed unphysical oscillations in non-ambipolar particle fluxes 
!
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, aj1, vperp2, g_work
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use run_parameters, only:has_phi
    use constants, only: zi
    use theta_grid, only: Rplot, Bpol, rhoc
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (-ntgrid:,:,:,:), intent (out) :: pflux_dist, vflux_par_dist, vflux_perp_dist, qflux_dist
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    call zero_array(pflux_dist) ; call zero_array(vflux_par_dist) ; call zero_array(vflux_perp_dist) ; call zero_array(qflux_dist)

    if (has_phi) then
       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do
       call get_flux_dist (g_work, phi, pflux_dist)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
       call get_flux_dist (g_work, phi, qflux_dist)

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             ! Not clear that this should use the grid-centre v|| (vpac)
             ! as this is not defined at ig = ntgrid and all other quantities
             ! are evaluated on the grid points.
             g_work(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do
       call get_flux_dist (g_work, phi, vflux_par_dist)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) = -zi*aky(ik)*gnew(:,isgn,iglo)*aj1(:,iglo) &
                  *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,iglo)*spec(is)%smz/(qval*shat*bmag**2)
          end do
       end do
       call get_flux_dist (g_work, phi, vflux_perp_dist)
    end if

 end subroutine flux_dist


!=============================================================================
! Density: Calculate Density perturbations
 !=============================================================================
 
  !> Calculate some heating quantities:
  !> - ion/electron heating
  !> - antenna power and B-field contributions to E and E_dot
  !> - gradient contributions to heating
  !> - hyperviscosity
  !> - hyperresistivity
  !>
  !> FIXME: This is too large and does too much, split into smaller routines
  subroutine get_heat (h, hk, phi, apar, bpar, phinew, aparnew, bparnew)
    use mp, only: proc0
    use constants, only: zi
    use kt_grids, only: ntheta0, naky, aky, akx, kperp2
#ifdef LOWFLOW
    use lowflow, only: hneoc
#endif
    use dist_fn_arrays, only: vpac, aj0, aj1, vperp2, g, gnew, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_heating, only: heating_diagnostics
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, ie_idx
    use le_grids, only: integrate_moment
    use species, only: spec, nspec,has_electron_species
    use theta_grid, only: jacob, delthet, ntgrid
    use run_parameters, only: fphi, fapar, fbpar, beta, tite
    use gs2_time, only: code_dt, tunits
    use nonlinear_terms, only: nonlin
    use antenna, only: antenna_apar, a_ext_data
    use hyper, only: D_v, D_eta, nexp, hypervisc_filter

    implicit none
    type (heating_diagnostics), intent(in out) :: h
    type (heating_diagnostics), dimension(:,:), intent(in out) :: hk
!    complex, dimension (-ntgrid:,:,:), pointer :: hh, hnew
    complex, dimension (-ntgrid:,:,:), intent(in) :: phi, apar, bpar, phinew, aparnew, bparnew
    complex, dimension(:,:,:,:), allocatable :: tot
!    complex, dimension (:,:,:), allocatable :: epar
    complex, dimension(:,:,:), allocatable :: bpardot, apardot, phidot, j_ext
    complex :: chi, havg
    complex :: chidot, j0phiavg, j1bparavg, j0aparavg
!    complex :: pstar, pstardot, gdot
    complex :: phi_m, apar_m, bpar_m, hdot
    complex :: phi_avg, bpar_avg, bperp_m, bperp_avg
    complex :: de, denew
    complex :: dgdt_hypervisc
    real, dimension (:), allocatable :: wgt
    real :: fac2, dtinv, akperp4
    integer :: isgn, iglo, ig, is, ik, it, ie

    g_work(ntgrid,:,:) = 0.

! ==========================================================================
! Ion/Electron heating------------------------------------------------------
! ==========================================================================

    allocate ( phidot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (apardot(-ntgrid:ntgrid, ntheta0, naky))
    allocate (bpardot(-ntgrid:ntgrid, ntheta0, naky))

    call dot ( phi,  phinew,  phidot, fphi)
    call dot (apar, aparnew, apardot, fapar)
    call dot (bpar, bparnew, bpardot, fbpar)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Next two calls make the variables g, gnew = h, hnew 
!!! until the end of this procedure!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call g_adjust (g,    phi,    bpar,    direction = from_g_gs2)
    call g_adjust (gnew, phinew, bparnew, direction = from_g_gs2)

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2
          
          do ig=-ntgrid, ntgrid-1
             
             chidot = aj0(ig,iglo)*(phidot(ig,it,ik) &
                  - vpac(ig,isgn,iglo) * spec(is)%stm * apardot(ig,it,ik)) &
                  + aj1(ig,iglo)*2.0*vperp2(ig,iglo)*bpardot(ig,it,ik)*spec(is)%tz
             
             hdot = fdot (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo), dtinv)
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))
             
! First term on RHS and LHS of Eq B-10 of H1:

             g_work(ig,isgn,iglo) = spec(is)%dens*conjg(havg)* &
                  (chidot*spec(is)%z-hdot*spec(is)%temp)

          end do
       end do
    end do

    deallocate (phidot, apardot, bpardot)

    allocate (tot(-ntgrid:ntgrid, ntheta0, naky, nspec))

    call integrate_moment (g_work, tot)

    if (proc0) then
       allocate (wgt(-ntgrid:ntgrid))
       wgt = 0.
       do ig=-ntgrid,ntgrid-1
!          wgt(ig) = delthet(ig)*jacob(ig)
! delthet is cell-centered, but jacob is given on grid
          wgt(ig) = delthet(ig)*(jacob(ig)+jacob(ig+1))*0.5
       end do
       wgt = wgt/sum(wgt)         

       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                    hk(it,ik) % heating(is) = hk(it,ik) % heating(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2 
                end do
                h % heating(is) = h % heating(is) + hk(it,ik) % heating(is)
             end do
          end do
       end do
    end if

! ==========================================================================
! Antenna Power and B-field contribution to E and E_dot---------------------
! ==========================================================================
    if (proc0) then
!       allocate (epar (-ntgrid:ntgrid, ntheta0, naky)) ; epar = 0.
       allocate (j_ext(-ntgrid:ntgrid, ntheta0, naky)) ; j_ext = 0.
       call antenna_apar (kperp2, j_ext)       
       
       if (beta > epsilon(0.)) then
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0
                
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                
                do ig=-ntgrid, ntgrid-1
                   
                   !GGH Time and space averaged estimate of d/dt(apar)
                   apar_m = fdot (apar   (ig  ,it,ik), &
                                  apar   (ig+1,it,ik), &
                                  aparnew(ig  ,it,ik), &
                                  aparnew(ig+1,it,ik), dtinv)*fapar

! J_ext.E when driving antenna only includes A_parallel:

                   hk(it,ik) % antenna = hk(it, ik) % antenna + real(conjg(j_ext(ig,it,ik))*apar_m)*wgt(ig)*fac2

                   !GGH Time and space averaged estimate of d/dt(bperp)
                   bperp_m = fdot (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                   aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                   aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), dtinv) * fapar

                   !GGH Time and space averaged estimate of d/dt(bpar)
                   bpar_m = fdot (bpar   (ig  ,it,ik), &
                                  bpar   (ig+1,it,ik), &
                                  bparnew(ig  ,it,ik), &
                                  bparnew(ig+1,it,ik), dtinv)*fbpar

                   !GGH Time and space averaged estimate of bperp
                   bperp_avg = favg (apar   (ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     apar   (ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik)), &
                                     aparnew(ig  ,it,ik)*sqrt(kperp2(ig  ,it,ik)), &
                                     aparnew(ig+1,it,ik)*sqrt(kperp2(ig+1,it,ik))) * fapar

                   !GGH Time and space averaged estimate of bpar
                   bpar_avg = favg (bpar   (ig  ,it,ik), &
                                    bpar   (ig+1,it,ik), &
                                    bparnew(ig  ,it,ik), &
                                    bparnew(ig+1,it,ik)) * fbpar

! 1/2 * d/dt B**2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        real(0.25 * conjg(bperp_m)*bperp_avg + conjg(bpar_m)*bpar_avg) &
                        * wgt(ig)*fac2*(2.0/beta)

! B**2/2
!! GGH: Bug fixed on 2/06; error was in relative weight of B_par**2 and B_perp**2   
                   hk(it,ik) % energy = hk(it,ik) % energy &
                        + 0.5*real((0.25*conjg(bperp_avg)*bperp_avg + conjg(bpar_avg)*bpar_avg)) &
                        * wgt(ig)*fac2*(2.0/beta)

                   !Eapar = int k_perp^2 A_par^2/(8 pi)                   
                   hk(it,ik) % eapar = hk(it,ik) % eapar &
                        + 0.5*real(0.25*conjg(bperp_avg)*bperp_avg) * wgt(ig)*fac2*(2.0/beta)

                   !Ebpar = int B_par^2/(8 pi)
                   hk(it,ik) % ebpar = hk(it,ik) % ebpar &
                        + 0.5*real(conjg(bpar_avg)*bpar_avg) * wgt(ig)*fac2*(2.0/beta)
                   
                end do
                h % antenna = h % antenna + hk(it, ik) % antenna
                h % eapar = h % eapar + hk(it, ik) % eapar
                h % ebpar = h % ebpar + hk(it, ik) % ebpar
             end do
          end do
       else
          hk % antenna = 0.
          h  % antenna = 0.
          hk % energy_dot = 0.
          hk % energy = 0.
          hk % eapar = 0.
          h  % eapar = 0.
          hk % ebpar = 0.
          h  % ebpar = 0.
       end if
       deallocate (j_ext)
    end if

! ==========================================================================
! Finish E_dot--------------------------------------------------------------
! ==========================================================================

!GGH Include response of Boltzmann species for single-species runs

    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   phi_m   = fdot (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik), dtinv)

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot + &
                        fphi * real(conjg(phi_avg)*phi_m) &
                        * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to E_dot for single species runs---------------------
 
!GGH New E_dot calc
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       dtinv = 1./(code_dt*tunits(ik))
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             !Calculate old fluctuating energy de
             havg = favg_x (g(ig  ,isgn,iglo), &
                            g(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phi(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phi(ig  ,it,ik), &
                               phi(ig+1,it,ik)) * fphi * spec(is)%zt

             de=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)

            !Calculate new fluctuating energy denew
             havg = favg_x (gnew(ig  ,isgn,iglo), &
                            gnew(ig+1,isgn,iglo))

             j0phiavg = favg_x (aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg_x (phinew(ig  ,it,ik), &
                               phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             denew=0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg) 

             !Set g_work as the change of energy (denew-de)/dt
             g_work(ig,isgn,iglo) = fdot_t(de,denew,dtinv)

          end do
       end do
    end do
    !GGH -END new e_dot calc

   call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy_dot = hk(it,ik) % energy_dot  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy_dot = h % energy_dot + hk(it,ik) % energy_dot
          end do
       end do
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ==========================================================================
! Gradient Contributions to Heating-----------------------------------------
! ==========================================================================

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       ie = ie_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle

       do isgn=1,2
          do ig=-ntgrid, ntgrid-1

             chi = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phi   (ig+1,it,ik),  &
                         aj0(ig  ,iglo)*phinew(ig  ,it,ik),  &
                         aj0(ig+1,iglo)*phinew(ig+1,it,ik)) &
                         *fphi

!!GGH Bug fix: The apar part should be subtracted (because chi= phi - v|| A|| + B||)
             chi = chi - &
                  favg (aj0(ig  ,iglo)*apar   (ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*apar   (ig+1,it,ik)*vpac(ig+1,isgn,iglo),  &
                        aj0(ig  ,iglo)*aparnew(ig  ,it,ik)*vpac(ig  ,isgn,iglo),  &
                        aj0(ig+1,iglo)*aparnew(ig+1,it,ik)*vpac(ig+1,isgn,iglo)) &
                        *spec(is)%stm*fapar
                
             chi = chi + &
                  favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                        aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                        aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                        *spec(is)%tz*fbpar

             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

#ifdef LOWFLOW
             g_work(ig,isgn,iglo) = zi * wstar(ik,ie,is)*hneoc(ig,isgn,iglo)/code_dt * conjg(havg)*chi &
#else
             g_work(ig,isgn,iglo) = zi * wstar(ik,ie,is)/code_dt * conjg(havg)*chi &
#endif
                  * spec(is)%dens
            
          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % gradients(is) = hk(it,ik) % gradients(is) &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
                h % gradients(is) = h % gradients(is) + hk(it,ik) % gradients(is)
             end do
          end do
       end do
    end if
! ==========================================================================
! Hyperviscosity------------------------------------------------------------
! ==========================================================================

    if (D_v > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
!          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0phiavg = favg (aj0(ig  ,iglo) * phi(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phi(ig+1,it,ik), &
                                 aj0(ig  ,iglo) * phinew(ig  ,it,ik), &
                                 aj0(ig+1,iglo) * phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                j1bparavg= favg (aj1(ig  ,iglo)*2.0*bpar   (ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bpar   (ig+1,it,ik)*vperp2(ig+1,iglo),  &
                                 aj1(ig  ,iglo)*2.0*bparnew(ig  ,it,ik)*vperp2(ig  ,iglo),  &
                                 aj1(ig+1,iglo)*2.0*bparnew(ig+1,it,ik)*vperp2(ig+1,iglo)) &
                                 *fbpar

                dgdt_hypervisc = 0.5*((1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo) &
                     + (1.0-1./hypervisc_filter(ig+1,it,ik))*gnew(ig+1,isgn,iglo))/code_dt

!Set g_work for hyperviscous heating
!                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_v*akperp4* &
!                     ( conjg(havg)*havg - conjg(havg)*j0phiavg - &
!                     conjg(havg)*j1bparavg)
                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*conjg(havg)*dgdt_hypervisc

             end do
          end do
       end do

       call integrate_moment (g_work, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hypervisc(is) = hk(it,ik) % hypervisc(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hypervisc(is) = h % hypervisc(is) + hk(it,ik) % hypervisc(is)
                end do
             end do
          end do
       end if

    end if !End Hyperviscous Heating Calculation


! ==========================================================================
! Hyperresistivity------------------------------------------------------------
! ==========================================================================
 
    if (D_eta > epsilon(0.)) then

       do iglo=g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo, iglo)
          it = it_idx(g_lo, iglo)
          ik = ik_idx(g_lo, iglo)
          if (nonlin .and. it == 1 .and. ik == 1) cycle
          akperp4 = (aky(ik)**2 + akx(it)**2)**nexp
           do isgn=1,2
             do ig=-ntgrid, ntgrid-1
                
                havg = favg (g   (ig  ,isgn,iglo), &
                             g   (ig+1,isgn,iglo), &
                             gnew(ig  ,isgn,iglo), &
                             gnew(ig+1,isgn,iglo)) 

                j0aparavg = favg (aj0(ig  ,iglo) * apar(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * apar(ig+1,it,ik), &
                                 aj0(ig  ,iglo)  * aparnew(ig  ,it,ik), &
                                 aj0(ig+1,iglo)  * aparnew(ig+1,it,ik)) & 
                                 * fapar * spec(is)%zstm * vpac(ig,isgn,iglo)

!Set g_work for hyperresistive heating
                g_work(ig,isgn,iglo) = spec(is)%dens*spec(is)%temp*D_eta*akperp4* &
                     conjg(havg)*j0aparavg

             end do
          end do
       end do

       call integrate_moment (g_work, tot)
       if (proc0) then
          do ik = 1, naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             do it = 1, ntheta0
                if (nonlin .and. it == 1 .and. ik == 1) cycle
                do is = 1, nspec
                   do ig = -ntgrid, ntgrid-1
                      hk(it,ik) % hyperres(is) = hk(it,ik) % hyperres(is) &
                           + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                   end do
                   h % hyperres(is) = h % hyperres(is) + hk(it,ik) % hyperres(is)
                end do
             end do
          end do
       end if

    end if !End Hyperresistivity Heating Calculation

!==========================================================================
!Finish Energy-------------------------------------------------------------
!==========================================================================

!GGH Calculate hs2-------------------------------------------------------------
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             g_work(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg)
          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % hs2(is) = hk(it,ik) % hs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                end do
             end do
             h % hs2(:) = h % hs2(:) + hk(it,ik) % hs2(:)
          end do
       end do
    end if

!Calculate phis2-------------------------------------------------------------
    if (proc0) then
       do ik=1,naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1,ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do ig=-ntgrid, ntgrid-1
                do is = 1, nspec
                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik)) * fphi * spec(is)%zt

                   !hs2 = int_r int_v T/F0 hs^2/2
                   hk(it,ik) % phis2(is) = hk(it,ik) % phis2(is)  &
                        +0.5*spec(is)%temp*spec(is)%dens*real(conjg(phi_avg)*phi_avg) &
                        * wgt(ig) * fac2
                enddo
             end do
             h % phis2(:) = h % phis2(:) + hk(it,ik) % phis2(:)
          end do
       end do
    endif

! Calculate delfs2 (rest of energy)-----------------------------------------------

!GGH  Include response of Boltzmann species for single species runs
    if (.not. has_electron_species(spec)) then
       if (proc0) then
          !NOTE: It is assumed here that n0i=n0e and zi=-ze
          do ik=1,naky
             fac2 = 0.5
             if (aky(ik) < epsilon(0.0)) fac2 = 1.0
             dtinv = 1./(code_dt*tunits(ik))
             do it = 1,ntheta0

                if (nonlin .and. it == 1 .and. ik == 1) cycle

                do ig=-ntgrid, ntgrid-1

                   phi_avg = favg (phi   (ig  ,it,ik), &
                        phi   (ig+1,it,ik), &
                        phinew(ig  ,it,ik), &
                        phinew(ig+1,it,ik))

                   !NOTE: Adiabatic (Boltzmann) species has temperature
                   !       T = spec(1)%temp/tite
                   hk(it,ik) % energy = hk(it,ik) % energy + &
                        fphi * real(conjg(phi_avg)*phi_avg) &
                        * 0.5 * spec(1)%dens * spec(1)%z * spec(1)%z * (tite/spec(1)%temp) &
                        * wgt(ig)*fac2

                end do
             end do
          end do
       endif
    endif !END Correction to energy for single species runs---------------------

    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       is = is_idx(g_lo, iglo)
       it = it_idx(g_lo, iglo)
       ik = ik_idx(g_lo, iglo)
       if (nonlin .and. it == 1 .and. ik == 1) cycle
       do isgn=1,2

          do ig=-ntgrid, ntgrid-1
             
             havg = favg (g   (ig  ,isgn,iglo), &
                          g   (ig+1,isgn,iglo), &
                          gnew(ig  ,isgn,iglo), &
                          gnew(ig+1,isgn,iglo))

             j0phiavg = favg (aj0(ig  ,iglo)*phi   (ig  ,it,ik), &
                              aj0(ig+1,iglo)*phi   (ig+1,it,ik), &  
                              aj0(ig  ,iglo)*phinew(ig  ,it,ik), &  
                              aj0(ig+1,iglo)*phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             phi_avg = favg (phi   (ig  ,it,ik), &
                             phi   (ig+1,it,ik), &
                             phinew(ig  ,it,ik), &
                             phinew(ig+1,it,ik)) * fphi * spec(is)%zt

             g_work(ig,isgn,iglo) = 0.5*spec(is)%temp*spec(is)%dens*(conjg(havg)*havg &
                  + conjg(phi_avg)*phi_avg &
                  - conjg(j0phiavg)*havg &
                  - conjg(havg)*j0phiavg)
          end do
       end do
    end do

    call integrate_moment (g_work, tot)

    if (proc0) then
       do ik = 1, naky
          fac2 = 0.5
          if (aky(ik) < epsilon(0.0)) fac2 = 1.0
          do it = 1, ntheta0
             if (nonlin .and. it == 1 .and. ik == 1) cycle
             do is = 1, nspec             
                do ig = -ntgrid, ntgrid-1
                   hk(it,ik) % energy = hk(it,ik) % energy  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2

                   !Delfs2 = int_r int_v T/F0 dfs^2/2
                   hk(it,ik) % delfs2(is) = hk(it,ik) % delfs2(is)  &
                        + real(tot(ig,it,ik,is))*wgt(ig)*fac2
                end do
             end do
             h % energy = h % energy + hk(it,ik) % energy
             h % delfs2(:) = h % delfs2(:) + hk(it,ik) % delfs2(:)
          end do
       end do
       deallocate (wgt)
    end if

    deallocate (tot)

!!
!! Put g, gnew back to their usual meanings
!!
    call g_adjust (g,    phi,    bpar,    direction = to_g_gs2)
    call g_adjust (gnew, phinew, bparnew, direction = to_g_gs2)

  end subroutine get_heat

  !> Error estimate obtained by comparing standard integral with less-accurate integral
  !>
  !> Estimate of the (1) absolute and (2) relative errors resulting from velocity space
  !> integrals in the calculation of the following quantities in the given dimensions: (1)
  !> k phi, energy (2) k phi, untrapped pitch angles (3) k phi, trapped pitch angles, (4)
  !> k apar, energy, (5) k apar, untrapped angles. Relative errors should be < 0.1.
  subroutine get_verr (errest, erridx, phi, bpar)
    use mp, only: broadcast
    use le_grids, only: integrate_species
    use le_grids, only: eint_error, trap_error, lint_error, wdim
    use le_grids, only: ng2, nlambda, new_trap_int, grid_has_trapped_particles
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky, aky, akx
    use species, only: nspec, spec
    use dist_fn_arrays, only: gnew, aj0, vpa, g_adjust, g_work, to_g_gs2, from_g_gs2
    use run_parameters, only: has_phi, has_apar, beta
    use gs2_layouts, only: g_lo
    use collisions, only: init_lorentz, init_ediffuse, init_lorentz_conserve
    use collisions, only: etol, ewindow, etola, ewindowa, init_diffuse_conserve
    use collisions, only: vnmult, vary_vnew
    use array_utils, only: zero_array
    implicit none

    integer :: ig, it, ik, iglo, isgn, ntrap

    !> Indices of maximum error.
    !> FIXME: Make this fixed size
    integer, dimension (:,:), intent (out) :: erridx
    !> Estimated error.
    !> FIXME: Make this fixed size.
    real, dimension (:,:), intent (out) :: errest
    !> Electrostatic potential and parallel magnetic field
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar

    complex, dimension (:,:,:), allocatable :: phi_app, apar_app
    complex, dimension (:,:,:,:), allocatable :: phi_e, phi_l, phi_t
    complex, dimension (:,:,:,:), allocatable :: apar_e, apar_l, apar_t
    real, dimension (:), allocatable :: wgt, errtmp
    integer, dimension (:), allocatable :: idxtmp

    real :: vnmult_target
    logical :: trap_flag=.true. !Value doesn't really matter as just used in present() checks

    allocate(wgt(nspec))
    allocate(errtmp(2))
    allocate(idxtmp(3))

    !This should probably be something like if(first) then ; call broadcast ; first=.false. ; endif
    !as we only deallocate kmax during finish_dist_fn.
    !Really wdim should just be broadcast when we first calculate it -- currently slightly complicated
    !as calculation triggered in gs2_diagnostics and only called by proc0.
    if (.not. allocated(kmax)) call broadcast (wdim)  ! Odd-looking statement.  BD

    if (has_phi) then
       allocate(phi_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(phi_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(phi_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

    if (has_apar) then
       allocate(apar_app(-ntgrid:ntgrid,ntheta0,naky))
       allocate(apar_e(-ntgrid:ntgrid,ntheta0,naky,wdim))
       allocate(apar_l(-ntgrid:ntgrid,ntheta0,naky,ng2))
    end if

! first call to g_adjust converts gyro-averaged dist. fn. (g)
! into nonadiabatic part of dist. fn. (h)

    call g_adjust (gnew, phi, bpar, direction = from_g_gs2)

! take gyro-average of h at fixed total position (not g.c. position)
    if (has_phi) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          ! TEMP FOR TESTING -- MAB
!          il = il_idx(g_lo,iglo)
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                ! TEMP FOR TESTING -- MAB
                g_work(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
!                g_work(ig,isgn,iglo) = sqrt(max(0.0,1.0-bmag(ig)*al(il)))
!                g_work = 1.
             end do
          end do
       end do

       wgt = spec%z*spec%dens

       call integrate_species (g_work, wgt, phi_app)

! integrates dist fn of each species over v-space
! after dropping an energy grid point and returns
! phi_e, which contains the integral approximations
! to phi for each point dropped

       call eint_error (g_work, wgt, phi_e)

! integrates dist fn of each species over v-space
! after dropping an untrapped lambda grid point and returns phi_l.
! phi_l contains ng2 approximations for the integral over lambda that
! come from dropping different pts from the gaussian quadrature grid

       call lint_error (g_work, wgt, phi_l)

! next loop gets error estimate for trapped particles, if there are any

       if (grid_has_trapped_particles() .and. new_trap_int) then
          ntrap = nlambda - ng2
          allocate(phi_t(-ntgrid:ntgrid,ntheta0,naky,ntrap))       
          call zero_array(phi_t)
          call trap_error (g_work, wgt, phi_t)
       end if

    end if

    if (has_apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             do ig=-ntgrid, ntgrid
                g_work(ig,isgn,iglo) = aj0(ig,iglo)*vpa(ig,isgn,iglo)*gnew(ig,isgn,iglo)
             end do
          end do
       end do
       

       wgt = 2.0*beta*spec%z*spec%dens*sqrt(spec%temp/spec%mass)
       call integrate_species (g_work, wgt, apar_app)

       call eint_error (g_work, wgt, apar_e)
       call lint_error (g_work, wgt, apar_l)
       if (grid_has_trapped_particles() .and. new_trap_int) then
          ntrap = nlambda - ng2
          allocate(apar_t(-ntgrid:ntgrid,ntheta0,naky,ntrap))       
          call zero_array(apar_t)
          call trap_error (g_work, wgt, apar_t)
       end if
       
    end if

! second call to g_adjust converts from h back to g

    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    if (.not. allocated(kmax)) then
       allocate (kmax(ntheta0, naky))
       do ik = 1, naky
          do it = 1, ntheta0
             kmax(it,ik) = max(akx(it),aky(ik))
          end do
       end do
    end if
    
    errest = 0.0
    erridx = 0
    
    if (has_phi) then

       call estimate_error (phi_app, phi_e, kmax, errtmp, idxtmp)
       errest(1,:) = errtmp
       erridx(1,:) = idxtmp
       
       call estimate_error (phi_app, phi_l, kmax, errtmp, idxtmp)
       errest(2,:) = errtmp
       erridx(2,:) = idxtmp

       if (grid_has_trapped_particles() .and. new_trap_int) then
          call estimate_error (phi_app, phi_t, kmax, errtmp, idxtmp, trap_flag)
          errest(3,:) = errtmp
          erridx(3,:) = idxtmp
          deallocate (phi_t)
       end if

    end if
    
    if (has_apar) then

       call estimate_error (apar_app, apar_e, kmax, errtmp, idxtmp)
       errest(4,:) = errtmp
       erridx(4,:) = idxtmp

       call estimate_error (apar_app, apar_l, kmax, errtmp, idxtmp)
       errest(5,:) = errtmp
       erridx(5,:) = idxtmp
    end if

    if (vary_vnew) then
       vnmult_target = vnmult(2)
       
       if (errest(1,2) > etol + ewindow .or. errest(4,2) > etola + ewindowa) then
          call get_vnewk (vnmult(2), vnmult_target, increase)
       else if (errest(1,2) < etol - ewindow .and. errest(4,2) < etola - ewindowa) then
          call get_vnewk (vnmult(2), vnmult_target, decrease)
       end if

       call init_ediffuse (vnmult_target)
       call init_diffuse_conserve
    end if
       
    if (vary_vnew) then
       vnmult_target = vnmult(1)
       
       if (errest(2,2) > etol + ewindow .or. errest(3,2) > etol + ewindow &
            .or. errest(5,2) > etola + ewindowa) then
          call get_vnewk (vnmult(1), vnmult_target, increase)
       else if (errest(2,2) < etol - ewindow .and. errest(5,2) < etola - ewindowa .and. &
            (errest(3,2) < etol - ewindow .or. .not.trap_flag)) then
          call get_vnewk (vnmult(1), vnmult_target, decrease)
       end if
       
       call init_lorentz (vnmult_target)
       call init_lorentz_conserve
    end if
    
    deallocate (wgt, errtmp, idxtmp)
    if (has_phi) deallocate(phi_app, phi_e, phi_l)
    if (has_apar) deallocate(apar_app, apar_e, apar_l)

  end subroutine get_verr

  !> FIXME : Add documentation    
  subroutine get_vnewk (vnm, vnm_target, incr)
    use collisions, only: vnfac, vnslow

    implicit none

    logical, intent (in) :: incr
    real, intent (in) :: vnm
    real, intent (out) :: vnm_target

    if (incr) then
       vnm_target = vnm * vnfac
    else 
!       vnmult_target =  vnmult * ((1.0-1.0/vnslow)/vnfac + 1.0/vnslow)
       vnm_target =  vnm * vnslow
    end if

  end subroutine get_vnewk

  !> FIXME : Add documentation  
  subroutine estimate_error (app1, app2, kmax_local, errest, erridx, trap)
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid
    use le_grids, only: ng2, jend, il_is_trapped

    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: app1
    complex, dimension (-ntgrid:,:,:,:), intent (in) :: app2
    real, dimension (:,:), intent (in) :: kmax_local
    logical, intent (in), optional :: trap
    real, dimension (:), intent (out) :: errest
    integer, dimension (:), intent (out) :: erridx

    integer :: ik, it, ig, ipt, end_idx
    integer :: igmax, ikmax, itmax, gpcnt
    real :: gdsum, gdmax, gpavg, gnsum, gsmax, gpsum, gptmp

    igmax = 0; ikmax = 0; itmax = 0
    gdsum = 0.0; gdmax = 0.0; gpavg = 0.0; gnsum = 0.0; gsmax = 1.0
    do ik = 1, naky
       do it = 1, ntheta0
          do ig=-ntgrid,ntgrid
             if (.not.present(trap) .or. il_is_trapped(jend(ig))) then
                gpcnt = 0; gpsum = 0.0

                   if (present(trap)) then
                      end_idx = jend(ig)-ng2
                   else
                      end_idx = size(app2,4)
                   end if

                   do ipt=1,end_idx
                      
                      gptmp = kmax_local(it,ik)*cabs(app1(ig,it,ik) - app2(ig,it,ik,ipt))
                      gpsum = gpsum + gptmp
                      gpcnt = gpcnt + 1
                      
                   end do
                   
                   gpavg = gpsum/gpcnt
                   
                   if (gpavg > gdmax) then
                      igmax = ig
                      ikmax = ik
                      itmax = it
                      gdmax = gpavg
                      gsmax = kmax_local(it,ik)*cabs(app1(ig,it,ik))
                   end if
                   
                   gnsum = gnsum + gpavg
                   gdsum = gdsum + kmax_local(it,ik)*cabs(app1(ig,it,ik))

!                end if
             end if
          end do
       end do
    end do
       
    gdmax = gdmax/gsmax
       
    erridx(1) = igmax
    erridx(2) = ikmax
    erridx(3) = itmax
    errest(1) = gdmax

    ! Guard against divide by zero. This can commonly occur on the first time step
    ! in runs with fapar /= 0 due to our initial conditions typically ensuring that
    ! the initial apar is identically zero.
    if (abs(gdsum) < epsilon(0.0)) then
       errest(2) = 0.0
    else
       errest(2) = gnsum/gdsum
    end if

  end subroutine estimate_error

  !> Error estimate based on monitoring amplitudes of legendre polynomial coefficients.
  !>
  !> FIXME: finish documentation
  subroutine get_gtran (geavg, glavg, gtavg, phi, bpar)
    use le_grids, only: legendre_transform, nlambda, ng2, nesub, jend, grid_has_trapped_particles, il_is_trapped
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use dist_fn_arrays, only: gnew, aj0, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_layouts, only: g_lo
    use mp, only: proc0, broadcast
    use array_utils, only: zero_array
    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, intent (out) :: geavg, glavg, gtavg

    real, dimension (:), allocatable :: gne2, gnl2, gnt2
    complex, dimension (:,:,:,:,:), allocatable :: getran, gltran, gttran

    real :: genorm, gemax, genum, gedenom
    real :: glnorm, glmax, glnum, gldenom
    real :: gtnorm, gtmax, gtnum, gtdenom
    integer :: ig, it, ik, is, ie, il, iglo, isgn

    allocate(gne2(0:nesub-1))
    allocate(gnl2(0:ng2-1))

    genorm = 0.0 ; glnorm = 0.0
    gne2  = 0.0 ; gnl2 = 0.0
    gemax = 0.0 ; glmax = 0.0
    genum = 0.0 ; gedenom = 0.0
    glnum = 0.0 ; gldenom = 0.0
    gtnum = 0.0 ; gtdenom = 0.0

    allocate(getran(0:nesub-1,-ntgrid:ntgrid,ntheta0,naky,nspec))
    allocate(gltran(0:ng2-1,-ntgrid:ntgrid,ntheta0,naky,nspec))
 
    call zero_array(getran); call zero_array(gltran)

    if (grid_has_trapped_particles()) then
       allocate(gnt2(0:2*(nlambda-ng2-1)))
       allocate(gttran(0:2*(nlambda-ng2-1),-ntgrid:ntgrid,ntheta0,naky,nspec))      
       gtnorm = 0.0 ; gnt2 = 0.0 ; gtmax = 0.0 ; call zero_array(gttran)
    end if

! transform from g to h
    call g_adjust (gnew, phi, bpar, direction = from_g_gs2)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          do ig=-ntgrid, ntgrid
             g_work(ig,isgn,iglo) = aj0(ig,iglo)*gnew(ig,isgn,iglo)
          end do
       end do
    end do

! perform legendre transform on dist. fn. to obtain coefficients
! used when expanding dist. fn. in legendre polynomials 
    if (allocated(gttran)) then
       call legendre_transform (g_work, getran, gltran, gttran)
    else
       call legendre_transform (g_work, getran, gltran)
    end if

! transform from h back to g
    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    if (proc0) then
       do is=1,nspec
          do ik=1,naky
             do it=1,ntheta0
                do ig=-ntgrid,ntgrid
                   do ie=0,nesub-1
                      gne2(ie) = real(getran(ie,ig,it,ik,is)*conjg(getran(ie,ig,it,ik,is)))
                   end do
                   do il=0,ng2-1
                      gnl2(il) = real(gltran(il,ig,it,ik,is)*conjg(gltran(il,ig,it,ik,is)))
                   end do
                   genorm = maxval(gne2)
                   if (nesub < 3) then
                      gemax = gne2(size(gne2)-1)
                   else
                      gemax = maxval(gne2(nesub-3:nesub-1))
                   end if
                   glnorm = maxval(gnl2)
                   glmax = maxval(gnl2(ng2-3:ng2-1))

                   genum = genum + gemax
                   gedenom = gedenom + genorm
                   glnum = glnum + glmax
                   gldenom = gldenom + glnorm

                   if (grid_has_trapped_particles()) then
                      do il=0,2*(jend(ig)-ng2-1)
                         gnt2(il) = real(gttran(il,ig,it,ik,is)*conjg(gttran(il,ig,it,ik,is)))
                      end do
                      gtnorm = maxval(gnt2(0:2*(jend(ig)-ng2-1)))
                      if (il_is_trapped(jend(ig))) then
                         gtmax = maxval(gnt2(2*(jend(ig)-ng2-1)-2:2*(jend(ig)-ng2-1)))
                      else
                         gtmax = gnt2(0)
                      end if
                      
                      gtnum = gtnum + gtmax
                      gtdenom = gtdenom + gtnorm
                   end if

                end do
             end do
          end do
       end do
       geavg = genum/gedenom
       glavg = glnum/gldenom
       if (grid_has_trapped_particles()) gtavg = gtnum/gtdenom
    end if

    call broadcast (geavg)
    call broadcast (glavg)
    if (grid_has_trapped_particles()) then
       call broadcast (gtavg)
       deallocate (gnt2, gttran)
    end if

    deallocate(gne2, gnl2)    
    deallocate(getran, gltran)

  end subroutine get_gtran

  !> Write \(g(v_\parallel,v_\perp)\) at `ik=it=is=1, ig=0` to text file `<runname>.dist`
  subroutine write_f (last)
    use mp, only: proc0, send, receive
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: g_lo, ik_idx, it_idx, is_idx, il_idx, ie_idx
    use gs2_layouts, only: idx_local, proc_id
    use le_grids, only: al, energy, forbid, negrid, nlambda
    use egrid, only: zeroes, x0
    use theta_grid, only: bmag
    use dist_fn_arrays, only: gnew
    use species, only: nspec

    integer :: iglo, ik, it, is
    integer :: ie, il, ig
    integer, save :: unit
    real :: vpa, vpe
    complex, dimension(2) :: gtmp
    !> If true, also close the file
    logical, intent(in)  :: last 

    real, dimension (:,:), allocatable, save :: xpts
    real, dimension (:), allocatable, save :: ypts

    if (.not. allocated(xpts)) then
       allocate(xpts(negrid,nspec))
       allocate(ypts(nlambda))

    ! should really just get rid of xpts and ypts
       xpts(1:negrid-1,:) = zeroes(1:negrid-1,:)
       xpts(negrid,:) = x0(:)
       ypts = 0.0
       ! change argument of bmag depending on which theta you want to write out
       do il=1,nlambda
          if (1.0-al(il)*bmag(0) .gt. 0.0) ypts(il) = sqrt(1.0-al(il)*bmag(0))
       end do

       if (proc0) then
          call open_output_file (unit, ".dist")
          write(unit, *) negrid*nlambda
       end if
    endif

    do iglo = g_lo%llim_world, g_lo%ulim_world
       ik = ik_idx(g_lo, iglo) ; if (ik /= 1) cycle
       it = it_idx(g_lo, iglo) ; if (it /= 1) cycle
       is = is_idx(g_lo, iglo) ; if (is /= 1) cycle
       ie = ie_idx(g_lo, iglo) 
       ig = 0
       il = il_idx(g_lo, iglo)
       if (idx_local (g_lo, ik, it, il, ie, is)) then
          if (proc0) then 
             gtmp = gnew(ig,:,iglo)
          else
             call send (gnew(ig,:,iglo), 0)
          endif
       else if (proc0) then
          call receive (gtmp, proc_id(g_lo, iglo))
       endif
       if (proc0) then
          if (.not. forbid(ig,il)) then
             vpa = sqrt(energy(ie,is)*max(0.0, 1.0-al(il)*bmag(ig)))
             vpe = sqrt(energy(ie,is)*al(il)*bmag(ig))
             write (unit, "(8(1x,e12.6))") vpa, vpe, energy(ie,is), al(il), &
                  xpts(ie,is), ypts(il), real(gtmp(1)), real(gtmp(2))
          end if
       end if
    end do
    if (proc0) write (unit, *)
    if (last .and. proc0) call close_output_file (unit)
    
  end subroutine write_f

  !> Write distribution function (\(g\) or possibly \(f\)?) in real space to text file
  !> "<runname>.yxdist"
  subroutine write_fyx (phi,bpar,last)
    use mp, only: proc0, send, receive, barrier
    use file_utils, only: open_output_file, close_output_file
    use gs2_layouts, only: il_idx, ig_idx, ik_idx, it_idx, is_idx, isign_idx, ie_idx
    use gs2_layouts, only: idx_local, proc_id, yxf_lo, accelx_lo, g_lo
    use le_grids, only: al, energy=>energy_maxw, forbid, negrid, nlambda
    use theta_grid, only: bmag, ntgrid
    use dist_fn_arrays, only: gnew, g_adjust, g_work, to_g_gs2, from_g_gs2
    use nonlinear_terms, only: accelerated
    use gs2_transforms, only: transform2, init_transforms
    use array_utils, only: zero_array, copy
#ifdef SHMEM
    use shm_mpi3, only: shm_alloc, shm_free
#endif
    !> Electrostatic potential and parallel magnetic field
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    !> If true, also close the file
    logical, intent (in) :: last

    real, dimension (:,:), allocatable :: grs, gzf
#ifndef SHMEM
    real, dimension (:,:,:), allocatable :: agrs, agzf
#else
    real, dimension (:,:,:), pointer, contiguous :: agrs => null(), agzf => null()
    complex, dimension(:,:,:), pointer, contiguous :: g0_ptr => null()
#endif
    real, dimension (:), allocatable :: agp0, agp0zf
    real :: gp0, gp0zf
    integer :: ig, it, ik, il, ie, is, iyxlo, isign, ia, iaclo, iglo, acc
    integer, save :: unit

    if (accelerated) then
#ifndef SHMEM
       allocate (agrs(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))          
       allocate (agzf(2*ntgrid+1, 2, accelx_lo%llim_proc:accelx_lo%ulim_alloc))          
#else
       call shm_alloc(agrs, (/ 1, 2*ntgrid+1, 1, 2, &
            accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
       call shm_alloc(agzf, (/ 1, 2*ntgrid+1, 1, 2, &
            accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
        call shm_alloc(g0_ptr, (/ 1, 2*ntgrid+1, 1, 2, &
            accelx_lo%llim_proc, accelx_lo%ulim_alloc /))
#endif
       allocate (agp0(2), agp0zf(2))
       call zero_array(agrs); call zero_array(agzf); agp0 = 0.0; agp0zf = 0.0
    else
#ifdef SHMEM
       g0_ptr => g_work
#endif
       allocate (grs(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
       allocate (gzf(yxf_lo%ny, yxf_lo%llim_proc:yxf_lo%ulim_alloc))
       call zero_array(grs); call zero_array(gzf); gp0 = 0.0; gp0zf = 0.0
    end if

    if (.not. fyxinit) then
       if (proc0) then
          acc = 0
          call open_output_file (unit,".yxdist")
          if (accelerated) acc = 1
          write(unit,*) 2*negrid*nlambda, bmag(0), acc
       end if

       fyxinit = .true.
    end if

    call g_adjust (gnew, phi, bpar, direction = from_g_gs2)

#ifndef SHMEM
    call copy(gnew, g_work)
#else
    g0_ptr = gnew
#endif

#ifndef SHMEM
    if (accelerated) then
       call transform2 (g_work, agrs, ia)
    else
       call transform2 (g_work, grs)
    end if

    call zero_array(g_work)
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       if (ik == 1) g_work(:,:,iglo) = gnew(:,:,iglo)
    end do

    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    if (accelerated) then
       call transform2 (g_work, agzf, ia)
    else
       call transform2 (g_work, gzf)
    end if

#else

    if (accelerated) then
       call transform2 (g0_ptr, agrs, ia)
    else
       call transform2 (g0_ptr, grs)
    end if

    g0_ptr = 0.0
    do iglo=g_lo%llim_proc, g_lo%ulim_proc
       ik = ik_idx(g_lo,iglo)
       if (ik == 1) g0_ptr(:,:,iglo) = gnew(:,:,iglo)
    end do

    call g_adjust (gnew, phi, bpar, direction = to_g_gs2)

    if (accelerated) then
       call transform2 (g0_ptr, agzf, ia)
    else
       call transform2 (g0_ptr, gzf)
    end if

#endif

    if (accelerated) then
       do iaclo=accelx_lo%llim_world, accelx_lo%ulim_world
          it = it_idx(accelx_lo, iaclo)
          ik = ik_idx(accelx_lo, iaclo)
          if (it == 1 .and. ik == 1) then
             il = il_idx(accelx_lo, iaclo)
             ig = 0
             if (.not. forbid(ig,il)) then
                ie = ie_idx(accelx_lo, iaclo)
                is = is_idx(accelx_lo, iaclo)

                if (proc0) then
                   if (.not. idx_local(accelx_lo, ik, it, il, ie, is)) then
                      call receive (agp0, proc_id(accelx_lo, iaclo))
                      call receive (agp0zf, proc_id(accelx_lo, iaclo))
                   else
                      agp0 = agrs(ig+ntgrid+1, :, iaclo)
                      agp0zf = agzf(ig+ntgrid+1, :, iaclo)
                   end if
                else if (idx_local(accelx_lo, ik, it, il, ie, is)) then
                   call send (agrs(ig+ntgrid+1, :, iaclo), 0)
                   call send (agzf(ig+ntgrid+1, :, iaclo), 0)
                end if
                
                if (proc0) then
                   write (unit, "(6(1x,e13.6))") energy(ie), al(il), &
                        agp0(1), agp0(2), agp0zf(1), agp0zf(2)
                end if
             end if
          end if
          call barrier
       end do
       deallocate(agp0, agp0zf)
#ifndef SHMEM
       deallocate(agrs, agzf)
#else
       call shm_free(agrs)
       call shm_free(agzf)
       call shm_free(g0_ptr)
#endif
    else
       do iyxlo=yxf_lo%llim_world, yxf_lo%ulim_world
          
          ig = ig_idx(yxf_lo, iyxlo)
          it = it_idx(yxf_lo, iyxlo)
          if (ig == 0 .and. it == 1) then
             il = il_idx(yxf_lo, iyxlo)
             if (.not. forbid(ig,il)) then
                ik = 1
                ie = ie_idx(yxf_lo, iyxlo)
                is = is_idx(yxf_lo, iyxlo)
                isign = isign_idx(yxf_lo, iyxlo)
                
                if (proc0) then
                   if (.not. idx_local(yxf_lo, ig, isign, it, il, ie, is)) then
                      call receive (gp0, proc_id(yxf_lo, iyxlo))
                      call receive (gp0zf, proc_id(yxf_lo, iyxlo))
                   else
                      gp0 = grs(ik, iyxlo)
                      gp0zf = gzf(ik, iyxlo)
                   end if
                else if (idx_local(yxf_lo, ig, isign, it, il, ie, is)) then
                   call send (grs(ik, iyxlo), 0)
                   call send (gzf(ik, iyxlo), 0)                   
                end if
                
                if (proc0) then
                   write (unit, "(4(1x,e13.6),i8)") energy(ie), al(il), &
                        gp0, gp0zf, isign
                end if
             end if
          end if
          call barrier
       end do
       deallocate (grs, gzf)
    end if

    if (proc0) flush (unit)

    if (proc0) then
       write(unit,*)
       if (last) call close_output_file (unit)
    end if

  end subroutine write_fyx

  !> Calculate the collision error and write to text file `<runname>.cres`
  subroutine collision_error (phi, bpar, last)
    use mp, only: proc0, send, receive, barrier
    use le_grids, only: ng2, jend, nlambda, lambda_map
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, g_adjust, g_work, to_g_gs2, from_g_gs2
    use gs2_layouts, only: lz_lo, ig_idx, idx_local, proc_id
    use gs2_layouts, only: ik_idx, ie_idx, is_idx, it_idx, il_idx
    use collisions, only: dtot, fdf, fdb
    use redistribute, only: gather, scatter
    use file_utils, only: open_output_file, close_output_file
    use gs2_time, only: user_time
    use array_utils, only: zero_array
    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    !> If true, close the file
    logical, intent (in) :: last

    integer :: je, te, ig, il, ip, ilz, ie, is, ik, it
    integer :: igmax, ikmax, itmax, iemax, ilmax, ismax
    integer, save :: unit
    complex, dimension (:), allocatable :: ltmp, ftmp
    complex, dimension (:,:), allocatable :: lcoll, fdcoll, glze
!    logical :: first = .true.
    real :: etmp, emax, etot, eavg, edenom, ltmax
    real :: time

    allocate (ltmp(2*nlambda), ftmp(2*nlambda))
    allocate (lcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (fdcoll(2*nlambda,lz_lo%llim_proc:lz_lo%ulim_alloc))
    allocate (glze(max(2*nlambda,2*ng2+1),lz_lo%llim_proc:lz_lo%ulim_alloc))

    call zero_array(lcoll); call zero_array(fdcoll); call zero_array(glze)
    ltmp = 0.0; ftmp = 0.0
    etmp = 0.0; emax = 0.0; etot = 0.0; eavg = 0.0; edenom = 1.0; ltmax = 1.0

!    if (first .and. proc0) then
    if (.not.cerrinit .and. proc0) then
       call open_output_file (unit,".cres")
!       first = .false.
       cerrinit = .true.
    end if

! convert gnew from g to h
    call g_adjust(gnew,phi,bpar, direction = from_g_gs2)

    g_work = gnew

! convert gnew from h back to g
    call g_adjust(gnew,phi,bpar, direction = to_g_gs2)

! map from g_work(ig,isgn,iglo) to glze(il,ilz)
    call gather (lambda_map, g_work, glze)

! loop over ig, isign, ik, it, ie, is
    do ilz = lz_lo%llim_proc, lz_lo%ulim_proc
       ig = ig_idx(lz_lo,ilz)
       
       if (jend(ig) == 0) then      ! no trapped particles
! je = number of + vpa grid pts
! te = number of + and - vpa grid pts
          je = ng2
          te = 2*ng2
          
! find d/d(xi) ((1+xi**2)( d g(xi)/ d(xi) )) at each xi
! using lagrange (lcoll) and finite difference (fdcoll)
          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,ng2
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=ng2+1, 2*ng2-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2*ng2-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = 2*ng2-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = 2*ng2
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! deal with xi from 1-eps -> eps
          do il=2,ng2
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
! deal with xi from -eps -> -1+eps
          do il=ng2+1, 2*ng2-1
             fdcoll(il,ilz) = fdb(ig,2*ng2-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,2*ng2-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do

          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(2*ng2,ilz) = -fdf(ig,1)*(glze(2*ng2,ilz) - glze(2*ng2-1,ilz))

       else       ! trapped particle runs          
          je = jend(ig)
          te = 2*je - 1

          il = 1
          do ip = il, il+2
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip)*glze(ip,ilz)
          end do

          il = 2
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+2)*glze(ip,ilz)
          end do

          do il=3,je
             do ip=il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,il,ip-il+3)*glze(ip,ilz)
             end do
          end do

          do il=je+1, te-2
             do ip = il-2,il+2
                lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,te-il+1,il-ip+3)*glze(ip,ilz)
             end do
          end do

          il = te-1
          do ip = il-1, il+1
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,2,il-ip+2)*glze(ip,ilz)
          end do

          il = te
          do ip = il-2, il
             lcoll(il,ilz) = lcoll(il,ilz) + dtot(ig,1,il-ip+1)*glze(ip,ilz)
          end do

! is il=je handled correctly here?
          do il=2,je
             fdcoll(il,ilz) = fdf(ig,il)*(glze(il+1,ilz) - glze(il,ilz)) - fdb(ig,il)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          do il=je+1,te-1
             fdcoll(il,ilz) = fdb(ig,te-il+1)*(glze(il+1,ilz) - glze(il,ilz)) - fdf(ig,te-il+1)*(glze(il,ilz) - glze(il-1,ilz))
          end do
          
          fdcoll(1,ilz) = fdf(ig,1)*(glze(2,ilz) - glze(1,ilz))
          fdcoll(te,ilz) = -fdf(ig,1)*(glze(te,ilz) - glze(te-1,ilz))
          
       end if
    end do

    time = user_time

    do ilz=lz_lo%llim_world, lz_lo%ulim_world
       ig = ig_idx(lz_lo, ilz)
       ik = ik_idx(lz_lo, ilz)
       it = it_idx(lz_lo, ilz)
       ie = ie_idx(lz_lo, ilz)
       is = is_idx(lz_lo, ilz)
       je = jend(ig)

       if (je == 0) then
          te = 2*ng2
       else
          te = 2*je-1
       end if

       if (idx_local (lz_lo, ilz)) then
          if (proc0) then 
             ltmp = lcoll(:,ilz)
             ftmp = fdcoll(:,ilz)
          else
             call send (lcoll(:,ilz), 0)
             call send (fdcoll(:,ilz), 0)
          endif
       else if (proc0) then
          call receive (ltmp, proc_id(lz_lo, ilz))
          call receive (ftmp, proc_id(lz_lo, ilz))
       endif
       call barrier

       do il=1,te
          etmp = cabs(ltmp(il) - ftmp(il))
          
          if (etmp > emax) then
             emax = etmp
             ltmax = cabs(ltmp(il))
             ikmax = ik
             itmax = it
             iemax = ie
             ismax = is
             ilmax = il
             igmax = ig
          end if
          
          etot = etot + etmp
          edenom = edenom + cabs(ltmp(il))
       end do
    end do

    eavg = etot/edenom
    emax = emax/ltmax

    if (proc0) then
       write(unit,"((1x,e13.6),6(i8),2(1x,e13.6))") time, &
            igmax, ikmax, itmax, iemax, ilmax, ismax, emax, eavg
       if (last) then
          call close_output_file (unit)
       end if
    end if

    deallocate (lcoll, fdcoll, glze, ltmp, ftmp)
    
  end subroutine collision_error

  !> FIXME : Add documentation
  subroutine boundary(linked)
    implicit none
    logical, intent(out) :: linked
    call init_dist_fn
    linked = boundary_option_switch == boundary_option_linked
  end subroutine boundary

  !> FIXME : Add documentation
  subroutine find_leftmost_link (iglo, iglo_left, ipleft)
    use gs2_layouts, only: it_idx,ik_idx,g_lo,il_idx,ie_idx,is_idx,idx,proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_left, ipleft
    integer :: iglo_star
    integer :: it_cur,ik,it,il,ie,is
    iglo_star = iglo
    it_cur=it_idx(g_lo,iglo)
    it=it_cur
    ik=ik_idx(g_lo,iglo)

    !Now get the leftmost it
    it_cur=get_leftmost_it(it,ik)

    !If we're at the same it then don't need to do much
    if(it.eq.it_cur)then
       iglo_left=iglo
       ipleft=proc_id(g_lo,iglo)
       return
    endif

    !If not then we need to calculate iglo_left and ipleft
    il=il_idx(g_lo,iglo)
    ie=ie_idx(g_lo,iglo)
    is=is_idx(g_lo,iglo)
    iglo_left=idx(g_lo,ik,it_cur,il,ie,is)
    ipleft=proc_id(g_lo,iglo_left)

  end subroutine find_leftmost_link

  !> FIXME : Add documentation  
  subroutine find_rightmost_link (iglo, iglo_right, ipright)
    use gs2_layouts, only: it_idx,ik_idx,g_lo,il_idx,ie_idx,is_idx,idx,proc_id
    implicit none
    integer, intent (in) :: iglo
    integer, intent (out) :: iglo_right, ipright
    integer :: iglo_star
    integer :: it_cur,ik,it,il,ie,is
    iglo_star = iglo
    it_cur=it_idx(g_lo,iglo)
    ik=ik_idx(g_lo,iglo)
    it=it_cur

    !Now get the rightmost it
    it_cur=get_rightmost_it(it,ik)

    !If we're at the same it then don't need to do much
    if(it.eq.it_cur)then
       iglo_right=iglo
       ipright=proc_id(g_lo,iglo)
       return
    endif

    !If not then we need to calculate iglo_left and ipleft
    il=il_idx(g_lo,iglo)
    ie=ie_idx(g_lo,iglo)
    is=is_idx(g_lo,iglo)
    iglo_right=idx(g_lo,ik,it_cur,il,ie,is)
    ipright=proc_id(g_lo,iglo_right)

  end subroutine find_rightmost_link


  !> Get a theta-centered and time-centered estimate of the time derivative 
  !! of a field.
  !! 
  !! tunits(ky) == 1. unless the "omega_*" units are chosen.
  !! omega_* units normalize time by an additional factor of ky.
  subroutine dot (a, anew, adot, fac)
    use gs2_time, only: code_dt, tunits
    use kt_grids, only: naky, ntheta0
    use theta_grid, only: ntgrid

    implicit none
    complex, intent (in), dimension (-ntgrid:,:,:) :: a, anew
    complex, intent (out), dimension (-ntgrid:,:,:) :: adot
    real, intent (in) :: fac
    real :: dtinv
    integer :: ig, it, ik

    do ik=1,naky
       dtinv = 1./(code_dt*tunits(ik))
       do it=1,ntheta0
          do ig=-ntgrid,ntgrid-1
             adot(ig,it,ik) = 0.5*fac*(anew(ig+1,it,ik)+anew(ig,it,ik) - &
                  (a(ig+1,it,ik)+a(ig,it,ik)))*dtinv
          end do
       end do
    end do
    
  end subroutine dot

  !> FIXME : Add documentation
  complex function fdot (fl, fr, fnewl, fnewr, dtinv)
    complex, intent (in) :: fl, fr, fnewl, fnewr
    real, intent (in) :: dtinv
    fdot = 0.5*(fnewl+fnewr-fl-fr)*dtinv
  end function fdot

  !> Construct time and space-centered quantities
  !! (should use actual bakdif and fexpr values?)
  complex function favg (fl, fr, fnewl, fnewr)
    complex, intent (in) :: fl, fr, fnewl, fnewr
    favg = 0.25*(fnewl+fnewr+fl+fr)
  end function favg

  !> FIXME : Add documentation
  complex function fdot_t (f,fnew, dtinv)
    complex, intent (in) :: f, fnew
    real, intent (in) :: dtinv
    fdot_t = (fnew-f)*dtinv
  end function fdot_t

  !> FIXME : Add documentation  
  complex function favg_x (fl, fr)
    complex, intent (in) :: fl, fr
    favg_x = 0.5*(fl+fr)
  end function favg_x

  !> FIXME : Add documentation
  !>
  !> The returned variables are used if write_full_moments_notgc=T
  !> or ginit_options='recon3'.
  subroutine init_mom_coeff(mom_coeff, mom_coeff_npara, mom_coeff_nperp, mom_coeff_tpara, &
       mom_coeff_tperp, mom_shift_para, mom_shift_perp)
    use gs2_layouts, only: g_lo
    use species, only: nspec, spec
    use kt_grids, only: nakx => ntheta0, naky, kperp2
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: vpa, vperp2, aj0, aj1
    use le_grids, only: integrate_moment
    use array_utils, only: zero_array
    implicit none
    real, dimension(:, :, :, :), allocatable, intent(in out) :: mom_coeff
    real, dimension(:, :, :), allocatable, intent(in out) :: mom_coeff_npara, mom_coeff_nperp
    real, dimension(:, :, :), allocatable, intent(in out) :: mom_coeff_tpara, mom_coeff_tperp
    real, dimension(:, :, :), allocatable, intent(in out) :: mom_shift_para, mom_shift_perp
    integer :: i, isgn, iglo, it, ik, is
    real :: bsq
    complex, dimension(:, :, :, :), allocatable :: coeff0
    complex, dimension(:, :, :), allocatable :: gtmp
    real, dimension(:), allocatable :: wgt
    logical, parameter :: analytical = .true.
    integer, parameter :: ncnt_mom_coeff=8

    if (allocated(mom_coeff)) return

    if(.not.allocated(mom_coeff)) &
         & allocate(mom_coeff(nakx,naky,nspec,ncnt_mom_coeff))
    if(.not.allocated(mom_coeff_npara)) &
         & allocate(mom_coeff_npara(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_nperp)) &
         & allocate(mom_coeff_nperp(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_tpara)) &
         & allocate(mom_coeff_tpara(nakx,naky,nspec))
    if(.not.allocated(mom_coeff_tperp)) &
         & allocate(mom_coeff_tperp(nakx,naky,nspec))
    if(.not.allocated(mom_shift_para)) &
         & allocate(mom_shift_para(nakx,naky,nspec))
    if(.not.allocated(mom_shift_perp)) &
         & allocate(mom_shift_perp(nakx,naky,nspec))

    call zero_array(mom_coeff)
    call zero_array(mom_coeff_npara) ; call zero_array(mom_coeff_nperp)
    call zero_array(mom_coeff_tpara) ; call zero_array(mom_coeff_tperp)
    call zero_array(mom_shift_para) ; call zero_array(mom_shift_perp)

    allocate(wgt(-ntgrid:ntgrid))
    allocate(coeff0(-ntgrid:ntgrid,nakx,naky,nspec))
    allocate(gtmp(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    wgt = 0.
    call zero_array(coeff0)
    call zero_array(gtmp)

    if (analytical) then
       do it=1,nakx
          do ik=1,naky
             do is=1,nspec
                bsq=.25*spec(is)%smz**2*kperp2(0,it,ik)
                mom_coeff(it,ik,is,1) = exp(-bsq)
                mom_coeff(it,ik,is,2) = exp(-bsq) *.5
                mom_coeff(it,ik,is,3) = exp(-bsq) *(1.-bsq)
                mom_coeff(it,ik,is,4) = exp(-bsq) *.75
                mom_coeff(it,ik,is,5) = exp(-bsq) *(1.-bsq)*.5
                mom_coeff(it,ik,is,6) = exp(-bsq) *.5
                mom_coeff(it,ik,is,7) = exp(-bsq) *.25
                mom_coeff(it,ik,is,8) = exp(-bsq) *(1.-.5*bsq)
             end do
          end do
       end do
    else
       do i = 1, ncnt_mom_coeff
          do iglo = g_lo%llim_proc, g_lo%ulim_proc
             do isgn = 1,2
                if(i==1) wgt(:)=aj0(:,iglo)
                if(i==2) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**2
                if(i==3) wgt(:)=aj0(:,iglo)*vperp2(:,iglo)
                if(i==4) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**4
                if(i==5) wgt(:)=aj0(:,iglo)*vpa(:,isgn,iglo)**2*vperp2(:,iglo)
                if(i==6) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)
                if(i==7) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)*vpa(:,isgn,iglo)**2
                if(i==8) wgt(:)=vperp2(:,iglo)*aj1(:,iglo)*vperp2(:,iglo)
                gtmp(-ntgrid:ntgrid,isgn,iglo) = wgt(-ntgrid:ntgrid)*cmplx(1.,0.)
             end do
          end do
          call integrate_moment(gtmp,coeff0,.true.,full_arr=.true.)
          where(real(coeff0(0,1:nakx,1:naky,1:nspec)) == 0.)
             mom_coeff(1:nakx,1:naky,1:nspec,i)=1.
          elsewhere
             mom_coeff(1:nakx,1:naky,1:nspec,i)= &
                  & coeff0(0,1:nakx,1:naky,1:nspec)
          end where
       end do
    endif

    !<DD>Currently below could include divide by zero if analytical=.true.
    mom_shift_para=mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1)
    mom_shift_perp=mom_coeff(:,:,:,3)/mom_coeff(:,:,:,1)

    mom_coeff_npara=2.*mom_coeff(:,:,:,2)/mom_coeff(:,:,:,1)
    mom_coeff_nperp=2.*mom_coeff(:,:,:,6)/mom_coeff(:,:,:,1)

    mom_coeff_tperp= &
         & (mom_coeff(:,:,:,5)-mom_shift_perp*mom_coeff(:,:,:,2)) / &
         & (mom_coeff(:,:,:,8)-mom_shift_perp*mom_coeff(:,:,:,6))
    mom_coeff_tpara= &
         & (mom_coeff(:,:,:,7)-mom_shift_para*mom_coeff(:,:,:,6)) / &
         & (mom_coeff(:,:,:,4)-mom_shift_para*mom_coeff(:,:,:,2))

    deallocate(gtmp,coeff0,wgt)
  end subroutine init_mom_coeff
  
  !> Calculate the flux surface average term for the adiabatic response.
  !>
  !> If using adiabatic electrons and the option
  !> dist_fn_knobs
  !>   adiabatic_option = "iphi00=2"
  !> /
  !> the field solve should subtract the flux surface average
  !> from the perturbed electron distribution function,
  !> i.e. include the term <phi> in
  !>
  !>   f1e = q( phi - <phi>) f0e / Te
  !>
  !> This function calculates <phi> and returns it as fl_avg.
  !>
  !> Joseph Parker, STFC
  !>
  !> @note We might want to rename this routine to reflect
  !> that it only calculates a very specific flux surface average
  !> and not a general average of any field shaped array (e.g. as
  !> gamtot3 is only allocated/set in a special case and the gamtot
  !> factor is not required in the general flux surface average).
  !>
  !> @note antot/gamtot is the expression for phi when beta = 0
  !> and this reflects what is averaged in this routine. As such we
  !> appear to be assuming beta = 0 here. This is probably fine as
  !> we should only use the result for adiabatic electrons (where
  !> one would expect to be ignoring EM perturbations) but perhaps
  !> there should be a consistency check? This is something we should
  !> still flag and could rectify with a little effort.
  subroutine calculate_flux_surface_average (fl_avg,antot)
    use kt_grids, only: naky, ntheta0, aky, kwork_filter
    use run_parameters, only: tite
    use theta_grid, only: ntgrid, theta, jacob
    use integration, only: trapezoidal_integration
    implicit none
    complex, dimension (ntheta0, naky), intent (out) :: fl_avg
    complex, dimension (-ntgrid:ntgrid,ntheta0,naky), intent (in) :: antot
    integer :: ik, it

    fl_avg = 0.

    if (.not. allocated(awgt)) then
       allocate (awgt(ntheta0, naky))
       awgt = 0.
       do ik = 1, naky
          if (aky(ik) > epsilon(0.0)) cycle
          do it = 1, ntheta0
             ! We might not want this cycle given that this is a one-off setup of awgt
             ! and kwork_filter may change during the run (e.g. between initialisation
             ! and time advance). I suspect kwork_filter will be false the first
             ! time we call this routine, but there is no promise.
             if(kwork_filter(it,ik)) cycle
             awgt(it,ik) = 1.0/trapezoidal_integration(theta, jacob*gamtot3(:,it,ik))
          end do
       end do
    endif

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(ik, it) &
    !$OMP SHARED(naky, ntheta0, kwork_filter, fl_avg, tite, antot, &
    !$OMP theta, jacob, gamtot, awgt) &
    !$OMP SCHEDULE(static)
    do ik = 1, naky
       ! Note awg == 0 for aky > 0 so we might be able to cycle/limit the loop.
       do it = 1, ntheta0
          if(kwork_filter(it,ik)) cycle
          fl_avg(it,ik) = tite*trapezoidal_integration(theta, jacob*antot(:,it,ik)/gamtot(:,it,ik))*awgt(it,ik)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine calculate_flux_surface_average

  !> Set the module level config types
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_dist_fn_config(dist_fn_config_in, source_config_in, dist_fn_species_config_in)
    use mp, only: mp_abort
    type(dist_fn_config_type), intent(in), optional :: dist_fn_config_in
    type(source_config_type), intent(in), optional :: source_config_in
    type(dist_fn_species_config_type), intent(in), dimension(:), allocatable, optional :: dist_fn_species_config_in
 
    if (initialized) then
       call mp_abort("Trying to set dist_fn_config when already initialized.", to_screen = .true.)
    end if
    if (present(dist_fn_config_in)) then
       dist_fn_config = dist_fn_config_in
    end if
    if (present(source_config_in)) then
       source_config = source_config_in
    end if
    if (present(dist_fn_species_config_in)) then
       dist_fn_species_config = dist_fn_species_config_in
    end if    
  end subroutine set_dist_fn_config   

  !> Get the module level config instance
  function get_dist_fn_config()
    type(dist_fn_config_type) :: get_dist_fn_config
    get_dist_fn_config = dist_fn_config
  end function get_dist_fn_config

  !> Get the module level config instance
  function get_dist_fn_species_config()
    type(dist_fn_species_config_type), allocatable :: get_dist_fn_species_config(:)
    if (allocated(dist_fn_species_config)) then
      get_dist_fn_species_config = dist_fn_species_config
    else
      allocate(get_dist_fn_species_config(0))
    end if
  end function get_dist_fn_species_config

  !> Get the module level config instance
  function get_source_config()
    type(source_config_type) :: get_source_config
    get_source_config = source_config
  end function get_source_config

  !---------------------------------------
  ! Following is for the dist_fn config_type
  !---------------------------------------
  
  !> Reads in the dist_fn_knobs namelist and populates the member variables
  subroutine read_dist_fn_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(dist_fn_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: boundary_option
    character(len = 30) :: adiabatic_option
    integer :: g_exb_start_timestep
    logical :: def_parity, esv, even, gf_lo_integrate, hyper_in_initialisation, lf_decompose, lf_default, mult_imp, nonad_zero, opt_source, start_from_previous_solution, zero_forbid
    real :: afilter, apfac, btor_slab, driftknob, g_exb, g_exb_error_limit, g_exb_start_time, g_exbfac, gridfac, mach, omprimfac, poisfac, tpdriftknob, vparknob, wfb

    namelist /dist_fn_knobs/ adiabatic_option, afilter, apfac, boundary_option, btor_slab, def_parity, driftknob, esv, even, g_exb, g_exb_error_limit, g_exb_start_time, g_exb_start_timestep, g_exbfac, &
         gf_lo_integrate, gridfac, hyper_in_initialisation, lf_decompose, lf_default, mach, mult_imp, nonad_zero, omprimfac, opt_source, poisfac, &
         start_from_previous_solution, tpdriftknob, vparknob, wfb, zero_forbid

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    adiabatic_option = self%adiabatic_option
    afilter = self%afilter
    apfac = self%apfac
    boundary_option = self%boundary_option
    btor_slab = self%btor_slab
    def_parity = self%def_parity
    driftknob = self%driftknob
    esv = self%esv
    even = self%even
    g_exb = self%g_exb
    g_exb_error_limit = self%g_exb_error_limit
    g_exb_start_time = self%g_exb_start_time
    g_exb_start_timestep = self%g_exb_start_timestep
    g_exbfac = self%g_exbfac
    gf_lo_integrate = self%gf_lo_integrate
    gridfac = self%gridfac
    hyper_in_initialisation = self%hyper_in_initialisation
    lf_decompose = self%lf_decompose
    lf_default = self%lf_default
    mach = self%mach
    mult_imp = self%mult_imp
    nonad_zero = self%nonad_zero
    omprimfac = self%omprimfac
    opt_source = self%opt_source
    poisfac = self%poisfac
    start_from_previous_solution = self%start_from_previous_solution
    tpdriftknob = self%tpdriftknob
    vparknob = self%vparknob
    wfb = self%wfb
    zero_forbid = self%zero_forbid

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = dist_fn_knobs)

    ! Now copy from local variables into type members
    self%adiabatic_option = adiabatic_option
    self%afilter = afilter
    self%apfac = apfac
    self%boundary_option = boundary_option
    self%btor_slab = btor_slab
    self%def_parity = def_parity
    self%driftknob = driftknob
    self%esv = esv
    self%even = even
    self%g_exb = g_exb
    self%g_exb_error_limit = g_exb_error_limit
    self%g_exb_start_time = g_exb_start_time
    self%g_exb_start_timestep = g_exb_start_timestep
    self%g_exbfac = g_exbfac
    self%gf_lo_integrate = gf_lo_integrate
    self%gridfac = gridfac
    self%hyper_in_initialisation = hyper_in_initialisation
    self%lf_decompose = lf_decompose
    self%lf_default = lf_default
    self%mach = mach
    self%mult_imp = mult_imp
    self%nonad_zero = nonad_zero
    self%omprimfac = omprimfac
    self%opt_source = opt_source
    self%poisfac = poisfac
    self%start_from_previous_solution = start_from_previous_solution
    self%tpdriftknob = tpdriftknob
    self%vparknob = vparknob
    self%wfb = wfb
    self%zero_forbid = zero_forbid

    self%exist = exist
  end subroutine read_dist_fn_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_dist_fn_config(self, unit)
    implicit none
    class(dist_fn_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("adiabatic_option", self%adiabatic_option, unit_internal)
    call self%write_key_val("afilter", self%afilter, unit_internal)
    call self%write_key_val("apfac", self%apfac, unit_internal)
    call self%write_key_val("boundary_option", self%boundary_option, unit_internal)
    call self%write_key_val("btor_slab", self%btor_slab, unit_internal)
    call self%write_key_val("def_parity", self%def_parity, unit_internal)
    call self%write_key_val("driftknob", self%driftknob, unit_internal)
    call self%write_key_val("esv", self%esv, unit_internal)
    call self%write_key_val("even", self%even, unit_internal)
    call self%write_key_val("g_exb", self%g_exb, unit_internal)
    call self%write_key_val("g_exb_error_limit", self%g_exb_error_limit, unit_internal)
    call self%write_key_val("g_exb_start_time", self%g_exb_start_time, unit_internal)
    call self%write_key_val("g_exb_start_timestep", self%g_exb_start_timestep, unit_internal)
    call self%write_key_val("g_exbfac", self%g_exbfac, unit_internal)
    call self%write_key_val("gf_lo_integrate", self%gf_lo_integrate, unit_internal)
    call self%write_key_val("gridfac", self%gridfac, unit_internal)
    call self%write_key_val("hyper_in_initialisation", self%hyper_in_initialisation, unit_internal)
    call self%write_key_val("lf_decompose", self%lf_decompose, unit_internal)
    call self%write_key_val("lf_default", self%lf_default, unit_internal)
    call self%write_key_val("mach", self%mach, unit_internal)
    call self%write_key_val("mult_imp", self%mult_imp, unit_internal)
    call self%write_key_val("nonad_zero", self%nonad_zero, unit_internal)
    call self%write_key_val("omprimfac", self%omprimfac, unit_internal)
    call self%write_key_val("opt_source", self%opt_source, unit_internal)
    call self%write_key_val("poisfac", self%poisfac, unit_internal)
    call self%write_key_val("start_from_previous_solution", self%start_from_previous_solution, unit_internal)
    call self%write_key_val("tpdriftknob", self%tpdriftknob, unit_internal)
    call self%write_key_val("vparknob", self%vparknob, unit_internal)
    call self%write_key_val("wfb", self%wfb, unit_internal)
    call self%write_key_val("zero_forbid", self%zero_forbid, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_dist_fn_config

  !> Resets the config object to the initial empty state
  subroutine reset_dist_fn_config(self)
    class(dist_fn_config_type), intent(in out) :: self
    type(dist_fn_config_type) :: empty
    select type (self)
    type is (dist_fn_config_type)
       self = empty
    end select
  end subroutine reset_dist_fn_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_dist_fn_config(self)
    use mp, only: broadcast
    implicit none
    class(dist_fn_config_type), intent(in out) :: self
    call broadcast(self%adiabatic_option)
    call broadcast(self%afilter)
    call broadcast(self%apfac)
    call broadcast(self%boundary_option)
    call broadcast(self%btor_slab)
    call broadcast(self%def_parity)
    call broadcast(self%driftknob)
    call broadcast(self%esv)
    call broadcast(self%even)
    call broadcast(self%g_exb)
    call broadcast(self%g_exb_error_limit)
    call broadcast(self%g_exb_start_time)
    call broadcast(self%g_exb_start_timestep)
    call broadcast(self%g_exbfac)
    call broadcast(self%gf_lo_integrate)
    call broadcast(self%gridfac)
    call broadcast(self%hyper_in_initialisation)
    call broadcast(self%lf_decompose)
    call broadcast(self%lf_default)
    call broadcast(self%mach)
    call broadcast(self%mult_imp)
    call broadcast(self%nonad_zero)
    call broadcast(self%omprimfac)
    call broadcast(self%opt_source)
    call broadcast(self%poisfac)
    call broadcast(self%start_from_previous_solution)
    call broadcast(self%tpdriftknob)
    call broadcast(self%vparknob)
    call broadcast(self%wfb)
    call broadcast(self%zero_forbid)

    call broadcast(self%exist)
  end subroutine broadcast_dist_fn_config
  
  !> Gets the default name for this namelist
  function get_default_name_dist_fn_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_dist_fn_config
    get_default_name_dist_fn_config = "dist_fn_knobs"
  end function get_default_name_dist_fn_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_dist_fn_config()
    implicit none
    logical :: get_default_requires_index_dist_fn_config
    get_default_requires_index_dist_fn_config = .false.
  end function get_default_requires_index_dist_fn_config
  
  !---------------------------------------
  ! Following is for the source config_type
  !---------------------------------------
  
  !> Reads in the source_knobs namelist and populates the member variables
  subroutine read_source_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(source_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    character(len = 20) :: source_option
    real :: gamma0, omega0, phi_ext, source0, t0

    namelist /source_knobs/ gamma0, omega0, phi_ext, source0, source_option, t0

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    gamma0 = self%gamma0
    omega0 = self%omega0
    phi_ext = self%phi_ext
    source0 = self%source0
    source_option = self%source_option
    t0 = self%t0

    ! Now read in the main namelist
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = source_knobs)

    ! Now copy from local variables into type members
    self%gamma0 = gamma0
    self%omega0 = omega0
    self%phi_ext = phi_ext
    self%source0 = source0
    self%source_option = source_option
    self%t0 = t0

    self%exist = exist
  end subroutine read_source_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_source_config(self, unit)
    implicit none
    class(source_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("gamma0", self%gamma0, unit_internal)
    call self%write_key_val("omega0", self%omega0, unit_internal)
    call self%write_key_val("phi_ext", self%phi_ext, unit_internal)
    call self%write_key_val("source0", self%source0, unit_internal)
    call self%write_key_val("source_option", self%source_option, unit_internal)
    call self%write_key_val("t0", self%t0, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_source_config

  !> Resets the config object to the initial empty state
  subroutine reset_source_config(self)
    class(source_config_type), intent(in out) :: self
    type(source_config_type) :: empty
    select type (self)
    type is (source_config_type)
       self = empty
    end select
  end subroutine reset_source_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_source_config(self)
    use mp, only: broadcast
    implicit none
    class(source_config_type), intent(in out) :: self
    call broadcast(self%gamma0)
    call broadcast(self%omega0)
    call broadcast(self%phi_ext)
    call broadcast(self%source0)
    call broadcast(self%source_option)
    call broadcast(self%t0)

    call broadcast(self%exist)
  end subroutine broadcast_source_config

  !> Gets the default name for this namelist
  function get_default_name_source_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_source_config
    get_default_name_source_config = "source_knobs"
  end function get_default_name_source_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_source_config()
    implicit none
    logical :: get_default_requires_index_source_config
    get_default_requires_index_source_config = .false.
  end function get_default_requires_index_source_config

  !---------------------------------------
  ! Following is for the dist_fn_species config_type
  !---------------------------------------
  
  !> Reads in the dist_fn_species_knobs namelist and populates the member variables
  subroutine read_dist_fn_species_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(dist_fn_species_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    real :: bakdif, fexpi, fexpr

    namelist /dist_fn_species_knobs/ bakdif, fexpi, fexpr

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    bakdif = self%bakdif
    fexpi = self%fexpi
    fexpr = self%fexpr

    ! Now read in the main namelist
    call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index, exist)
    if (exist) read(in_file, nml = dist_fn_species_knobs)
    close(unit = in_file)

    ! Now copy from local variables into type members
    self%bakdif = bakdif
    self%fexpi = fexpi
    self%fexpr = fexpr

    self%exist = exist
  end subroutine read_dist_fn_species_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_dist_fn_species_config(self, unit)
    implicit none
    class(dist_fn_species_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("bakdif", self%bakdif, unit_internal)
    call self%write_key_val("fexpi", self%fexpi, unit_internal)
    call self%write_key_val("fexpr", self%fexpr, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_dist_fn_species_config

  !> Resets the config object to the initial empty state
  subroutine reset_dist_fn_species_config(self)
    class(dist_fn_species_config_type), intent(in out) :: self
    type(dist_fn_species_config_type) :: empty
    select type (self)
    type is (dist_fn_species_config_type)
       self = empty
    end select
  end subroutine reset_dist_fn_species_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_dist_fn_species_config(self)
    use mp, only: broadcast
    implicit none
    class(dist_fn_species_config_type), intent(in out) :: self
    call broadcast(self%bakdif)
    call broadcast(self%fexpi)
    call broadcast(self%fexpr)

    call broadcast(self%exist)
  end subroutine broadcast_dist_fn_species_config

  !> Gets the default name for this namelist
  function get_default_name_dist_fn_species_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_dist_fn_species_config
    get_default_name_dist_fn_species_config = "dist_fn_species_knobs"
  end function get_default_name_dist_fn_species_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_dist_fn_species_config()
    implicit none
    logical :: get_default_requires_index_dist_fn_species_config
    get_default_requires_index_dist_fn_species_config = .true.
  end function get_default_requires_index_dist_fn_species_config
end module dist_fn
