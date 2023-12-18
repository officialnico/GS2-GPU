!> Deals with most types of equilibria to generate GS2's expected geometrical parameters
!> Input:
!> rhoc :: minor radius of surface of interest (see irho below)
!> use_large_aspect :: default false.  Choose true to get "s-alpha" equilibria (see local_eq)
!> Equilibrium options : Choose 1 of:
!>          local_eq : Analytic parameterization of local equilibrium
!>          ppl_eq = .true. :: use PPL style NetCDF equilibrium (Menard)
!>          transp_eq = .true. :: use TRXPL style NetCDF equilibrium (Menard)
!>          gen_eq = .true. :: use GA style NetCDF equilibrium (TOQ)
!>          chs_eq = .true. ::  Use CHEASE generated equilibrium
!>          efit_eq = .true. :: use EFIT  equilibrium
!>          dfit_eq = .true. :: use dipole equilibrium
!> -------------------> [note: if any of the above except local_eq are true set eqfile = input file name]
!> irho :: choose flux surface label
!>         case (1) :: sqrt(toroidal flux)/sqrt(toroidal flux of LCFS)
!>         case (2) :: diameter/(diameter of LCFS).  recommended
!>         case (3) :: poloidal flux/(poloidal flux of LCFS)
!>         case (4) :: rho_mid? Only used for the vacuum ring dipole, dfit_eq = T
!>
!> equal_arc :: .true. makes grad_par coefficient a constant.
!>
!> New feature:  Given an actual equilibrium, allow s_hat and d pressure/d rho
!> to vary self-consistently (a la Greene and Chance, Bishop).  We use the
!> Bishop formalism [C. M. Bishop, et al., NF, Vol. 24, 1579 (1984).]
!>
!> bishop :: uses Bishop relations to generate coefficients.
!>        case(0) :: do not use Bishop relations -- primarily for testing
!>        case(1) :: use Bishop relations with actual equilibrium shat, p'
!>        case(2) :: use Bishop relations with new shat, and p' from p_prime_input
!>        case(3) :: use Bishop relations with new shat, L_p from below
!> s_hat_input :: new magnetic shear, rho/q*dq/drho
!> invLp_input :: new -1/pressure * d pressure/d rho
!>        case(4) :: use Bishop relations with new shat, beta' from below
!> s_hat_input :: new magnetic shear, rho/q*dq/drho
!> beta_prime_input :: new d beta/d rho
!>        case(5) :: use Bishop relations with new shat, alpha from below
!> s_hat_input :: new magnetic shear, rho/q*dq/drho
!> alpha_input :: new alpha = -q**2 * rmaj * d pressure/d rho
!>        case(6) :: use Bishop relations with equilibrium shat, beta' from below
!> beta_prime_input :: new d beta/d rho
!>        case(7) :: use Bishop relations with equilibrium shat, beta' from below
!> dp_mult :: d pressure/d rho = equilibrium gradient * dp_mult
!>         (In each case, the definition of rho is determined by irho above)
!>
!> nperiod :: number of cells of width 2 pi in theta.
!>
!> geoType :: determines the analytic geometry specification
!>         case (0) :: the tilted Miller geometry: generalization of Miller et al. PoP 1998
!>           (see Chapter 3.2.1 of J. Ball MIT Masters thesis)
!>         case (1) :: the global MHD solution geometry: calculated by inverting the large
!>           aspect ratio, constant-current MHD solution
!>           (see J. Ball "Optimized up-down asymmetry to drive fast intrinsic rotation in tokamak reactors")
!>         case (2) :: the generalized ellipticity geometry: the polar formula for an ellipse generalized
!>           to arbitrary shaping effect mode number (see Chapter 5.1.2 of J. Ball Oxford PhD thesis)
!>         case (3) :: the Fourier series geometry: Fourier decomposition of the minor radius as a function
!>           of poloidal angle (see Chapter 5.1.1 of J. Ball Oxford PhD thesis)
!> rmaj :: major radius of magnetic axis (in units of aminor) if not local_eq
!>      or major radius of center of flux surface of interest (units of aminor) otherwise
!> r_geo :: major radius of last closed flux surface or
!>       if local_eq: sets the reference magnetic field to be the value of the toroidal
!>       field at R=r_geo on the flux surface of interest
!> aSurf :: the real-space flux surface label of the surface with deltam and deltan shaping as well as
!>   a major radius of Rmaj (only used when geoType={1})
!> shift :: d rmaj/drho.  (Typically < 0.) (only used when geoType={0,2,3})
!> -> Alternative definitions:
!>   Raxis :: major radial location of the magnetic axis. (only used when geoType={1})
!> shiftVert :: d zmaj/drho. (< 0 for a vertical shift) (only used when geoType={0,2,3})
!> -> Alternative definitions:
!>   Zaxis :: axial location of the magnetic axis. (only used when geoType={1})
!> mMode :: the first poloidal mode number (only used when geoType={2,3})
!> nMode :: the second poloidal mode number (only used when geoType={2,3})
!> deltam :: magnitude of the first poloidal shaping effect (only used when geoType={2,3})
!> -> Alternative definitions:
!>   delta2 :: identical to akappa (only used when geoType={1})
!>   akappa :: flux surface elongation (only used when geoType={0})
!> deltan :: magnitude of the second poloidal shaping effect (only used when geoType={2,3})
!> -> Alternative definitions:
!>   delta3 :: magnitude of the triangulaity-like shaping effect (only used when geoType={1})
!>   tri :: flux surface triangulaity (only used when geoType={0})
!> deltampri :: d(deltam)/d(rhoc) (only used when geoType={2,3})
!> -> Alternative definitions:
!>   akappapri :: d(akappa)/d(rhoc) (only used when geoType={0})
!> deltanpri :: d(deltan)/d(rhoc) (only used when geoType={2,3})
!> -> Alternative definitions:
!>   tripri :: d(tri)/d(rhoc) (only used when geoType={0})
!> thetam :: tilt angle of the first poloidal shaping effect (only used when geoType={2,3})
!> -> Alternative definitions:
!>   theta2 :: tilt angle of the flux surface elongation (only used when geoType={1})
!>   thetak :: tilt angle of the flux surface elongation (only used when geoType={0})
!> thetan :: tilt angle of the first poloidal shaping effect (only used when geoType={2,3})
!> -> Alternative definitions:
!>   theta3 :: tilt angle of the triangularity-like shaping effect (only used when geoType={1})
!>   thetad :: tilt angle of the flux surface triangularity (only used when geoType={0})
!> qinp :: safety factor q at surface of interest (local_eq = T only)
!> shat :: d q/drho at surface of interest
!>
!> delrho :: numerical parameter.  Should be "small enough".  Typically 0.001 ok.
!>
!> force_sym :: true -> assume up-down symmetric equilibrium.
!> writelots :: .true. for more output
!>
!>  Use:
!> Include this module, set input variables appropriately and call eikcoefs.
!>
!> I do not recommend calling other routines, but a few are available for
!> diagnostic purposes.
module geometry
  use geo_utils, only: abstract_geo_type, EQFILE_LENGTH, flux_surface_type
  implicit none

  private

  public :: beta_a_fun, eikcoefs, geofax, pbarofrho, Rpos, rfun, diameter
  public :: qfun, rpofrho, finish_geometry, verb, big, pbarfun, geom
  public :: nperiod, eqfile, eqnormfile, irho, surf, delrho
  public :: ppl_eq, transp_eq, chs_eq, efit_eq, gs2d_eq, idfit_eq, dfit_eq, local_eq, gen_eq
  public :: force_sym, writelots, equal_arc, bishop, use_large_aspect
  public :: p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
  public :: eikcoefs_output_type, run_eikcoefs_and_resample, get_theta_eqarc

  class(abstract_geo_type), allocatable :: geom
  ! Start of inputs
  type(flux_surface_type) :: surf
  real :: delrho, p_prime_input, invLp_input, beta_prime_input, alpha_input, dp_mult
  integer :: nperiod, irho, bishop
  logical :: gen_eq, efit_eq, ppl_eq, local_eq, chs_eq, dfit_eq, idfit_eq, gs2d_eq, transp_eq
  logical :: writelots, equal_arc, force_sym, use_large_aspect
  character(len = EQFILE_LENGTH) :: eqfile, eqnormfile
  ! End of inputs
  !> We always set rpmin = psi_0 and rpmax = psi_a so we probably don't need both.
  real :: psi_0, psi_a, rpmin, rpmax
  real :: rmaj, r_geo, rhoc, qinp !< From surf
  integer :: big = 8

  !> The single-period theta domain "half-grid" size. That is, the number of points in
  !> \((0, \pi]\). `theta(-nth:nth)` is the single central period \([-\pi, \pi]\).  For
  !> EFIT and local equilbria, this is set from
  !> [[theta_grid_parameters_config_type:ntheta]]; for all other numerical equilbria, it
  !> is determined by the data file. Always even.
  integer :: nth
  !> The full theta domain "half-grid" size. That is, the total size of the
  !> [[geometry:theta]] grid is `(2 * ntgrid) + 1`, with array bounds
  !> `theta(-ntgrid:ntgrid)`. Equal to `(2*nperiod - 1)*nth`
  integer :: ntgrid
  !> Equal to `2 * nth`
  !> FIXME: Remove and replace uses with calls to [[periodic_copy]]
  integer :: ntheta
  !> Verbosity of print statements | Should we use runtime_tests:verbosity?
  integer :: verb = 2
  !> If debug then print out lots of debug info.
  !> debug is true if (verb >= 3)
  logical :: debug = .false.

  !> A type to hold a collection of geometrical quantities as calculated by
  !> [[eikcoefs]]
  type eikcoefs_output_type
     integer :: ntheta
     real :: s_hat_new
     real :: drhodpsi
     real :: kxfac
     real :: qsf
     real :: surfarea
     real :: dvdrhon
     real :: shat
     real :: dbetadrho
     real :: rhoc
     !> theta :: the theta grid of the results.  The final grid has
     !>          pi - 2*pi*nperiod <= theta <= 2*pi*nperiod - pi
     real, dimension(:), allocatable :: theta
     !> bmag :: |B|
     real, dimension(:), allocatable :: bmag
     !> gradpar :: coefficient of b_hat.grad operator
     real, dimension(:), allocatable :: gradpar
     !> grho :: grad(rho)
     real, dimension(:), allocatable :: grho
     !> cdrift :: part of coriolis drift operator independent of theta_0
     real, dimension(:), allocatable :: cdrift
     !> cvdrift  :: coefficient of v_par**2/Omega b_hat x (b_hat.grad b_hat).grad
     real, dimension(:), allocatable :: cvdrift
     !> gbdrift  :: coefficient of mu/Omega b_hat x (grad B).grad operator
     real, dimension(:), allocatable :: gbdrift
     !> cdrift0 :: part of coriolis drift operator proportional to theta_0
     real, dimension(:), allocatable :: cdrift0
     !> cvdrift0 :: part of total curvature drift operator proportional to theta_0
     real, dimension(:), allocatable :: cvdrift0
     !> gbdrift0 :: part of total grad B drift operator proportional to theta_0
     real, dimension(:), allocatable :: gbdrift0
     !> gds2 :: part of grad_perp**2 operator independent of theta_0
     real, dimension(:), allocatable :: gds2
     !> gds21 :: part of grad_perp**2 operator proportional to theta_0
     real, dimension(:), allocatable :: gds21
     !> gds22 :: part of grad_perp**2 operator proportional to theta_0**2
     real, dimension(:), allocatable :: gds22
     !> gds23 :: part of v_E . grad theta operator independent of theta_0
     real, dimension(:), allocatable :: gds23
     !> gds24 :: part of v_E . grad theta operator proportional to theta_0
     real, dimension(:), allocatable :: gds24
     !> gds24_noq :: gds24/dqdrp
     real, dimension(:), allocatable :: gds24_noq
     !> rplot :: R(theta) - major radius on theta grid
     real, dimension(:), allocatable :: rplot
     !> zplot :: Z(theta) - height on theta grid
     real, dimension(:), allocatable :: zplot
     !> aplot :: a(theta) - alpha-phi on theta grid
     real, dimension(:), allocatable :: aplot
     !> rprime :: derivative wrt rho of R(theta) - major radius on theta grid
     real, dimension(:), allocatable :: rprime
     !> zprime :: derivative  wrt rho Z(theta) - height on theta grid
     real, dimension(:), allocatable :: zprime
     !> aprime :: derivative  wrt rho a(theta) - alpha-phi on theta grid
     real, dimension(:), allocatable :: aprime
     !> bpol :: Poloidal magnetic field strength
     real, dimension(:), allocatable :: bpol
  end type eikcoefs_output_type

contains

  !> Run eikcoefs with `ntheta = ntheta_input` and then resample onto
  !> a theta grid with resolution `ntheta_target`.
  subroutine run_eikcoefs_and_resample(ntheta_input, ntheta_target, results, job_id)
    implicit none
    !> This is the number of theta grid points per two pi period used in the
    !> geometry calculation (actually one less than this).
    integer, intent(in out) :: ntheta_input
    !> This is the desired number of theta grid points per two pi period
    !> in the resulting geometrical coefficients used elsewhere.
    integer, intent(in) :: ntheta_target
    type(eikcoefs_output_type), intent(out) :: results
    integer, intent(in), optional :: job_id
    type(eikcoefs_output_type) :: full_results
    real, dimension(:), allocatable :: input_index_space, target_index_space
    integer :: ig, ntgrid_target, ntgrid_input
    real :: spacing

    ! First run eikcoefs using requested resolution
    call eikcoefs(ntheta_input, full_results, job_id)

    ntgrid_input = (size(full_results%theta) - 1) / 2
    ntgrid_target =(ntheta_target*(2*nperiod-1)) / 2

    ! Set scalar outputs
    results%ntheta = ntheta_target
    results%s_hat_new = full_results%s_hat_new
    results%drhodpsi = full_results%drhodpsi
    results%kxfac = full_results%kxfac
    results%qsf = full_results%qsf
    results%surfarea = full_results%surfarea
    results%dvdrhon = full_results%dvdrhon
    results%shat = full_results%shat
    results%dbetadrho = full_results%dbetadrho
    results%rhoc = full_results%rhoc

    ! Make an index space array for the original theta
    ! grid. Spans [0,1] in size(theta)=2*ntgrid_input points
    allocate(input_index_space(-ntgrid_input:ntgrid_input))
    spacing = 1.0/(2*ntgrid_input)
    do ig = -ntgrid_input, ntgrid_input
       input_index_space(ig) = (ig+ntgrid_input) * spacing
    end do

    ! Make new theta grid
    allocate(target_index_space(-ntgrid_target:ntgrid_target))

    ! First make target index space grid
    spacing = 1.0/(2*ntgrid_target)
    do ig = -ntgrid_target, ntgrid_target
       target_index_space(ig) = (ig+ntgrid_target) * spacing
    end do

    ! Now interpolate theta grid onto target index space grid
    call resample( input_index_space, full_results%theta, target_index_space, results%theta)

    ! Now interpolate other output quantities onto target theta grid
    call resample( full_results%theta, full_results%bmag, results%theta, results%bmag)
    call resample( full_results%theta, full_results%gradpar, results%theta, results%gradpar)
    call resample( full_results%theta, full_results%grho, results%theta, results%grho)
    call resample( full_results%theta, full_results%cdrift, results%theta, results%cdrift)
    call resample( full_results%theta, full_results%cvdrift, results%theta, results%cvdrift)
    call resample( full_results%theta, full_results%gbdrift, results%theta, results%gbdrift)
    call resample( full_results%theta, full_results%cdrift0, results%theta, results%cdrift0)
    call resample( full_results%theta, full_results%cvdrift0, results%theta, results%cvdrift0)
    call resample( full_results%theta, full_results%gbdrift0, results%theta, results%gbdrift0)
    call resample( full_results%theta, full_results%gds2, results%theta, results%gds2)
    call resample( full_results%theta, full_results%gds21, results%theta, results%gds21)
    call resample( full_results%theta, full_results%gds22, results%theta, results%gds22)
    call resample( full_results%theta, full_results%gds23, results%theta, results%gds23)
    call resample( full_results%theta, full_results%gds24, results%theta, results%gds24)
    call resample( full_results%theta, full_results%gds24_noq, results%theta, results%gds24_noq)
    call resample( full_results%theta, full_results%rplot, results%theta, results%rplot)
    call resample( full_results%theta, full_results%zplot, results%theta, results%zplot)
    call resample( full_results%theta, full_results%aplot, results%theta, results%aplot)
    call resample( full_results%theta, full_results%rprime, results%theta, results%rprime)
    call resample( full_results%theta, full_results%zprime, results%theta, results%zprime)
    call resample( full_results%theta, full_results%aprime, results%theta, results%aprime)
    call resample( full_results%theta, full_results%bpol, results%theta, results%bpol)

  contains
    subroutine resample(x_orig, y_orig, x_new, y_new)
      use splines, only: new_spline, delete_spline, splint, spline
      real, dimension(:), intent(in) :: x_orig, y_orig
      real, dimension(-ntgrid_target:), intent(in) :: x_new
      real, dimension(:), allocatable, intent(out) ::  y_new
      type(spline) :: the_spline
      integer :: i, lower, upper

      ! Spline the old data
      call new_spline(size(x_orig), x_orig, y_orig, the_spline)

      lower = lbound(x_new, 1)
      upper = ubound(x_new, 1)

      allocate(y_new(lower:upper))

      do i = lower, upper
         y_new(i) = splint(x_new(i), the_spline)
      end do

      ! Get rid of the spline
      call delete_spline(the_spline)
    end subroutine resample

  end subroutine run_eikcoefs_and_resample

  !> FIXME : Add documentation
  subroutine eikcoefs (ntheta_input, outputs, job_id_in)
    use constants, only: twopi
    use file_utils, only: open_output_file, close_output_file
    use file_utils, only: futils_initialized => initialized, init_file_utils
    use file_utils, only: finish_file_utils
    use optionals, only: get_option_with_default
    implicit none
    !> For most numerical equilibria, **returns** the value of [[geometry:nth]];
    !> for EFIT and local equilibria, instead **sets** [[geometry:nth]].
    !>
    !> If setting [[geometry:nth]], then even values are used directly, while odd values
    !> set `nth = ntheta_input - 1`
    integer, intent (inout) :: ntheta_input
    !> A type containing all of the outputs of eikcoefs
    type(eikcoefs_output_type), intent(out) :: outputs
    !> job_id is the job id in a  Trinity or list mode run. Used
    !> for debug statements
    integer, intent(in), optional :: job_id_in
    logical :: start_file_utils, dummy_flag
    integer :: job_id
    ! Local variables:
    real, allocatable, dimension(:)   :: seik, rgrid, bpolmag, invR, trip, bmod, ds, jacob
    real, allocatable, dimension(:)   :: dsdthet, rmajor, dummy, th_bish, invRpol, ltheta
    real, allocatable, dimension(:,:) :: thgrad, rpgrad, crpgrad, grads
    real, allocatable, dimension(:,:) :: gradrptot, gradstot, gradztot, gradthtot
    real :: shat, dbetadrho, dp_new, di_new, s_hat_new, dqdrp, B_T0, qval, dum, bi
    real :: rp, pbar, drhodrp, drhodrhod, drhodpsi, dpsidrp, dpsidrho
    ! Used to translate definition of the poloidal angle such that the location
    ! of the maximium magnetic field is at pi (see section 3.2.3 of Ball MIT
    ! Masters thesis or section 2.4 of "GS2 analytic geometry specification")
    real :: thetaShift
    !> For outputting geom values
    integer :: i, iBmax, fort11_unit

    ! Unpack from surf into module level for now
    rhoc = surf%r ; qinp = surf%q
    rmaj = surf%rmaj ; r_geo = surf%r_geo
    surf%nt = ntheta_input ; surf%dr = delrho
    !Set output verbosity levels
    job_id = get_option_with_default(job_id_in, -1)
    start_file_utils = .not. futils_initialized

    if (start_file_utils) then
      ! Initialize file utils with a run name we provide
      ! and by setting trin_run prevent it from searching
      ! the command line arguments.
      write(*,*) 'WARNING: Automatically initializing file_utils&
      & in eikcoefs. This is a hack to &
      & support older programs like eiktest. You should&
      & initialize file utils before calling eikcoefs'
      call init_file_utils(dummy_flag, name='geometry_module', trin_run=.true., input=.false.)
    end if
    !Create module level output file and get unit
    call open_output_file(fort11_unit,".fort.11")

    !Compute the initial constants
    debug = (verb > 2)
    if (debug) writelots = .true.

    if (debug) write(6,*) "eikcoefs: local_eq=",local_eq, 'jid=',job_id

    if (local_eq) then
       if (irho /= 2) write(*,*) 'Forcing irho = 2'
       irho = 2
       if (bishop == 0) then
          write(*,*) 'Forcing bishop = 1'
          bishop = 1 !Not clear that this is a good choice for Miller
       end if
    end if

    if (debug) write(6,*) "eikcoefs: call check", 'jid=',job_id
    call check(gen_eq, efit_eq, ppl_eq, local_eq, dfit_eq, idfit_eq, chs_eq, gs2d_eq, transp_eq, bishop)

    thetaShift = 0.0
    if (.not. allocated(geom)) call allocate_geom(geom)
    call initialise_geometry_and_get_bmod

    ! Start theta shift calc. This is used to modify the definition of
    ! theta such that the maxiumum magnetic field is at pi. This is
    ! necessary to correctly calculate bounce points. Can only shift with local_eq
    if (local_eq) then
       ! Locate the index of the maximum magnetic field in the range [-pi,pi)
       iBmax = maxloc(Bmod(-nth:nth-1), dim = 1) - (nth + 1)

       ! If max is not at -pi and theta then shift theta.
       if (iBmax /= -nth) then
          thetaShift = outputs%theta(iBmax) - outputs%theta(-nth)
          call geom%finalise !Clean up geometry to allow reinitialisation
          write(*,*) "Bmax is not located at +/-pi, translating theta by",thetaShift," so that it is..."
          call initialise_geometry_and_get_bmod
       end if
    end if

    !     compute |grad rho|:
    allocate(outputs%grho(-ntgrid:ntgrid))
    do i = -nth, nth
       outputs%grho(i) = sqrt(rpgrad(i,1)**2 + rpgrad(i,2)**2) * abs(drhodrp)
    end do
    call periodic_copy(outputs%grho, 0., nperiod)
    if (force_sym) call sym(outputs%grho, 0, ntgrid)

    if (writelots) write(fort11_unit,*) 'dpsidrp= ',dpsidrp
    drhodpsi = drhodrp / dpsidrp ; dpsidrho = 1.0 / drhodpsi

    !     compute derivatives of eikonal
    if (bishop /= 0) then
       call bishop_grads
    else
       call bishop_large_aspect_ratio_approximation
    end if

    if (writelots) then
       write(fort11_unit,*) 'i,theta,grads1,2=  '
       do i = -ntgrid, ntgrid
          write(fort11_unit,*) i, outputs%theta(i), grads(i, 1), grads(i, 2)
       end do
    end if
    !     compute total grad S including toroidal part
    call gradstottgrid(rmajor, grads, gradstot)

    gradrptot(-nth:nth, 1:2) = rpgrad(-nth:nth, :)
    call periodic_copy(gradrptot(:, 1), 0., nperiod)
    call periodic_copy(gradrptot(:, 2), 0., nperiod)
    gradrptot(:, 3) = 0.

    ! gradthtot is grad theta since ballooning, theta goes from
    ! -ntheta*nperiod -> ntheta*nperiod, but grad theta is periodic in theta
    gradthtot(-nth:nth, 1:2) = thgrad(-nth:nth, :)
    call periodic_copy(gradthtot(:, 1), 0., nperiod)
    call periodic_copy(gradthtot(:, 2), 0., nperiod)
    gradthtot(:, 3) = 0.

    !     compute the magnitude squared of grad S
    allocate(outputs%gds2(-ntgrid:ntgrid), outputs%gds21(-ntgrid:ntgrid), &
         outputs%gds22(-ntgrid:ntgrid), outputs%gds23(-ntgrid:ntgrid), &
         outputs%gds24(-ntgrid:ntgrid), outputs%gds24_noq(-ntgrid:ntgrid))
    call calculate_metric_terms(gradstot, gradrptot, gradthtot, dpsidrho, dqdrp, &
         qval, rhoc, bmod, force_sym, outputs%gds2, outputs%gds21, outputs%gds22, &
         outputs%gds23, outputs%gds24, outputs%gds24_noq)

    !     compute curvature, grad b, and coriolis drift terms
    allocate(outputs%gbdrift(-ntgrid:ntgrid), outputs%gbdrift0(-ntgrid:ntgrid))
    allocate(outputs%cvdrift(-ntgrid:ntgrid), outputs%cvdrift0(-ntgrid:ntgrid))
    allocate(outputs%cdrift(-ntgrid:ntgrid), outputs%cdrift0(-ntgrid:ntgrid))
    call drift(rgrid, rp, gradstot, gradrptot, gradztot, dqdrp, dpsidrho, dpsidrp,&
         Bmod, Bpolmag, invRpol, invR, th_bish, ltheta, outputs%theta, rmajor, &
         dp_new, bi, force_sym, outputs%gbdrift, outputs%gbdrift0, outputs%cvdrift, &
         outputs%cvdrift0, outputs%cdrift, outputs%cdrift0)

    !     compute the magnetic field along theta
    allocate(outputs%bmag(-ntgrid:ntgrid))
    if (use_large_aspect) then
       outputs%bmag(-nth:nth)= 1 / (1 + rgrid(-nth:nth) * cos(outputs%theta(-nth:nth)) / rmaj)
       call periodic_copy(outputs%bmag, 0., nperiod)
    else
       outputs%bmag = Bmod
    end if
    if (force_sym) call sym(outputs%bmag, 0, ntgrid)

    !     compute the gradparallel coefficient
    allocate(outputs%gradpar(-ntgrid:ntgrid))
    if (dfit_eq) then
       outputs%gradpar = trip / outputs%bmag
    else
       outputs%gradpar(-nth:nth) = -bi / (rmajor(-nth:nth)**2 * outputs%bmag(-nth:nth) * dsdthet(-nth:nth))
       call periodic_copy(outputs%gradpar, 0., nperiod)
    end if
    if (force_sym) call sym(outputs%gradpar, 0, ntgrid)

    ! Note : We should call this before we potentially redefine the theta grid
    ! due to equal_arc below as some specific equilibrium types rely on theta
    ! having a consistent definition with that used to calculate rgrid initially.
    ! In particular this is important for Efit and GS2D equilibria (in eeq). Other
    ! equilibrium types use interpolation so shouldn't depend on the specific theta grid
    ! (although it's possible there could still be some interaction, given that one of the
    ! outputs of the below 'aplot' is given by `seik` which is evaluated on the original
    ! grid rather than the "new" one defined below).
    ! @note Check if restriction on cartesian geo types still holds.
    allocate(outputs%rplot(-ntgrid:ntgrid), outputs%rprime(-ntgrid:ntgrid))
    allocate(outputs%zplot(-ntgrid:ntgrid), outputs%zprime(-ntgrid:ntgrid))
    allocate(outputs%aplot(-ntgrid:ntgrid), outputs%aprime(-ntgrid:ntgrid))
    call plotdata (rgrid, seik, grads, dpsidrho, outputs%theta, &
         outputs%rplot, outputs%zplot, outputs%aplot, &
         outputs%rprime, outputs%zprime, outputs%aprime)

    ! In the lines below we redefine theta to be theta', and gradpar
    ! to be gradpar'. This means that the values of theta in the gs2
    ! output file are not values of the poloidal angle but values of theta_prime
    if (equal_arc) then
       call arclength (ntheta, nperiod, outputs%theta, outputs%gradpar, dummy)
       outputs%theta = dummy
    end if

    call calculate_surface_area(outputs%theta, outputs%gradpar, outputs%bmag, outputs%grho, &
         dpsidrho, rhoc, nth, outputs%surfarea, outputs%dvdrhon)
    if (writelots) write (fort11_unit,*) 'surfarea=', outputs%surfarea

    !Close output file
    call close_output_file(fort11_unit)

    outputs%kxfac = abs(qval / (rhoc * drhodpsi))

    if (qval == 0. .and. (dfit_eq .or. idfit_eq)) then
       shat = max(1.e-7, shat)
       !RN> I think kxfac should include the sign. However, from GS2 convention,
       !RN> eventually kxfac = abs(kxfac). See get_grids in theta_grid. Even if the sign is
       !RN> effective, probably it has no effect due to symmetry. I just leave the sign here.
       outputs%kxfac = - 1 / drhodpsi
       outputs%gds22 = outputs%grho**2 / drhodpsi**2
    end if

    if (start_file_utils) call finish_file_utils

    ! Pack module/routine level data into output type
    outputs%ntheta = ntheta_input
    outputs%bpol = bpolmag

    outputs%s_hat_new = s_hat_new
    outputs%drhodpsi = drhodpsi
    outputs%qsf = qval
    outputs%shat = shat
    outputs%dbetadrho = dbetadrho
    outputs%rhoc = rhoc

  contains
    !> Initialises the specific geometry backend
    subroutine initialise_geometry_and_get_bmod
      use geo_utils, only: geo_input_type
      type(geo_input_type) :: inputs
      real :: avgrmid, rbar
      character(1) :: char
      integer :: fort78_unit, fort99_unit
      ! Initialise the specific geometry types (mostly reading data
      ! from file and calculating derived quantities).
      if (debug) write(6,*) "eikcoefs: chs_eq=", chs_eq, 'jid=',job_id

      if (local_eq) then
         inputs%surf = surf ; inputs%thShift = thetaShift
      else
         inputs%eqfile = eqfile
         inputs%big = big
         inputs%eqnormfile = eqnormfile
      end if
      call geom%initialise(inputs, psi_0, psi_a, rmaj, B_T0, avgrmid)
      rpmin = psi_0 ; rpmax = psi_a
      if (debug) write(6,*) "eikcoefs: rpmin,rpmax=", rpmin, rpmax
      call init_uniform_theta_grid(ntheta_input, outputs%theta)
      if (debug) write(6,*) 'eikcoefs: alloc_local_arrays'
      if (.not. allocated(seik)) call alloc_local_arrays(ntgrid)

      ! Calculate grad(psi), grad(BT) and grad(theta) for the
      ! specific geometry class we're using.
      call geom%finish_setup

      if (debug) write(6,*) "eikcoefs: if (writelots)"
      if (writelots) then
         write(fort11_unit,*) 'All lengths normalized to the avg midplane minor radius'
         if (irho /= 2) then
            write(fort11_unit,*) 'except rho, which always goes from 0 to 1.'
            write(fort11_unit,*) 'The rho normalization is determined by irho.'
         end if
         write(fort11_unit,*) 'avgrmid = ',avgrmid,' meters'
         if (.not. local_eq) write(fort11_unit,*) 'R_mag = ',rmaj
         if (.not. local_eq) write(fort11_unit,*) 'B_T0 = ',B_T0
      end if

      if (debug) write(6,*) "eikcoefs: find rp", 'jid=',job_id

      if (local_eq) then
         ! Because irho==2 if local_eq and diameter is effectively analytic
         ! one would expect rpofrho(rho) == rho in this case. So here we
         ! effectively just have an optimisation (I think)
         rp = rhoc
      else
         rp = rpofrho(rhoc, fort11_unit)
         R_geo = rcenter(rpmax) !Only used for writelots/checks
      end if
      if (debug) write(6,*) "eikcoefs: rp=",rp, 'jid=',job_id

      if (debug) write(6,*) "eikcoefs: find rgrid"
      if (efit_eq .or. dfit_eq .or. gs2d_eq) then !Cartesian like
         call rtg(rgrid, rp, outputs%theta)
      else
         rgrid = rp !This is likely just an optimisation, rtg should work here as well?
      end if

      pbar = pbarfun(rgrid(0), outputs%theta(0))
      drhodrp = drho_drp(rp)

      if ((.not. local_eq) .and. R_geo > Rmaj + 1.e-5) then
         write(6,*)
         write(6,*) 'Warning: the center of the LCFS is outboard of the '
         write(6,*) 'center of the reference surface.  This is probably wrong.'
         write(6,*)
      end if

      if (writelots) then
         write(fort11_unit,*) 'rp of rho = ',rp
         write(fort11_unit,*) 'Rcenter(rho) = ',rcenter(rp)
         write(fort11_unit,*) 'R_geo = ',R_geo,' !or ',btori(pbar)

         call open_output_file(fort78_unit,".fort.78")
         write(fort78_unit,*) '# theta_pol,  R,   Z'
         do i = -nth, nth
            write(fort78_unit,*) outputs%theta(i), Rpos(rgrid(i),outputs%theta(i)), &
                 Zpos(rgrid(i),outputs%theta(i))
         end do
         call close_output_file(fort78_unit)

         call open_output_file(fort99_unit,".fort.99")
         do i = -nth, nth
            write(fort99_unit,*) 'i=',i,' thet= ',outputs%theta(i),' rgrid= ',rgrid(i)
         end do
         call close_output_file(fort99_unit)

         rbar = 0.
         do i = -nth, nth - 1
            rbar = rbar + rgrid(i) * (outputs%theta(i + 1) - outputs%theta(i))
         end do
         rbar = rbar / twopi
         write (*,*) 'r_bar = ', rbar
         write (*,*)

         write(fort11_unit,*) 'drhodrp=', drhodrp
         if (irho == 1) then
            drhodrhod = drho_drhod(rp, drhodrp)
            write(fort11_unit,*) 'drho_drhod=', drhodrhod
         end if
      end if

      if (debug) write(fort11_unit,*) 'eikcoefs: various logs', gen_eq, ppl_eq, chs_eq, efit_eq, dfit_eq, idfit_eq, gs2d_eq, transp_eq

      if (debug) write(6,*) "eikcoefs: call rmajortgrid"
      call rmajortgrid(rgrid, outputs%theta, rmajor)
      call periodic_copy(rmajor, 0., nperiod)
      call invRtgrid(rgrid, outputs%theta, invR)
      call periodic_copy(invR, 0., nperiod)

      ! Get the toroidal field function (i.e. fp = R B_T)
      bi = btori(pbar)

      !     compute gradient of rp
      char = 'P'
      ! note that rpgrad is given in (R,Z,phi) coordinates if bishop=0
      ! and (rho,l,phi) coordinates if bishop>0
      ! also, rpgrad is grad psi if not local_eq and grad rho if local_eq -- MAB
      ! When bishop > 0
      ! rpgrad(:,1) is |grad rp|
      ! rpgrad(:,2) is 0
      ! crpgrad(:,1) is drp/dR
      ! crpgrad(:,2) is drp/dZ
      !
      ! bishop = 0 is high aspect ratio coefficients for debugging -- EGH
      call geom%gradient(rgrid, outputs%theta, crpgrad, char, dum, nth, ntgrid)
      call periodic_copy(crpgrad(:,1), 0., nperiod)
      call periodic_copy(crpgrad(:,2), 0., nperiod)

      if (bishop == 0) then
         rpgrad = crpgrad
      else
         call geom%bgradient(rgrid, outputs%theta, rpgrad, char, dum, nth, ntgrid)
      end if

      !     compute  coordinate gradients
      if (debug) write(6,*) "eikcoefs: call thetagrad"
      call thetagrad(rgrid, outputs%theta, thgrad)
      if (debug) write(6,*) "eikcoefs: done thetagrad"
      call tripprod2dtgrid(rpgrad, thgrad, invR, trip)
      call periodic_copy(trip, 0., nperiod)

      !     compute eikonal S
      if (local_eq) qval = qfun(pbarofrho(rhoc, fort11_unit)) !leq ignores input pbar
      if (debug) write(6,*) "eikcoefs: call eikonal"
      call eikonal(outputs%theta, invR, trip, qval, bi, seik, dsdthet, dpsidrp)
      if (local_eq) trip = trip * dpsidrp
      if (debug) write(6,*) "eikcoefs: done eikonal"
      if (writelots) write(fort11_unit,*) 'q= ', qval

      if (bishop > 0) then !Following bits _could_ be calculated for all bishop
         call B_pol(invR, crpgrad, Bpolmag, nth) !<Should this use rpgrad for local?
         if (local_eq) Bpolmag(-nth:nth) = Bpolmag(-nth:nth) * dpsidrp
         call periodic_copy(Bpolmag, 0., nperiod)
         call B_mod(invR, Bpolmag, Bmod, bi, nth)
         call periodic_copy(Bmod, 0., nperiod)
      end if !Note bmod not set for bishop == 0
    end subroutine initialise_geometry_and_get_bmod

    subroutine bishop_grads
      logical, parameter :: save_bishop = .true.
      integer :: bishop_unit
      real :: a_b, b_b, c_b, dp, di, s_hat, pressure, tmp
      real, allocatable, dimension(:) :: a_bish, b_bish, c_bish, da_bish, db_bish, dc_bish, &
           zoftheta, dzdl, dsdl

      call th_bishop(crpgrad, th_bish, nth)
      call periodic_copy(th_bish, 0., nperiod)
      call loftheta(rgrid, outputs%theta, ltheta)
      call periodic_copy(ltheta, ltheta(nth) - ltheta(-nth), nperiod)
      call inv_Rpol(outputs%theta, th_bish, ltheta, invRpol, nth)
      call periodic_copy(invRpol, 0., nperiod)

      if (debug) write(6,*) "eikcoefs:  gradients"
      if (debug) write(6,*) 'drhodpsi', drhodpsi
      if (debug) write(6,*) 'bi', bi
      if (verb > 5) write(6,*) 'rmajor', rmajor(:)
      if (verb > 5) write(6,*) 'Bpolmag', Bpolmag(:)
      if (verb > 5) write(6,*) rpgrad(-nth:nth,1)
      if (verb > 5) write(6,*) rpgrad(-nth:nth,2)
      if (verb > 5) write(6,*) thgrad(-nth:nth,1)
      if (verb > 5) write(6,*) thgrad(-nth:nth,2)
      if (verb > 5) write(6,*) 1.0 / trip(-nth:nth)
      if (verb > 5) write(6,*) -1.0 / (bpolmag(-nth:nth) * invRpol(-nth:nth))
      if (debug) write(6,*) "eikcoefs: end gradients"

      allocate(a_bish(-ntgrid:ntgrid), b_bish(-ntgrid:ntgrid), c_bish(-ntgrid:ntgrid), &
           da_bish(-ntgrid:ntgrid), db_bish(-ntgrid:ntgrid), dc_bish(-ntgrid:ntgrid), &
           zoftheta(-ntgrid:ntgrid), dzdl(-ntgrid:ntgrid), dsdl(-ntgrid:ntgrid))

      da_bish = (1 / rmajor**2) + (bi / (Bpolmag * rmajor**2))**2
      db_bish = bi / (Bpolmag * rmajor)**2
      dc_bish = bi * (sin(th_bish) + rmajor * invRpol) / (-Bpolmag * rmajor**4)

      da_bish = da_bish / trip
      db_bish = db_bish / trip
      dc_bish = dc_bish / trip

      call integrate(da_bish, outputs%theta, a_bish, ntgrid)
      call integrate(db_bish, outputs%theta, b_bish, ntgrid)
      call integrate(dc_bish, outputs%theta, c_bish, ntgrid)

      a_b = a_bish(nth) - a_bish(-nth)
      b_b = b_bish(nth) - b_bish(-nth)
      c_b = c_bish(nth) - c_bish(-nth)

      a_bish = - a_bish * Bpolmag * rmajor
      b_bish = - b_bish * Bpolmag * rmajor
      c_bish = - c_bish * Bpolmag * rmajor

      ! find shat, pressure' and report them
      dp = dpfun(pbar) ; di = dbtori(pbar)

      if (abs(qval) > epsilon(0.)) then
         tmp = rhoc / (qval * twopi * drhodpsi)
         s_hat = tmp * (a_b * di + b_b * dp + c_b * 2)
         if (local_eq) s_hat = surf%shat !Can't determine this by calculation
         if (debug) write(6,*) "a_b (f'), b_b (p'), c_b terms =", a_b*di, b_b*dp, c_b*2
         if (debug) write(6,*) 'qval, rho, drhodpsi, s_hat =', qval, rhoc, drhodpsi, s_hat
      else
         tmp = 1. ; s_hat = 0.
      end if

      write(fort11_unit,*)
      write(fort11_unit,*) 'Quantity, equilibrium value, value used'

      ! Decide how we want to set s_hat_new in the default case
      s_hat_new = s_hat
      select case (bishop)
      case (2, 3, 4, 5, 8, 9)
         s_hat_new = surf%shat
      case (6, 7)
         s_hat_new = s_hat
      end select

      ! Decide how we want to set dp_new
      dp_new = dp
      select case (bishop)
      case (2)
         dp_new = p_prime_input
      case (3)
         pressure = pfun(pbar)
         dp_new = -invLp_input * pressure * drhodpsi
      case (4, 6, 9)
         dp_new = 0.5 * beta_prime_input * drhodpsi
      case (5)
         dp_new = -alpha_input * drhodpsi / (qval**2 * rmaj)
      case (7)
         dp_new = dp * dp_mult
      case (8)
         dp_new = 0.5 * beta_prime_input * drhodpsi * dp_mult
      end select

      !Report on updated dp
      select case(bishop)
      case (2)
         write(fort11_unit,*) 'p_prime = ',dp,', ',dp_new
         if (writelots) write(*,*) 'p_prime = ',dp,', ',dp_new
      case (3)
         write(fort11_unit,*) '1/L_p = ',-dp / (drhodpsi * pressure), ', ', &
              -dp_new / (drhodpsi *pressure)
         if (writelots) write(*,*) '1/L_p = ',-dp / (drhodpsi * pressure), ', ', &
              -dp_new / (drhodpsi * pressure)
      case (5)
         write(fort11_unit,*) 'alpha = ', -qval**2 * rmaj * dp / drhodpsi,', ', &
              -qval**2 * rmaj * dp_new / drhodpsi
         if (writelots) write(*,*) 'alpha = ', -qval**2 * rmaj * dp / drhodpsi,', ', &
              -qval**2 * rmaj * dp_new / drhodpsi
      case default
         write(fort11_unit,*) 'd beta/d rho = ',2 * dp / drhodpsi,', ',2 * dp_new / drhodpsi
         if (writelots) write(*,*) 'd beta/d rho = ',2 * dp / drhodpsi,', ',2 * dp_new / drhodpsi
      end select

      if ((.not. local_eq) .and. (bishop == 1 .or. bishop > 9)) then
         di_new = di
      else
         di_new = (s_hat_new / tmp - 2 * c_b - dp_new * b_b) / a_b
      end if

      write(fort11_unit,*) 's_hat ',s_hat,', ',s_hat_new
      if (writelots) write(*,*) 's_hat ',s_hat,', ',s_hat_new
      write(fort11_unit,*)

      shat = s_hat_new
      dqdrp = s_hat_new * qval * drhodrp / rhoc
      dbetadrho = 2 * dp_new * dpsidrho

      ! Calculate grad(S)
      ! here the coordinates are rho and l
      call gradl(ltheta, seik, dSdl, twopi * qval, nth)
      call periodic_copy(dSdl, 0., nperiod)
      grads(:,1) = a_bish * di_new + b_bish * dp_new + c_bish * 2
      grads(:,2) = dSdl

      ! Calculate grad(Z)
      ! obtain Z(theta)
      call Ztgrid(rgrid, outputs%theta, Zoftheta)
      call periodic_copy(Zoftheta, 0., nperiod)
      ! calculate dZ/dl
      call gradl(ltheta, Zoftheta, dZdl, 0., nth)
      call periodic_copy(dZdl, 0., nperiod)

      ! gradztot(:,1) is grad Z . rhohat
      gradztot(:,1) = crpgrad(:,2) / outputs%grho
      ! gradztot(:,2) is grad Z . lhat
      gradztot(:,2) = dZdl
      ! gradztot(:,3) is grad Z . phihat
      gradztot(:,3) = 0.0

      if (save_bishop) then
         call open_output_file(bishop_unit, '.bishop')
         write(bishop_unit, fmt="('#di, di_new = ',F12.7,' ',F12.7)") di, di_new
         write(bishop_unit, fmt="('#dp, dp_new = ',F12.7,' ',F12.7)") dp, dp_new
         write(bishop_unit, fmt="('#bi =', F12.7)") bi
         write(bishop_unit, fmt="('#drdpsi =', F12.7)") drhodpsi
         write(bishop_unit, fmt="('#a_b (fp)=', F20.10, ' ,b_b (pp)=', F20.10,' ,c_b  =', F20.10,' ,s_hat =', F20.10)") a_b*di, b_b*dp, c_b, s_hat
         write(bishop_unit, fmt ="('#  theta, a_bish,  b_bish,  c_bish, rmajor, Z, Bpol, invRpol, th_bish')")
         do i = -ntgrid, ntgrid
            write(bishop_unit, '( 10(2X, F25.17) )') outputs%theta(i), a_bish(i), b_bish(i), &
                 c_bish(i), rmajor(i), Zoftheta(i),  Bpolmag(i), invRpol(i), th_bish(i)
         end do
         call close_output_file(bishop_unit)
      end if
    end subroutine bishop_grads

    subroutine bishop_large_aspect_ratio_approximation
      real :: delrp, rp1, rp2, qval1, qval2, bi1, bi2, dqdrho, dsdrptot, pbar_1, pbar_2
      real, dimension(:), allocatable :: seik1, seik2, dummy2, dsdrp
      integer :: itot, k
      !     rp derivative
      delrp = delrho / drhodrp
      rp1 = rp - delrp ; rp2 = rp + delrp

      ! Note for local_eq we need to ensure qval is set before
      !calling seikon. We force bishop > 0 for local_eq so shouldn't get here.
      if (local_eq) then
         dqdrho = dqdrhofun(rhoc, shat, fort11_unit)
         qval1 = qval - dqdrho * delrho ; qval2 = qval + dqdrho * delrho
      end if

      pbar_1 = pbarfun(rp1, outputs%theta(0)) ; pbar_2 = pbarfun(rp2, outputs%theta(0))
      ! Calculate modified btori -- note we use the fact that rgrid == rp in
      ! seikon to avoid having to form rgrid1/rgrid2 here.
      bi1 = btori(pbar_1) ; bi2 = btori(pbar_2)

      allocate(seik1(-ntgrid:ntgrid), seik2(-ntgrid:ntgrid), dummy2(-ntgrid:ntgrid))
      allocate(dsdrp(-nth:nth))
      ! Note seikon sets rgrid = rp. This is only the same as the regular
      ! rgrid for equilibrium types other than efit or dfit. For those cases
      ! this _might_ not do what is expected.
      call seikon(rp1, outputs%theta, invR, qval1, bi1, seik1, dummy2, dummy, dum)
      call seikon(rp2, outputs%theta, invR, qval2, bi2, seik2, dummy2, dummy, dum)

      do i = -nth, nth
         dsdrp(i) = (seik2(i) - seik1(i)) / (2 * delrp)
      end do

      if (.not. local_eq) then
         dqdrho = (qval2 - qval1) / (2 * delrho)
         shat = dqdrho * rhoc / qval
      end if
      dqdrp = dqdrho * drhodrp ! Output

      if (writelots) then
         write (fort11_unit,*) 'dqdrho=',dqdrho,' qval1,2= ',qval1,qval2
         write (fort11_unit,*) 's_hat= ',shat
         write (fort11_unit,*) 'i,rgrid,rgrid1,rgrid2:'
         do i = -nth, nth
            write (fort11_unit,*) i, rgrid(i), rp1, rp2
         end do
      end if

      s_hat_new = shat ! Output
      dbetadrho = 2 * dpfun(pbar) / drhodrp ! Output

      !     compute gradient of S
      ! this is really gradient of alpha
      ! note that coordinates are R and Z (not rho and l, as in bishop case)
      ! Note we don't use dsdthet1/2 (dummy2 above)
      do k = -nperiod+1, nperiod-1
         do i = -nth, nth
            itot = i + k * ntheta
            ! By calling periodic_copy(dsdrp, twopi * dqdrp) we should
            ! be able to use dsdrp directly here?
            dsdrptot = dsdrp(i) - k * twopi * dqdrp
            grads(itot, :) = rpgrad(i, :) * dsdrptot + thgrad(i, :) * dsdthet(i)  ! Output
         end do
      end do
      ! gradztot(:,1) is grad Z . Rhat ! Output
      gradztot(:,1) = 0.
      ! gradztot(:,2) is grad Z . Zhat ! Output
      gradztot(:,2) = 1.
      ! gradztot(:,3) is grad Z . phihat ! Output
      gradztot(:,3) = 0.
    end subroutine bishop_large_aspect_ratio_approximation

    !> FIXME : Add documentation
    subroutine alloc_local_arrays(n)
      implicit none
      integer, intent(in) :: n

      allocate(seik   (-n:n), &
           dsdthet    (-n:n), &
           th_bish    (-n:n), &
           invRpol    (-n:n), &
           invR       (-n:n), &
           trip       (-n:n), &
           ltheta     (-n:n), &
           rgrid      (-n:n), &
           dummy      (-n:n), &
           Bpolmag    (-n:n), &
           Bmod       (-n:n), &
           rmajor     (-n:n))

      allocate(thgrad (-n:n,2), &
           rpgrad     (-n:n,2), &
           crpgrad    (-n:n,2), &
           grads      (-n:n,2))

      allocate(gradrptot (-n:n,3), &
           gradstot   (-n:n,3), &
           gradztot   (-n:n,3), &
           gradthtot  (-n:n,3))  ! MAB
    end subroutine alloc_local_arrays
  end subroutine eikcoefs

  !> Allocate the geom instance to the correct type
  pure subroutine allocate_geom(geom)
    use geq, only: geq_type
    use leq, only: leq_type
    use ceq, only: ceq_type
    use peq, only: peq_type, teq_type
    use eeq, only: eeq_type, gs2d_type
    use deq, only: deq_type
    use ideq, only: ideq_type
    use geo_utils, only: abstract_geo_type
    implicit none
    class(abstract_geo_type), allocatable, intent(in out) :: geom
    if (gen_eq) then
       allocate(geq_type :: geom)
    else if (chs_eq) then
       allocate(ceq_type :: geom)
    else if (gs2d_eq) then
       allocate(gs2d_type :: geom)
    else if (efit_eq) then
       allocate(eeq_type :: geom)
    else if (transp_eq) then
       allocate(teq_type :: geom)
    else if (ppl_eq) then
       allocate(peq_type :: geom)
    else if (dfit_eq) then
       allocate(deq_type :: geom)
    else if (idfit_eq) then
       allocate(ideq_type :: geom)
    else if (local_eq) then
       allocate(leq_type :: geom)
    end if
  end subroutine allocate_geom

  pure subroutine calculate_metric_terms(gradstot, gradrptot, gradthtot, dpsidrho, dqdrp, &
       qval, rhoc, bmod, force_sym, gds2, gds21, gds22, gds23, gds24, gds24_noq)
    real, dimension(-ntgrid:, :), intent(in) :: gradstot, gradrptot, gradthtot
    real, intent(in) :: dpsidrho, dqdrp, qval, rhoc
    real, dimension(-ntgrid:), intent(in) :: bmod
    logical, intent(in) :: force_sym
    real, dimension(-ntgrid:), intent(out) :: gds2, gds21, gds22, gds23, gds24, gds24_noq
    real, dimension(:), allocatable :: gdsdum1, gdsdum2, gdsdum3, gdsdum4, gdsdum5
    allocate(gdsdum1(-ntgrid:ntgrid), gdsdum2(-ntgrid:ntgrid), gdsdum3(-ntgrid:ntgrid))
    call dottgridf(gradstot, gradstot, gdsdum1)
    call dottgridf(gradstot, gradrptot, gdsdum2)
    call dottgridf(gradrptot, gradrptot, gdsdum3)
    gds2  = gdsdum1            *dpsidrho**2
    gds21 = gdsdum2*dqdrp      *dpsidrho**2
    gds22 = gdsdum3*dqdrp*dqdrp*dpsidrho**2

    ! v_E . grad theta = i*ky*phi*(rho_{r}*v_{t,r}/a**2)*(gds23 + gsd24*theta0)
    ! calculation of gds23 and gds24 assumes rp=rho, which is not necessarily true
    ! for numerical equilibria
    allocate(gdsdum4(-ntgrid:ntgrid), gdsdum5(-ntgrid:ntgrid))
    call dottgridf(gradrptot, gradthtot, gdsdum4) !grad rho dot grad theta
    call dottgridf(gradstot, gradthtot, gdsdum5) !grad alpha dot grad theta
    gds23 = (gdsdum4 * gds2 - gdsdum5 * gdsdum2 * dpsidrho**2) / bmod**2
    gds24_noq = (gdsdum4 * gdsdum2 - gdsdum5 * gdsdum3) * (dpsidrho / bmod)**2
    gds24 = gds24_noq * dqdrp
    gds24_noq = gds24_noq * qval / rhoc

    if (force_sym) then
       call sym(gds2, 0, ntgrid)
       call sym(gds21, 1, ntgrid)
       call sym(gds22, 0, ntgrid)
       ! not sure about the following yet -- MAB
       ! call sym(gds23, 0, ntgrid); call sym(gds24, 1, ntgrid)
    end if
  end subroutine calculate_metric_terms

  subroutine calculate_surface_area(theta, gradpar, bmag, grho, dpsidrho, rhoc, nth, &
       surfarea, dvdrhon)
    use constants, only: twopi
    use file_utils, only: open_output_file, close_output_file
    real, dimension(-ntgrid:), intent(in) :: theta, gradpar, bmag, grho
    real, intent(in) :: dpsidrho, rhoc
    integer, intent(in) :: nth
    real, intent(out) :: surfarea, dvdrhon
    real, dimension(:), allocatable :: jacob, ds, dummy
    integer :: i, fort25_unit
    ! compute grad rho*Jacobian of the transformation to the new coordinates:
    ! Axisymmetry assumed to get this (nu indep. of phi):
    allocate(jacob(-ntgrid:ntgrid), ds(-ntgrid:ntgrid), dummy(-ntgrid:ntgrid))
    do i = -nth, nth
       jacob(i) = abs(dpsidrho / (gradpar(i)*bmag(i)))
       ds(i) = grho(i) * jacob(i)
    end do

    ! compute surface area from A=Int(R sqrt(r**2+drdtheta**2) dtheta dphi)
    !  does not work with equal_arc grid
    ! compute the surface area from A=Int(J |grad rho| dtheta dalpha)
    call integrate(ds, theta, dummy, nth)
    surfarea = twopi * (dummy(nth) - dummy(-nth))

    ! compute dV/drhon = 2 pi Int(J dtheta)
    call integrate(jacob, theta, dummy, nth)
    dvdrhon = twopi * (dummy(nth) - dummy(-nth))

    call open_output_file(fort25_unit,".fort.25")
    write(fort25_unit,*) rhoc, f_trap(bmag(-nth:nth), theta(-nth:nth), jacob(-nth:nth))
    call close_output_file(fort25_unit)
  end subroutine calculate_surface_area

  !> FIXME : Add documentation
  pure subroutine bvectortgrid(invR, nth, gradrptot, dpsidrp, bi, bvector)
    integer, intent (in) :: nth
    real, dimension(-ntgrid:  ), intent (in) :: invR
    real, dimension(-ntgrid:,:), intent (in) :: gradrptot
    real, dimension(-ntgrid:,:), intent (out) :: bvector
    real, intent (in) :: dpsidrp, bi
    integer :: i
    do i = -nth, nth
       bvector(i,1) = -dpsidrp * gradrptot(i,2) * invR(i)
       bvector(i,2) = dpsidrp * gradrptot(i,1) * invR(i)
       bvector(i,3) = bi * invR(i)
    end do
  end subroutine bvectortgrid

  !> FIXME : Add documentation
  pure subroutine gradstottgrid (rmajor, grads, gradstot)
    real, dimension(-ntgrid:), intent (in) :: rmajor
    real, dimension(-ntgrid:, :), intent (in) :: grads
    real, dimension(-ntgrid:, :), intent (out) :: gradstot
    integer :: i
    gradstot(:, 1:2) = grads
    if (use_large_aspect) then
       gradstot(:, 3) = 0.
    else
       do i = -nth, nth
          gradstot(i, 3) = 1.0 / rmajor(i)
       end do
       call periodic_copy(gradstot(:, 3), 0.0, nperiod)
    end if
  end subroutine gradstottgrid

  !> Calculate the cross product c = a x b on the `[-nth, nth]` domain
  pure subroutine crosstgrid(a, b, c)
    real, dimension(-ntgrid:, :), intent (in) ::  a, b
    real, dimension(-ntgrid:, :), intent (out) ::  c
    integer :: i
    do i = -nth, nth
       c(i, 1) = a(i, 2) * b(i, 3) - b(i, 2) * a(i, 3)
       c(i, 2) = a(i, 3) * b(i, 1) - b(i, 3) * a(i, 1)
       c(i, 3) = a(i, 1) * b(i, 2) - b(i, 1) * a(i, 2)
    end do
  end subroutine crosstgrid

  !> Calculate the dot product c = a.b where a is on the `[-nth, nth]` domain
  !> and b and c are both on the `[-ntgrid, ntgrid]` domain. In other words a
  !> is considered periodic on 2pi whilst b may not be.
  pure subroutine dottgrid(a, b, c)
    real, dimension (-ntgrid:, :), intent(in) :: a, b
    real, dimension (-ntgrid:), intent(out) :: c
    integer :: i, k, itot

    do k = -nperiod+1, nperiod-1
       do i = -nth, nth
          itot = i + k * ntheta
          c(itot) = a(i, 1) * b(itot, 1) + a(i, 2) * b(itot, 2) + a(i, 3) * b(itot, 3)
       enddo
    enddo
  end subroutine dottgrid

  !> Calculate the dot product c = a.b where a, b and c are on the `[-ntgrid, ntgrid]` domain
  pure subroutine dottgridf(a, b, c)
    real, dimension(-ntgrid:,:), intent(in) :: a, b
    real, dimension(-ntgrid:), intent(out) :: c
    c = a(:, 1) * b(:, 1) + a(:, 2) * b(:, 2) + a(:, 3) * b(:, 3)
  end subroutine dottgridf

  !> Returns the major radius R corresponding to the passed r and theta
  real function Rpos(r, thet)
    implicit none
    real, intent (in) :: r, thet
    Rpos = geom%Rpos(r, thet)
  end function Rpos

  !> Returns the height Z corresponding to the passed r and theta
  real function Zpos(r, thet)
    implicit none
    real, intent (in) :: r, thet
    Zpos = geom%Zpos(r, thet)
  end function Zpos

  !> Returns the inverse of the major radius, 1/R, corresponding to the passed r and theta
  !> Should we implement this in terms of 1/Rpos?
  real function invRfun(r, thet)
    implicit none
    real, intent (in) :: r, thet
    if(use_large_aspect) then
       invRfun = 1. / rmaj
    else
       invRfun = geom%invR(r, thet)
    end if
  end function invRfun

  !> Returns the diameter of the surface labelled by rp.
  real function diameter(rp)
    implicit none
    real, intent(in) :: rp
    integer :: i
    !Do we really want to silently ignore requests out of range? Perhaps we should remove.
    if (rp <= rpmin .and. .not. (efit_eq .or. gs2d_eq)) then
       diameter = 0.
    else
       diameter = geom%diameter(rp)
    end if
  end function diameter

  !> FIXME : Add documentation
  real function btori(pbar)
    implicit none
    real, intent (in) :: pbar
    btori = geom%btori(pbar)
  end function btori

  !> FIXME : Add documentation
  !> returns dI/dpsi
  real function dbtori(pbar)
    implicit none
    real, intent (in) :: pbar
    dbtori = geom%dbtori(pbar)
  end function dbtori

  !> FIXME : Add documentation
  !> Only used with ~undocumented `irho = 4` option. Ideally would update to take pbar
  real function rhofun (rp)
    implicit none
    real, intent (in) :: rp
    real :: pbar
    pbar = psifun(rp)
    rhofun = geom%rhofun(pbar)
  end function rhofun

  !> FIXME : Add documentation
  real function qfun(pbar)
    implicit none
    real, intent(in) :: pbar
    qfun = geom%qfun(pbar)
  end function qfun

  !> FIXME : Add documentation
  real function pfun(pbar)
    implicit none
    real, intent (in) :: pbar
    pfun = geom%pfun(pbar)
  end function pfun

  !> FIXME : Add documentation
  real function dpfun(pbar)
    implicit none
    real, intent (in) :: pbar
    dpfun = geom%dpfun(pbar)
  end function dpfun

  !> FIXME : Add documentation
  real function beta_a_fun(pbar)
    implicit none
    real, intent (in) :: pbar
    beta_a_fun = geom%betafun(pbar)
  end function beta_a_fun

  real function rfun(rp, theta)
    implicit none
    real, intent(in) :: rp, theta
    rfun = geom%rfun(rp, theta)
  end function rfun

  !> FIXME : Add documentation
  real function psi(r, thet)
    implicit none
    real, intent (in) :: r, thet
    psi = geom%psi(r, thet)
  end function psi

  !> FIXME : Add documentation
  real function rcenter(rp)
    implicit none
    real, intent (in) :: rp
    rcenter = geom%rcenter(rp)
  end function rcenter

  !> Calculates the normalised radial coordinate
  pure real function psifun (rp)
    implicit none
    real, intent (in) :: rp
    psifun = min(1.0, max(0.0, (rp - psi_0) / (psi_a - psi_0)))
  end function psifun

  real function pbarfun(r, theta)
    implicit none
    real, intent(in) :: r, theta
    pbarfun = psifun(psi(r, theta))
  end function pbarfun

  !> FIXME : Add documentation
  subroutine drift(rgrid, rp, gradstot, gradrptot, gradztot, dqdrp, dpsidrho, &
       dpsidrp, Bmod, Bpolmag, invRpol, invR, th_bish, ltheta, theta, rmajor, dp_new, &
       bi, force_sym, gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0)
    implicit none
    real, dimension (-ntgrid:), intent (in) :: rgrid, Bmod, Bpolmag, invRpol, invR, th_bish
    real, dimension (-ntgrid:), intent (in) :: ltheta, theta, rmajor
    real, dimension (-ntgrid:, :), intent (in) :: gradstot, gradrptot, gradztot
    real, intent (in) :: rp, dqdrp, dpsidrho, dpsidrp, dp_new, bi
    logical, intent(in) :: force_sym
    real, dimension(-ntgrid:ntgrid), intent(out) :: gbdrift, gbdrift0, cvdrift, cvdrift0, cdrift, cdrift0
    real, dimension (-ntgrid:ntgrid, 3) :: bgradtot, pgradtot, dummy, curve, bvector
    real, dimension (-ntgrid:ntgrid, 2) :: pgrad, bgrad1
    real, dimension (-ntgrid:ntgrid) :: gbdrift1, gbdrift2, cvdrift1, cvdrift2, cdrift1, cdrift2
    character(1) :: char

    if (use_large_aspect) then
       gbdrift = (-sin(theta) * gradstot(:,1) - cos(theta) * gradstot(:,2)) / rmaj
       gbdrift = 2 * dpsidrho * gbdrift

       gbdrift0 = (-sin(theta) * gradrptot(:,1) - cos(theta) * gradrptot(:,2)) / rmaj
       gbdrift0 = 2 * dpsidrho * gbdrift0 * dqdrp

       cvdrift = gbdrift ; cvdrift0 = gbdrift0
    else
       !     compute magnetic field
       call bvectortgrid(invR, nth, gradrptot, dpsidrp, bi, bvector)

       if (bishop >= 1) then
          call bishop_gradB(Bmod, Bpolmag, invRpol, th_bish, ltheta, rmajor, dp_new, bi, bgrad1)
       else
          char = 'B'
          call geom%gradient(rgrid, theta, bgrad1, char, rp, nth, ntgrid)
       end if

       bgradtot(-nth:nth, 1:2) = bgrad1(-nth:nth, :)
       bgradtot(-nth:nth, 3) = 0.0

       call crosstgrid(bvector, bgradtot, dummy)
       call dottgrid(dummy, gradstot, gbdrift1)
       call dottgrid(dummy, gradrptot, gbdrift2)

       !     create curvature
       if (bishop == 0) then
          char='R'
          call geom%gradient(rgrid, theta, pgrad, char, rp, nth, ntgrid)
          pgradtot(-nth:nth, 1:2) = pgrad(-nth:nth, :)
          pgradtot(-nth:nth, 3) = 0.0
       else
          pgradtot = gradrptot * dpsidrp * dp_new
       end if

       curve(-nth:nth, 1) = bgradtot(-nth:nth, 1) / bmod(-nth:nth) + &
            pgradtot(-nth:nth, 1) / bmod(-nth:nth)**2
       curve(-nth:nth, 2) = bgradtot(-nth:nth, 2) / bmod(-nth:nth) + &
            pgradtot(-nth:nth, 2) / bmod(-nth:nth)**2
       curve(-nth:nth, 3) = 0.

       call crosstgrid(bvector, curve, dummy)
       call dottgrid(dummy, gradstot, cvdrift1)
       call dottgrid(dummy, gradrptot, cvdrift2)

       ! factor of 2 appears below because will later multiply
       ! by wunits, which has a factor of 1/2 built-in
       gbdrift = 2 * dpsidrho * gbdrift1 / bmod**3
       gbdrift0 = 2 * dpsidrho * gbdrift2 * dqdrp / bmod**3
       cvdrift = 2 * dpsidrho * cvdrift1 / bmod**2
       cvdrift0 = 2 * dpsidrho * cvdrift2 * dqdrp / bmod**2
    end if

    ! coriolis drift term is -2*vpa*om_phi/Omega * (bhat x (zhat x bhat)) . grad (g+q<phi>F_0/T)
    ! here we calculate 2*(zhat . grad (alpha+q*theta0))/(B/Bref)

    ! this is the Zhat . grad(alpha) part of Zhat . grad(S)
    ! note that a factor of k_y = (n0/a) (drho/dpsiN) has been factored out
    call dottgrid(gradztot, gradstot, cdrift1)
    ! this is the Zhat . grad(q) part of Zhat . grad(S)
    ! note that a factor of k_y * theta0 has been factored out
    call dottgrid(gradztot, gradrptot, cdrift2)

    ! factor of 4 (instead of 2) below because will later mutliply
    ! by wunits, which has a built-in factor of 1/2
    ! cdrift0 is part of coriolis drift dependent on theta0
    cdrift  = -4 * dpsidrho * cdrift1 / bmod
    cdrift0 = -4 * dpsidrho * cdrift2 * dqdrp / bmod

    if (force_sym) then
       call sym(gbdrift, 0, ntgrid)
       call sym(cvdrift, 0, ntgrid)
       call sym(gbdrift0, 1, ntgrid)
       call sym(cvdrift0, 1, ntgrid)
       ! should maybe add in cdrift and cdrift0 here -- MAB
    end if
  end subroutine drift

  !> Forces symmetry in theta in the passed array. If isign == 0 then
  !> `a(+/- theta)` is replaced by the average of `a(+/-theta)`, otherwise
  !> `a(+/- theta)` is replaced by `+/- (a(theta)-a(-theta))`.
  pure subroutine sym(a, isign, ntgrid)
    implicit none
    integer, intent(in) :: isign, ntgrid
    real, dimension(-ntgrid:), intent (in out) :: a
    integer :: i
    if (isign == 0) then
       do i = 1, ntgrid
          a(i) = 0.5 * (a(i) + a(-i))
          a(-i) = a(i)
       end do
    else
       do i = 1, ntgrid
          a(i) = 0.5 * (a(i) - a(-i))
          a(i) = -a(-i)
       end do
       a(0) = 0.
    end if
  end subroutine sym

  !> FIXME: Add documentation. Only used for bishop == 0 testing
  subroutine seikon(rp, theta, invR, qval, bi, seik, dsdthet, rgrid, dpsidrp)
    real, intent(in) :: rp
    real, dimension(-ntgrid:), intent(in) :: theta, invR
    real, intent(in out) :: qval
    real, intent(in) :: bi
    real, dimension(-ntgrid:), intent(out) :: seik, dsdthet, rgrid
    real, intent(out) :: dpsidrp
    real, dimension(-ntgrid:ntgrid, 2) :: thgrad, rpgrad
    real, dimension(-ntgrid:ntgrid) :: trip
    character(1), parameter :: char = 'P'
    rgrid = rp
    ! Note we only call this routine when bishop == 0
    if (bishop == 0)  then
       call geom%gradient(rgrid, theta, rpgrad, char, rp, nth, ntgrid)
    else
       call geom%bgradient(rgrid, theta, rpgrad, char, rp, nth, ntgrid)
    end if
    call thetagrad(rgrid, theta, thgrad)
    call tripprod2dtgrid(rpgrad, thgrad, invR, trip)
    call eikonal(theta, invR, trip, qval, bi, seik, dsdthet, dpsidrp)
  end subroutine seikon

  !> FIXME : Add documentation
  pure subroutine eikonal(theta, invR, trip, qval, bi, seik, dsdthet, dpsidrp)
    use constants, only: twopi
    implicit none
    real, dimension(-ntgrid:), intent(in) :: theta, invR, trip
    real, dimension(-ntgrid:), intent(out) :: seik, dsdthet
    real, intent(in out) :: qval
    real, intent(in) :: bi
    real, intent(out) :: dpsidrp
    integer :: i

    do i = -nth, nth
       dsdthet(i) = -bi * invR(i)**2 / trip(i)
    end do
    call periodic_copy(dsdthet, 0., nperiod)

    call integrate(dsdthet, theta, seik, ntgrid)

    if (local_eq) then
       dpsidrp = -(seik(nth) - seik(-nth)) / (twopi * qval)
       seik = seik / dpsidrp
       dsdthet = dsdthet / dpsidrp
    else
       dpsidrp = 1.
       qval = -(seik(nth) - seik(-nth)) / (twopi)
    end if
  end subroutine eikonal

  !> FIXME : Add documentation
  subroutine thetagrad(rgrid, theta, thgrad)
    implicit none
    real, dimension(-ntgrid:), intent (in) :: rgrid, theta
    real, dimension(-ntgrid:, :), intent (out) :: thgrad
    real :: dum
    character(1), parameter :: char = 'T'
    if (bishop == 0) then
       call geom%gradient(rgrid, theta, thgrad, char, dum, nth, ntgrid)
    else
       call geom%bgradient(rgrid, theta, thgrad, char, dum, nth, ntgrid)
    end if
    call periodic_copy(thgrad(:, 1), 0., nperiod)
    call periodic_copy(thgrad(:, 2), 0., nperiod)
  end subroutine thetagrad

  !> FIXME : Add documentation
  real function rpofrho(rho, fort11_unit)
    implicit none
    real, intent(in) :: rho
    integer, intent(in) :: fort11_unit
    real, parameter :: xerrbi = 1.e-4, xerrsec = 1.e-8
    real :: a, b, fval
    integer :: ier

    ier = 123456 ! Default to indicating error
    a = rpmin ; b = rpmax

    if (irho == 1) then
       fval = rho * rho * phi(b)
       call root(phi, fval, a, b, xerrbi, xerrsec, fort11_unit, ier, rpofrho)
    else if (irho == 2) then
       fval = rho * diameter(b)
       call root(diameter, fval, a, b, xerrbi, xerrsec, fort11_unit, ier, rpofrho)
    else if (irho == 3) then
       fval = rho * psifun(b)
       call root(psifun, fval, a, b, xerrbi, xerrsec, fort11_unit, ier, rpofrho)
    else if (irho == 4) then
       fval = rho * rhofun(b)
       call root(rhofun, fval, a, b, xerrbi, xerrsec, fort11_unit, ier, rpofrho)
    endif

    if (debug) then
       write (*,*) "Values in root: "
       write (*,*) 'rho:  ', rho
       write (*,*) 'fval: ', fval
       write (*,*) 'rpofrho: ', rpofrho
       write (*,*) 'rpmin: ', rpmin
       write (*,*) 'rpmax: ', rpmax
    end if
    if (ier > 0) write(fort11_unit,*) 'error in rpofrho, rho=', rho, rpofrho, ier
  end function rpofrho

  !> Computes something like x cross y dot grad zeta.
  pure subroutine tripprod2dtgrid(x, y, invR, val)
    implicit none
    real, dimension(-ntgrid:), intent (in) :: invR
    real, dimension(-ntgrid:), intent (out) :: val
    real, dimension(-ntgrid:, :), intent (in) :: x, y
    integer :: i
    ! factor of 1/R from grad zeta
    do i = -nth, nth
       val(i) = (x(i, 1) * y(i, 2) - x(i, 2) * y(i, 1)) * invR(i)
    end do
  end subroutine tripprod2dtgrid

  !> Integrates arg(grid) from i = 0 (generally theta = 0) to either
  !> ends of the grid (usually +/- pi) storing the cumulative integral
  !> in the output `ans`.
  pure subroutine integrate(arg, grid, ans, n)
    real, dimension(-ntgrid:), intent (in) :: arg, grid
    real, dimension(-ntgrid:), intent (out) :: ans
    integer, intent(in) :: n
    integer :: i
    ans = 0.
    do i = 1, n
       ans(i) = ans(i-1) + 0.5 * (grid(i) - grid(i-1)) * (arg(i) + arg(i-1))
    end do
    do i = -1, -n, -1
       ans(i) = ans(i+1) + 0.5 * (grid(i) - grid(i+1)) * (arg(i) + arg(i+1))
    end do
  end subroutine integrate

  !> Calculates rmajor = R(theta), where R is the distance
  !> to the central axis (rmajor not to be confused with rmaj, r_geo etc)
  !> rgrid is the submodule radial grid (r_circ for EFIT & deq, rp otherwise)
  subroutine rmajortgrid(rgrid, theta, rmajor)
    implicit none
    real, dimension(-ntgrid:), intent(in) :: rgrid, theta
    real, dimension(-ntgrid:), intent(out) :: rmajor
    integer :: i
    do i = -nth, nth
       rmajor(i) = Rpos(rgrid(i), theta(i))
    end do
  end subroutine rmajortgrid

  !> Calculates invR = 1/R(theta), where R is the distance
  !> to the central axis (rmajor not to be confused with rmaj, r_geo etc)
  !> rgrid is the submodule radial grid (r_circ for EFIT & deq, rp otherwise)
  subroutine invRtgrid(rgrid, theta, invR)
    implicit none
    real, dimension(-ntgrid:), intent(in) :: rgrid, theta
    real, dimension(-ntgrid:), intent(out) :: invR
    integer :: i
    do i = -nth, nth
       invR(i) = invRfun(rgrid(i), theta(i))
    end do
  end subroutine invRtgrid

  !> Calculates the vertical height, Z, on the theta grid
  subroutine Ztgrid(rgrid, theta, Zoftheta)
    implicit none
    real, dimension(-ntgrid:), intent(in) :: rgrid, theta
    real, dimension(-ntgrid:), intent(out) :: Zoftheta
    integer :: i
    do i = -nth, nth
       Zoftheta(i) = Zpos(rgrid(i), theta(i))
    end do
  end subroutine Ztgrid

  !> FIXME : Add documentation
  real function pbarofrho(rho, fort11_unit)
    implicit none
    real, intent (in) :: rho
    integer, intent(in) :: fort11_unit
    pbarofrho = psifun(rpofrho(rho, fort11_unit))
  end function pbarofrho

  !> FIXME : Add documentation
  real function dqdrhofun(rho, shat, fort11_unit)
    implicit none
    real, intent (in) :: rho, shat
    integer, intent(in) :: fort11_unit
    dqdrhofun = shat * qfun(pbarofrho(rho, fort11_unit)) / rho
  end function dqdrhofun

  !> Find soln where f(soln) ~= fval. Works best if a and b bracket the root.
  !> We can then first bisect the bracket to narrow in on the solution.
  !> Following this we use the secant method to find the solution.
  !> The input xerrbi influences the number of bisection iterations whilst
  !> xerrsec is used to set the stopping tolerance on the secant method.
  subroutine root(f, fval, a, b, xerrbi, xerrsec, fort11_unit, ier, soln)
    implicit none
    interface
      real function f(x_)
        real, intent(in) :: x_
      end function f
    end interface
    real, intent(in) :: fval, a, b, xerrbi, xerrsec
    real, intent(out) :: soln
    integer, intent(in) :: fort11_unit
    integer, intent(out) :: ier
    real :: a1, b1, f1, f2, f3, trial, aold
    integer :: i, niter
    logical :: skip_bisection
    ier = 0 ; a1 = a ; b1 = b
    f1 = f(a1) - fval ;  f2 = f(b1) - fval
    skip_bisection = xerrbi <= 0 !Could be xerrbi <= abs(b1 - a1)?
    if (.not. skip_bisection .and. f1 * f2 > 0.) then
       write(fort11_unit, *) 'f1 and f2 have same sign in bisection routine root'
       write(fort11_unit, *) 'a1,f1,b1,f2=', a1, f1, b1, f2
       write(fort11_unit, *) 'fval=', fval
       ier = 1
       skip_bisection = .true.
    end if

    if (.not. skip_bisection) then
       ! Presumably we require |b1 - a1| >= xerrbi here to avoid niter < 1?
       niter = 1 + int(log(abs(b1 - a1) / xerrbi) / log(2.))
       ! Repeatedly bisect between limits, swapping one limit
       ! for the bisection point each iteration. Helps to narrow
       ! down the bracketing region substantially.
       do i = 1, niter
          trial = 0.5 * (a1 + b1)
          f3 = f(trial) - fval
          if (f3 * f1 > 0.) then
             a1 = trial
             f1 = f3
          else
             b1 = trial
             f2 = f3
          end if
       end do
    end if

    !Swap so |f1| < |f2|
    if( abs(f1) > abs(f2) ) then
       f3 = f1 ; f1 = f2 ; f2 = f3
       aold = a1 ; a1 = b1 ; b1 = aold
    end if
    !     start secant method
    do i = 1, 10
       aold = a1 ; f3 = f2 - f1
       if (abs(f3) < 1.e-11)  f3 = 1.e-11 !Should we exit here? Saying f2~f1
       a1 = a1 - f1 * (b1 - a1) / f3 !Update guess for root.
       f2 = f1 ; b1 = aold ; f1 = f(a1) - fval
       if (abs(a1 - b1) < xerrsec) exit !Root found closely enough so exit
    end do
    soln = a1
    ! Detect instances where we have not converged to a solution
    if (abs(f1) > xerrsec) ier = ier + 2
  end subroutine root

  !> FIXME : Add documentation
  real function phi(rp)
    implicit none
    real, intent (in) :: rp
    integer, parameter :: nimax = 200
    real :: pbar, pb(nimax), dpb, qi, qnext
    integer :: i

    pbar = psifun(rp)
    ! Construct a uniformly space grid from 0 to pbar
    dpb = pbar / float(nimax - 1)
    do i = 1, nimax
       pb(i) = dpb * float(i-1)
    end do

    ! Integrate area under q(pb) curve from pb = 0 to pbar
    phi = 0. ; qi = qfun(pb(1))
    do i = 1, nimax - 1
       qnext = qfun(pb(i+1))
       phi = phi + 0.5 * dpb * (qnext + qi)
       qi = qnext
    end do
    phi = phi * (psi_a - psi_0)
  end function phi

  !> FIXME : Add documentation
  !>
  !> @note This routine is just used by eik, not GS2
  subroutine geofax(rho, theta, t, e, d, fort11_unit)
    implicit none
    real, intent (in) :: rho
    real, dimension(-ntgrid:ntgrid), intent(in) :: theta
    real, intent (out) :: t, e, d
    integer, intent(in) :: fort11_unit
    real, dimension(-ntgrid:ntgrid) :: rgrid
    real :: a, b, aa, bb, rp, hmaxu, hmaxl, rmaxu, rmaxl, rmid, Rlo, Rhi, Rinnr, Routr
    integer :: i
    logical, save :: first=.true.
    if (first) then
       write(*,'("WARNING: GEOFAX USES UNINITIALISED VALUES -- RESULT MAY NOT BE SENSIBLE")')
       first = .false.
    end if
    rp = rpofrho(rho, fort11_unit)

    if (efit_eq .or. dfit_eq .or. gs2d_eq) then
       call rtg(rgrid, rp, theta)
    else
       rgrid = rp
    end if

    ! do not treat up-down asymmetry completely for now
    hmaxu = -1.e6 ; hmaxl = 1.e6
    rmaxu = 0. ; rmaxl = 0.

    do i = 0, nth
       a = Rpos(rgrid(i), theta(i))
       b = Zpos(rgrid(i), theta(i))

       aa = Rpos(rgrid(-i), theta(-i))
       bb = Zpos(rgrid(-i), theta(-i))

       if (hmaxu < b ) then
          hmaxu = b ; rmaxu = a
       end if

       if (hmaxl > bb ) then
          hmaxl = bb ; rmaxl = aa
       end if

       Rlo = min(a, aa) ; Rhi = max(a, aa)
       !<DD>WARNING : Rinnr has not been defined on first loop through
       if (Rinnr > Rlo) then
          Rinnr = Rlo
       end if
       !<DD>WARNING : Routr has not been defined on first loop through
       if (Routr < Rhi) then
          Routr = Rhi
       end if
    end do
    a = (Rpos(rgrid(0), theta(0)) - Rpos(rgrid(nth), theta(nth)))
    e = (hmaxu - hmaxl) / a
    ! Ro=(Routr+Rinnr)/2. should = rcenter(rp)
    rmid=(Routr - Rinnr) / 2. ! midplane minor radius of this surface
    ! Using only the upper triangularity:
    ! wasn't properly normalized before 29 Feb 08, and don't use asin:
    !    t = asin(rcenter(rp)-rmax)
    t = (rcenter(rp) - rmaxu) / rmid
    d = rcenter(rp) - rcenter(psi_a)
  end subroutine geofax

  !> Given the input theta and gradpar (= dtheta/dl) calculate and return a theta grid
  !> such that gradpar is constant and theta_eqarc is a theta grid with uniform arc length.
  !> @note It is not clear that the gradpar_eqarc and theta_eqarc are actually consistent
  !> with one another?
  subroutine get_theta_eqarc(ntheta, nperiod, theta, gradpar, gradpar_eqarc, theta_eqarc)
    use splines, only: inter_cspl
    implicit none
    integer, intent (in) :: ntheta, nperiod
    real, dimension(-ntgrid:), intent (in) :: theta, gradpar
    real, dimension(-ntgrid:), intent (out) :: gradpar_eqarc, theta_eqarc
    real, dimension(:), allocatable :: theta_prime, theta_prime_eqarc

    ! In the GK equation we have dg/dl. Our parallel coordinate is
    ! theta, so the term is dtheta/dl dg/dtheta and gradpar =
    ! dtheta/dl.  We now define a new coordinate theta' such that
    ! gradpar' = dtheta'/dl is constant, and the term is now
    ! dtheta'/dl dg/dtheta'. This is the coordinate that is used if
    ! equal_arc = .true. and as a consequence the values of theta that
    ! are used, e.g. in calculating the matrices for invert_rhs is
    ! emphatically not the poloidal angle but this new coordinate
    ! theta'.  It is important to note that when you use equal_arc it
    ! means that gradpar is constant, NOT that there are equal lengths
    ! of arc between gridpoints. In fact the grid points do not move,
    ! it's just that we have changed the definition of our parallel
    ! coordinate, and consequently gradpar. In other words Delta l/
    ! Delta theta' is constant, but Delta l is not... it is the same
    ! as it was before.  In the section of code below, which is used
    ! primarly by GRYFX, we are doing something more. We calculate a
    ! set of theta values such that Delta theta' is constant (rather
    ! than, as is the case for most geometries, Delta theta). As a
    ! consequence, on this new grid, which we call theta_eqarc, Delta
    ! l is itself a constant as a function of grid index.
    gradpar_eqarc = gradpar
    theta_prime = theta
    call arclength (ntheta, nperiod, theta, gradpar_eqarc, theta_prime)
    ! theta'(i) is theta'(i) = theta'(theta(i)) such that
    ! dl / dtheta' is constant. At the moment Delta theta is
    ! constant and Delta theta' is not. We now calculate a new
    ! set of theta values, theta_eqarc, such that Delta theta'
    ! is constant but Delta theta is not. Since we have convenient
    ! equally spaced set of grid points in theta, we use that for
    ! the x values of the interpolation routine.
    ! The four important arguments of the call below
    ! are:
    !     2. An unequally spaced set of theta' values: theta_prime
    !     3. The values of theta that correspond to them: theta
    !     5. A set of equally spaced theta' values: theta_prime_eqarc
    !     6. The new set of unequally spaced theta values, theta_eqarc, which correspond
    !          to equally spaced theta' values and consequently constant
    !          values of Delta l, i.e. constant arc lengths between gridpoints.
    !
    theta_prime_eqarc = theta
    call inter_cspl(size(theta), theta_prime,       theta, &
                    size(theta), theta_prime_eqarc, theta_eqarc)
  end subroutine get_theta_eqarc

  !> Given an initial theta and gradpar (=dtheta/dl) compute a new
  !> theta grid, returned in arcl, which has gradpar = dg/dl  = constant.
  !> @note The input gradpar is modified on output.
  pure subroutine arclength (ntheta, nperiod, theta, gpar, arcl)
    use constants, only: twopi, pi
    implicit none
    integer, intent (in) :: ntheta, nperiod
    real, dimension(-ntgrid:), intent (out) :: arcl
    real, dimension(-ntgrid:), intent (in) :: theta
    real, dimension(-ntgrid:), intent (in out) :: gpar
    integer :: nth, j
    real :: arcfac
    nth = ntheta / 2
    arcl(-nth) = 0.
    do j = -nth, nth-1
       arcl(j+1) = arcl(j) + 0.5 * (theta(j+1) - theta(j)) * (1. / gpar(j) + 1. / gpar(j+1))
    end do
    arcfac = 1. / arcl(nth)

    do j = -nth, nth
       arcl(j) = arcl(j) * twopi * arcfac - pi
    end do

    call periodic_copy(arcl, twopi, nperiod)
    gpar = twopi * arcfac
  end subroutine arclength

  !> FIXME : Add documentation
  subroutine loftheta(rgrid, theta, ltheta)
    implicit none
    real, dimension (-ntgrid:), intent(in) :: rgrid, theta
    real, dimension (-ntgrid:), intent(out) :: ltheta
    integer :: i
    real :: rold, told, rpos_old, zpos_old, rpos_new, zpos_new
    real :: r, thet, rpos_orig, zpos_orig
    ltheta(0) = 0.
    rold = rgrid(0) ; told = theta(0)
    rpos_orig = Rpos(rold, told) ; zpos_orig = Zpos(rold, told)
    rpos_new = rpos_orig ; zpos_new = zpos_orig
    do i = 1, nth
       rpos_old = rpos_new ; zpos_old = zpos_new
       r = rgrid(i) ; thet = theta(i)
       rpos_new = Rpos(r, thet) ; zpos_new = Zpos(r, thet)
       ltheta(i) = ltheta(i-1) + sqrt((Rpos_new - Rpos_old)**2 + (Zpos_new - Zpos_old)**2)
    end do

    rpos_new = rpos_orig ; zpos_new = zpos_orig
    do i = -1, -nth, -1
       rpos_old = rpos_new ; zpos_old = zpos_new
       r = rgrid(i) ; thet = theta(i)
       rpos_new = Rpos(r, thet) ; zpos_new = Zpos(r, thet)
       ltheta(i) = ltheta(i+1) - sqrt((Rpos_new - Rpos_old)**2 + (Zpos_new - Zpos_old)**2)
    end do
  end subroutine loftheta

  !> Compute df / dltheta using a periodic central difference. The
  !> input `ext` can be used to control the effective jump in f
  !> between both ends of the "periodic" grid.
  pure subroutine gradl (ltheta, f, dfdl, ext, n)
    implicit none
    real, dimension (-ntgrid:), intent(in) :: ltheta, f
    real, dimension (-ntgrid:), intent(out) :: dfdl
    real, intent(in) :: ext
    integer, intent(in) :: n
    integer :: i
    ! Note here we're effectively assuming a uniformly (or symetric)
    ! spaced ltheta grid - we could probably use
    !(ltheta(-n) - ltheta(-n+1) + ltheta(n) - ltheta(n-1)) to avoid this
    dfdl(-n) = (f(n-1) + ext - f(-n+1)) / (2. * (ltheta(-n) - ltheta(-n+1)))
    do i = -n+1, n-1
       dfdl(i) = (f(i+1) - f(i-1)) / (ltheta(i+1) - ltheta(i-1))
    end do
    dfdl(n) = dfdl(-n)
  end subroutine gradl

  !> Finds the angle between the surface tangent vector in the poloidal plane
  !> and the direction of increasing major radius.
  pure subroutine th_bishop(rpgrad, th_bish, nth)
    use constants, only: twopi
    implicit none
    integer, intent (in) :: nth
    real, dimension (-ntgrid:, :), intent (in) :: rpgrad
    real, dimension (-ntgrid:), intent (out) :: th_bish
    real, dimension (-ntgrid:ntgrid) :: magrp
    real, dimension (-ntgrid:ntgrid, 2) :: tvec
    integer :: i

    do i = -nth, nth
       magrp(i) = sqrt(rpgrad(i, 1)**2 + rpgrad(i, 2)**2)
    end do

    ! tvec is the tangent vector in the surface of the poloidal plane
    do i = -nth, nth
       tvec(i, 1) = rpgrad(i, 2) / magrp(i)
       tvec(i, 2) = -rpgrad(i, 1) / magrp(i)
    end do

    ! need to find grad(R) and dot it with tvec.
    ! but grad(R) is simple -- it is just (1,0)
    ! so simply return acos(tvec(:,1))
    ! with corrections for left-hand quadrants
    !
    do i = -nth, nth
       tvec(i, 1) = min(1., max(-1., tvec(i, 1)))
       if (tvec(i, 2) > 0.) then
          th_bish(i) = twopi - acos(max(min(tvec(i, 1), 1.0), -1.0))
       else
          th_bish(i) = acos(max(min(tvec(i, 1), 1.0), -1.0))
       end if
    end do
  end subroutine th_bishop

  !> FIXME : Add documentation
  pure subroutine B_pol(invR, rpgrad, Bpolmag, nth)
    implicit none
    integer, intent (in) :: nth
    real, dimension (-ntgrid:), intent (in) :: invR
    real, dimension (-ntgrid:, :), intent (in) :: rpgrad
    real, dimension (-ntgrid:)   , intent (out) :: Bpolmag
    integer :: i
    do i = -nth, nth
       Bpolmag(i) = sqrt(rpgrad(i, 1)**2 + rpgrad(i, 2)**2) * invR(i)
    end do
  end subroutine B_pol

  !> FIXME : Add documentation
  pure subroutine B_mod(invR, Bpolmag, Bmod, bi, nth)
    implicit none
    integer, intent (in) :: nth
    real, dimension (-ntgrid:), intent (in) :: invR, Bpolmag
    real, dimension(-ntgrid:), intent (out) :: Bmod
    real, intent(in) :: bi
    integer :: i
    do i = -nth, nth
       Bmod(i) = sqrt((bi * invR(i))**2 + Bpolmag(i)**2)
    end do
  end subroutine B_mod

  !> General calculation of the poloidal magnetic field radius of curvature
  !> (see section 2.3 of "GS2 analytic geometry specification") ! JRB
  !> @Note we actually calculate 1/Rpol as this is what we actually use.
  pure subroutine inv_Rpol(theta, th_bish, ltheta, invRpol, nth)
    use constants, only: twopi, pi
    implicit none
    integer, intent (in) :: nth
    real, dimension(-ntgrid:), intent (in) :: th_bish, theta, ltheta
    real, dimension(-ntgrid:), intent (out) :: invRpol
    real, dimension(-ntgrid:ntgrid) :: dthdl
    integer :: i

    ! Calculate du/dtheta (where u is defined according to Miller PoP 5 973 (1998))
    ! Note that we must properly account for discontinuities in u
    do i = -nth+1, nth-1
       if (th_bish(i+1)-th_bish(i-1) > pi) then
          invRpol(i) = (th_bish(i+1) - th_bish(i-1) - twopi) / (theta(i+1) - theta(i-1))
       else if(th_bish(i+1)-th_bish(i-1) < -pi) then
          invRpol(i) = (th_bish(i+1) - th_bish(i-1) + twopi) / (theta(i+1) - theta(i-1))
       else
          invRpol(i) = (th_bish(i+1) - th_bish(i-1)) / (theta(i+1) - theta(i-1))
       end if
    end do

    invRpol(nth) = (th_bish(-nth+1) - th_bish(nth-1)) / (theta(-nth+1) - theta(nth-1) + twopi)
    invRpol(-nth) = invRpol(nth)

    ! Calculate dtheta/dl (where l is the poloidal arc length)
    dthdl(-nth) = (theta(nth-1) - theta(-nth+1) - twopi) &
         /(ltheta(nth-1) - ltheta(-nth+1) - ltheta(nth) + ltheta(-nth))

    do i = -nth+1, nth-1
       dthdl(i) = (theta(i+1) - theta(i-1)) / (ltheta(i+1) - ltheta(i-1))
    end do

    dthdl(nth) = dthdl(-nth)

    ! Scale by dtdl to complete calculation of  1/Rpol
    ! (the inverse radius of curvature of the poloidal magnetic field)
    do i = -nth, nth
       invRpol(i) = (invRpol(i) * dthdl(i))
    end do
  end subroutine inv_Rpol

  !> FIXME : Add documentation
  subroutine finish_geometry
    implicit none
    if (allocated(geom)) then
       call geom%finalise
       deallocate(geom)
    end if
  end subroutine finish_geometry

  !> Initialse [[geometry:theta]] to a uniform grid with `2*((2*nperiod - 1)*(nt/2)) + 1`
  !> points. Also set [[geometry:nth]], [[geometry:ntheta]] and [[geometry::ntgrid]]
  !> consistently; importantly [[geometry:nth]] is forced to be even. If `nt` is odd, then
  !> `nth == nt - 1`.
  !>
  !> `nt` is set equal to [[geometry:ntheta]] on exit
  subroutine init_uniform_theta_grid(nt, theta)
    use constants, only: pi
    implicit none
    integer, intent(in out) :: nt
    real, dimension(:), allocatable, intent(in out) :: theta
    integer :: i

    nth = nt / 2
    ! Make sure ntheta is even
    ntheta = nth * 2
    nt = ntheta
    ntgrid = (2 * nperiod - 1) * nth

    if (verb>3) write(6,*) "init_theta: allocated(theta),ntgrid,nt,nperiod=",&
      allocated(theta),ntgrid,nt,nperiod

    if (allocated(theta)) deallocate (theta)
    allocate(theta(-ntgrid:ntgrid))

    theta = [ (i * pi / real(nth), i=-ntgrid, ntgrid) ]
  end subroutine init_uniform_theta_grid

  !> Copy the `[-nth,nth] = [-pi,pi]` section of the passed array into
  !> the `theta > +/- pi` parts, adding on `k * ext` (with k an integer
  !> labelling the sections/periods). Here we expect `ext == array(nth) - array(-nth)`.
  pure subroutine periodic_copy(a, ext, nperiod)
    implicit none
    real, dimension(-ntgrid:ntgrid), intent(in out) :: a
    real, intent(in) :: ext
    integer, intent(in) :: nperiod
    integer :: i, k, itot
    if (nperiod <= 1) return
    ! Note we end up enforcing a(nth)-a(-nth) = ext, but we should probably check this is
    ! also true on input to the routine
    do k = -nperiod + 1, -1
       do i = -nth, nth
          itot = i + k * ntheta
          a(itot) = a(i) + k * ext
       end do
    end do

    do k = 1, nperiod - 1
       do i = -nth, nth
          itot = i + k * ntheta
          a(itot) = a(i) + k * ext
       end do
    end do
  end subroutine periodic_copy

  !> FIXME : Add documentation
  pure subroutine bishop_gradB(Bmod, Bpolmag, invRpol, th_bish, ltheta, rmajor, dp_new, bi, gradB)
    implicit none
    real, dimension(-ntgrid:), intent (in) :: Bmod, Bpolmag, invRpol
    real, dimension(-ntgrid:), intent (in) :: th_bish, ltheta, rmajor
    real, intent(in) :: dp_new, bi
    real, dimension(-ntgrid:,:), intent(out) :: gradB
    real :: bp
    integer :: i
    do i = -nth, nth
       bp = -Bpolmag(i)
       gradB(i, 1) = (bp / Bmod(i)) * (bp * invRpol(i) + dp_new * rmajor(i) &
            - bi**2 * sin(th_bish(i)) / (rmajor(i)**3 * bp)) !dBdrho
    end do
    call gradl(ltheta, Bmod, gradB(:, 2), 0., nth) !dBdl
  end subroutine bishop_gradB

  !> Perform minor checks on consistency of flags
  subroutine check(geq, eeq, peq, leq, deq, ideq, ceq, gs2deq, teq, bishop)
    use mp, only: mp_abort
    implicit none
    logical, intent(in) :: geq, eeq, peq, leq, deq, ideq, ceq, gs2deq, teq
    integer, intent(in) :: bishop
    integer :: nflag_set
    if ((eeq .or. deq .or. gs2deq) .and. bishop == 0) &
         call mp_abort('gs2_eq/efit_eq/dfit_eq are not compatible with bishop = 0',.true.)
    nflag_set = count([geq, ceq, eeq, gs2deq, peq, teq, leq, deq, ideq])
    if (nflag_set == 0) then
       call mp_abort('Must set at least one equilibrium flag to true', .true.)
    else if (nflag_set > 1) then
       call mp_abort('More than one equilibrium flag to true, not allowed', .true.)
    end if
  end subroutine check

  !> Calculate rgrid for the case of eeq and deq: in this case, rgrid
  !> is the distance to the magnetic axis as a function of theta.
  subroutine rtg(rgrid, rp, theta)
    implicit none
    real, dimension(-ntgrid:), intent (out) :: rgrid
    real, intent (in) :: rp
    real, dimension(-ntgrid:), intent (in) :: theta
    integer :: i
    do i = -nth, nth
       rgrid(i) = geom%rfun(rp, theta(i))
    end do
    call periodic_copy(rgrid, 0.0, nperiod)
  end subroutine rtg

  !> Calculate the derivative of the flux label rho
  !> with respect to the internal flux label rp (the
  !> normalised poloidal flux)
  real function drho_drp(rp)
    implicit none
    real, intent(in) :: rp
    real :: rp1, rp2, rho1, rho2
    !     compute d rho / d rp
    rp1 = rp * (1 - delrho) ; rp2 = rp * (1 + delrho)
    if (irho == 1) then
       rho1 = sqrt(abs(phi(rp1) / phi(rpmax))) ; rho2 = sqrt(abs(phi(rp2) / phi(rpmax)))
    else if (irho == 2) then
       rho1 = 0.5 * diameter(rp1) ; rho2 = 0.5 * diameter(rp2)
    else if (irho == 3) then
       rho1 = psifun(rp1) ; rho2 = psifun(rp2)
    else if (irho == 4) then
       rho1 = rhofun(rp1) ; rho2 = rhofun(rp2)
    end if
    drho_drp = (rho2 - rho1) / (rp2 - rp1)
  end function drho_drp

  !> FIXME : Add documentation | Only used for writelots output and if irho=1
  real function drho_drhod(rp, drhodrp)
    implicit none
    real, intent (in) :: drhodrp, rp
    real :: rp1, rp2,  rho1, rho2
    !     compute d rho / d rb
    rp1 = rp * (1 - delrho) ; rp2 = rp * (1 + delrho)
    rho1 = 0.5 * diameter(rp1) ; rho2 = 0.5 * diameter(rp2)
    drho_drhod = drhodrp / ((rho2 - rho1) / (rp2 - rp1))
  end function drho_drhod

  !> FIXME : Add documentation
  pure real function fluxavg(f, theta, jacob)
    implicit none
    real, dimension(-nth:nth), intent(in) :: f, theta, jacob
    real, dimension(-nth:nth) :: delth
    delth(-nth+1:nth-1) = 0.5*(theta(-nth+2:nth) - theta(-nth:nth-2))
    delth(-nth) = 0.5*(theta(-nth+1) - theta(-nth))
    delth(nth) = 0.5*(theta(nth) - theta(nth-1))

    ! We might want to write this to use trapezoidal_integration instead
    ! of this sum form. Currently only used in [[f_trap]], which itself
    ! is only used here for one piece of information written to screen.
    ! As it doesn't get into the simulated system we probably don't mind
    ! too much beyond readability concerns.
    fluxavg = sum(f * jacob * delth) / sum(jacob * delth)
  end function fluxavg

  !> FIXME : Add documentation
  pure real function f_trap(b_mag, theta, jacob)
    implicit none
    real, dimension(-nth:nth), intent (in) :: b_mag, theta, jacob
    real :: ftu, ftl, havg, h2avg, h(-nth:nth)
    ! use expressions from Y. R. Lin-Liu and Bob Miller, coding from Q. P. Liu
    h = min(b_mag/maxval(b_mag),1.0)
    havg = fluxavg(h, theta, jacob)
    h2avg = fluxavg(h**2, theta, jacob)
    ftu = 1.0 - (h2avg / havg**2) * (1.0 - sqrt(1.0 - havg) * (1.0 + 0.5 * havg))
    ftl = 1.0 - h2avg * fluxavg((1.0 - sqrt(1.0 - h) * (1.0 + 0.5 * h)) / h**2, theta, jacob)
    f_trap = 0.75 * ftu + 0.25 * ftl
  end function f_trap

  !> FIXME : Add documentation
  subroutine plotdata (rgrid, seik, grads, dpsidrho, theta, &
       rplot, zplot, aplot, rprime, zprime, aprime)
    implicit none
    real, dimension (-ntgrid:), intent (in) :: rgrid, seik, theta
    real, dimension (-ntgrid:,:), intent (in) :: grads
    real, dimension (-ntgrid:), intent(out) :: rplot, zplot, aplot, rprime, zprime, aprime
    real, intent (in) :: dpsidrho
    real :: rplus, rminus, dr
    integer :: i
    ! At the present flux surface, calculate these three functions
    ! as a function of theta, together with their derivates with
    ! respect to rho:
    !
    ! R, Z, alpha-phi
    !
    ! Call these variables
    !
    ! Rplot, Zplot, aplot
    !
    ! and their derivatives:
    !
    ! Rprime, Zprime, aprime
    !
    ! BD: bug??
    ! Is this correct for EFIT data?  Need to check
    do i = -ntgrid, ntgrid
       Rplot(i) = Rpos(rgrid(i), theta(i))
       Zplot(i) = Zpos(rgrid(i), theta(i))
    end do
    aplot = seik

    do i = -ntgrid, ntgrid
       rplus  = rgrid(i) * (1 + delrho)
       rminus = rgrid(i) * (1 - delrho)
       dr = 2 * delrho * rgrid(i)
       Rprime(i) = (Rpos(rplus, theta(i)) - Rpos(rminus, theta(i))) / dr
       Zprime(i) = (Zpos(rplus, theta(i)) - Zpos(rminus, theta(i))) / dr
       aprime(i) = grads(i, 1) * dpsidrho
    end do
  end subroutine plotdata
end module geometry
