#include "unused_dummy.inc"

!> A module to provide utilities used in geometry calculations.
!> Mainly to provide an abstract type used as a base for specific
!> geometry implementations.
module geo_utils
  use splines, only: spline, periodic_spline
  implicit none

  public :: abstract_geo_type, abstract_eqfile_geo_type, abstract_eqfile_cart_geo_type
  public :: geo_input_type, flux_surface_type, mod2pi, sort, find_lower_bound
  public :: bilinear_interpolate, mxh_max_moments

  integer, parameter :: EQFILE_LENGTH = 800

  integer, parameter :: geo_type_miller = 0
  integer, parameter :: geo_type_global = 1
  integer, parameter :: geo_type_generalized_elongation = 2
  integer, parameter :: geo_type_fourier_series = 3
  integer, parameter :: geo_type_miller_extended_harmonic = 4

  !> Maximum number of moments for the Miller Extended Harmonic equilbrium
  integer, parameter :: mxh_max_moments = 5

  !> Helper type to wrap up some flux surface specific variables for
  !> analytic equilibrium type.
  type :: flux_surface_type
     !> Flux surface major radius, \(R_{0N}\), input variable: `Rmaj`
     real :: Rmaj
     !> Location of reference magnetic field, \(R_{geoN}\), input variable: `R_geo`
     real :: R_geo
     !> Minor radius, \(r_{\psi N}\), input variable: `rhoc`
     real :: r
     !> Step size for radial derivatives, \(\Delta r_{\psi N}\), input variable: `delrho`
     real :: dr
     !> Minor radius of shaping surface, \(a_{\psi N}\), input variable: `aSurf`
     real :: aSurf
     !> Horizontal Shafranov shift, \(dR_{0N}/dr_{\psi N}\), input variable: `shift`, `Raxis`
     real :: sHorz
     !> Vertical Shafranov shift, \(dZ_{0N}/dr_{\psi N}\), input variable: `shiftVert`, `Zaxis`
     real :: sVert
     !> First shaping mode strength, \(\Delta_m\), input variable: `akappa`, `delta2`, `deltam`
     real :: delm
     !> Second shaping mode strength, \(\Delta_n\), input variable: `tri`, `delta3`, `deltan`
     real :: deln
     !> First shaping mode derivative, \(d\Delta_m/dr_{\psi N}\), input variable: `akappri`, `deltampri`
     real :: delmp
     !> Second shaping mode derivative, \(d\Delta_n/dr_{\psi N}\), input variable: `tripri`, `deltanpri`
     real :: delnp
     !> First shaping mode tilt angle, \(\theta_m\), input variable: `thetak`, `theta2`, `thetam`
     real :: thm
     !> Second shaping mode tilt angle, \(\theta_n\), input variable: `thetad`, `theta2`, `thetan`
     real :: thn
     !> Safety factor, \(q\), input variable: `qinp`
     real :: q
     !> Magnetic shear, \(\hat{s}\), input variable: `s_hat_input`
     real :: shat
     !> Number of points in theta, input variable: `ntheta`
     integer :: nt
     !> Geometry specification type, input variable: `geoType`
     integer :: geoType
     !> The first poloidal mode number, input variable: `mMode`
     integer :: mMode
     !> The second poloidal mode number, input variable: `nMode`
     integer :: nMode
     !> Number of moments for MXH equilibrium, input variable: `n_mxh`
     integer :: n_mxh
     !> Zeroth cosine moment for MXH equilibrium, input variable: `c0_mxh`
     real :: c0_mxh
     !> Cosine moments for MXH equilibrium, input variable: `c_mxh`
     real, dimension(mxh_max_moments) :: c_mxh
     !> Sine moments for MXH equilibrium, input variable: `s_mxh`
     real, dimension(mxh_max_moments) :: s_mxh
     !> Radial derivatives of cosine moments for MXH equilibrium, input variable: `dc_mxh_dr`
     real, dimension(mxh_max_moments) :: dc_mxh_dr
     !> Radial derivatives of sine moments for MXH equilibrium, input variable: `ds_mxh_dr`
     real, dimension(mxh_max_moments) :: ds_mxh_dr
  end type flux_surface_type

  interface flux_surface_type
    module procedure new_flux_surface_type
  end interface flux_surface_type

  !> A derived type used to hold the inputs to the geo initialisation routine.
  !> By encapsulating all possible inputs in this type we can avoid having to
  !> have multiple initialisation interfaces which all types have to define.
  type :: geo_input_type
     !> Filename to read in for equilibria based on reading a file
     character(len = EQFILE_LENGTH) :: eqfile = 'EQFILE_NOT_SET'
     !> Upsampling factor used in efit, dfit and gs2d cases
     integer :: big = -1
     !> File containing normalisation factors for dfit
     character(len = EQFILE_LENGTH) :: eqnormfile = ' EQNORMFILE_NOT_SET'
     !> Target theta grid used for ideq
     real, dimension(:), allocatable :: theta
     !> Used for analytic local equiliria
     type(flux_surface_type) :: surf
     !> Theta shift used in analytic local equiliria
     real :: thShift
  end type geo_input_type

  type, abstract :: abstract_geo_type
     !> Name of the specific implementation
     character(len=:), allocatable :: type_name
     !> Has this instance been initialised
     logical :: initialised = .false.
     !> Does this implementation have the full or half theta grid? If only half
     !> then we assume up-down symmetry and use this to set gradients etc.
     logical :: has_full_theta_range = .false.
     !> Typically the number of radial and theta grid points.
     integer :: nr, nt
     real :: psi_0, psi_a, B_T, beta_0, R_mag, Z_mag, aminor
     !> 2D map of co-ordinates to R, Z or B
     real, allocatable, dimension(:,:) :: R_psi, Z_psi, B_psi
     !> 2D maps of theta and eqpsi
     real, allocatable, dimension(:,:) :: eqth, eqpsi_2d
     !> Minor radius gradient in cylindrical coordinates \((\hat{e}_R,
     !> \hat{e}_Z, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dpcart
     !> Theta gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dtcart
     !> B gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dbcart
     !> Minor radius gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dpbish
     !> Theta gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dtbish
     !> B gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dbbish
     !> Equilibrium psi (not normalised) and pressure grids
     real, allocatable, dimension(:) :: eqpsi, pressure
     !> Normalised psi, fp, q, beta, diameter and R_centre values on grid
     real, allocatable, dimension(:) :: psi_bar, fp, qsf, beta, diam, rc
   contains
     procedure :: initialise
     procedure(read_and_set_interface), deferred :: read_and_set
     procedure(noargs_interface), deferred :: finalise, finish_setup
     procedure(rtheta_interface), deferred :: rpos, zpos, psi, rfun
     procedure(rp_interface), deferred :: rcenter, diameter
     procedure(pbar_interface), deferred :: btori, dbtori, qfun, pfun, dpfun, betafun
     !> Calculate the derivative of rp w.r.t. R and Z and return
     !> the modulus sqrt(drp/dR ^ 2 + drp/dZ^2). I.e. return |grad rp|.
     !> Parameters are:
     !>   - rgrid: Submodule radial grid (r_circ for EFIT, rp otherwise)
     !>   - theta: grid of theta
     !>   - gradf(:,1): |grad psi|
     !>   - gradf(:,2): 0.0
     !>   - char: if char = 'R', return |grad pressure| instead
     !>   - char: if char = 'T', then return theta gradient in bishop form
     !>    - gradf(:,1): (dtheta/dZ d rp/dR - dtheta/dR drp/dZ)/ |grad rp|
     !>    - gradf(:,2): (dtheta/dR d rp/dR + dtheta/dZ drp/dZ)/ |grad rp|
     !>   - rp: value of rp on the flux surface where we want the grad
     !>   - nth: number of theta points
     !>   - ntgrid: lower index of theta array is -ntgrid
     procedure :: gradient
     procedure(eqitem_interface), deferred :: eqitem
     procedure :: calculate_gradients
     procedure :: eqdcart, eqdbish, invR, derm, rhofun => rhofun_null
     procedure :: alloc_common_arrays, dealloc_common_arrays
     procedure :: alloc_arrays, dealloc_arrays, hahm_burrell
     procedure :: alloc_special_arrays => alloc_special_arrays_null
     procedure :: dealloc_special_arrays => dealloc_special_arrays_null
  end type abstract_geo_type

  type, abstract, extends(abstract_geo_type) :: abstract_eqfile_geo_type
     character(len = EQFILE_LENGTH) :: filename
     type(spline) :: btori_spl, qsf_spl, pressure_spl, beta_spl, diam_spl, rcenter_spl
   contains
     procedure :: finalise => finalise_eqfile
     procedure :: finish_setup => finish_setup_eqfile
     procedure :: setup_special_splines => setup_special_splines_eqfile_null
     procedure :: delete_special_splines => delete_special_splines_eqfile_null
     procedure :: setup_splines => setup_splines_eqfile
     procedure :: delete_splines => delete_splines_eqfile
     procedure :: btori => btori_eqfile, dbtori => dbtori_eqfile
     procedure :: qfun => qfun_eqfile, betafun => betafun_eqfile
     procedure :: pfun => pfun_eqfile, dpfun => dpfun_eqfile
     procedure :: diameter => diameter_eqfile, rcenter => rcenter_eqfile
     procedure :: psi => psi_eqfile, rfun => psi_eqfile
     procedure :: rpos => rpos_eqfile, zpos => zpos_eqfile
     procedure :: eqitem => eqitem_eqfile
  end type abstract_eqfile_geo_type

  type, abstract, extends(abstract_eqfile_geo_type) :: abstract_eqfile_cart_geo_type
     integer :: nw, nh
     real, allocatable, dimension(:) :: thetab, r_bound, eq_R, eq_Z
     type(periodic_spline) :: rbound_spl
   contains
     procedure :: bound => bound_eqfile_cart, rfun, zbrent, shared_setup => shared_setup_cart
     procedure :: psi => psi_eqfile_cart, rpos => rpos_eqfile_cart, zpos => zpos_eqfile_cart
     procedure :: derm => derm_cart, eqitem => eqitem_cart, diameter => diameter_eqfile_cart
  end type abstract_eqfile_cart_geo_type

  ! Define the procedure interfaces.
  interface
     real function rtheta_interface(self, r, theta)
       import abstract_geo_type
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function rtheta_interface

     real function rp_interface(self, rp)
       import abstract_geo_type
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: rp
     end function rp_interface

     real function pbar_interface(self, pbar)
       import abstract_geo_type
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function pbar_interface

     real function eqitem_interface(self, r_in, theta_in, f, char) result(fstar)
       import abstract_geo_type
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r_in, theta_in
       real, dimension(:, :), intent(in) :: f
       character(len=1), intent(in) :: char
     end function eqitem_interface

     subroutine read_and_set_interface(self, inputs)
       import abstract_geo_type, geo_input_type
       class(abstract_geo_type), intent(in out) :: self
       type(geo_input_type), intent(in) :: inputs
     end subroutine read_and_set_interface

     subroutine noargs_interface(self)
       import abstract_geo_type
       class(abstract_geo_type), intent(in out) :: self
     end subroutine noargs_interface
  end interface

contains

  subroutine initialise(self, inputs, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid)
    class(abstract_geo_type), intent(in out) :: self
    type(geo_input_type), intent(in) :: inputs
    real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
    if (.not. self%initialised) then
       call self%read_and_set(inputs) ; call self%finish_setup
    end if
    psi_a_out = self%psi_a ; psi_0_out = self%psi_0
    rmaj = self%R_mag ; B_T0 = abs(self%B_T) ; avgrmid = self%aminor
    self%initialised = .true.
  end subroutine initialise

  pure subroutine alloc_arrays(self)
    class(abstract_geo_type), intent(in out) :: self
    call self%alloc_common_arrays ; call self%alloc_special_arrays
  end subroutine alloc_arrays

  pure subroutine dealloc_arrays(self)
    class(abstract_geo_type), intent(in out) :: self
    call self%dealloc_common_arrays ; call self%dealloc_special_arrays
  end subroutine dealloc_arrays

  pure subroutine alloc_special_arrays_null(self)
    class(abstract_geo_type), intent(in out) :: self
    UNUSED_DUMMY(self)
  end subroutine alloc_special_arrays_null

  pure subroutine dealloc_special_arrays_null(self)
    class(abstract_geo_type), intent(in out) :: self
    UNUSED_DUMMY(self)
  end subroutine dealloc_special_arrays_null

  pure subroutine alloc_common_arrays(self)
    class(abstract_geo_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt), self%B_psi(nr, nt))
    allocate(self%dpcart(nr, nt, 2), self%dbcart(nr, nt, 2), self%dbbish(nr, nt, 2))
    allocate(self%dtcart(nr, nt, 2), self%dpbish(nr, nt, 2), self%dtbish(nr, nt, 2))
    allocate(self%eqpsi(nr), self%pressure(nr), self%eqth(nr, nt), self%eqpsi_2d(nr, nt))
    allocate(self%psi_bar(nr), self%fp(nr), self%qsf(nr), self%beta(nr), self%rc(nr), self%diam(nr))
  end subroutine alloc_common_arrays

  pure subroutine dealloc_common_arrays(self)
    class(abstract_geo_type), intent(in out) :: self
    deallocate(self%R_psi, self%Z_psi, self%B_psi, self%eqth, self%eqpsi, self%pressure)
    deallocate(self%dpcart, self%dpbish, self%dtcart, self%dtbish, self%dbcart, self%dbbish)
    deallocate(self%eqpsi_2d, self%psi_bar, self%fp, self%qsf, self%beta, self%rc, self%diam)
  end subroutine dealloc_common_arrays

  !> Given theta, R, Z, B and psi on 2D grids calculate the index space
  !> derivatives in the two grid dimensions and use these to find gradients
  !> in cartesian and Bishop space.
  pure subroutine calculate_gradients(self)
    implicit none
    class(abstract_geo_type), intent(in out) :: self
    ! The following quantities represent changes between adjacent grid
    ! points. The first two dimensions are usually minor radius and theta,
    ! respectively. The last dimension is the difference along either
    ! the (minor) radial or theta directions, respectively
    !> Major radius differences on theta, minor radius grids
    real, allocatable, dimension (:,:,:) :: drm
    !> Vertical differences on theta, minor radius grids
    real, allocatable, dimension (:,:,:) :: dzm
    !> Temporary version used to hold minor radius(psi), theta or B index space derivatives
    real, allocatable, dimension (:,:,:) :: dtmp_m

    ! Need drm and dzm for all eqdcart calls
    allocate(drm(self%nr, self%nt, 2), dzm(self%nr, self%nt, 2))
    call self%derm(self%R_psi, drm, 'E')
    call self%derm(self%Z_psi, dzm, 'O')

    ! grad(psi) in cartesian and bishop form
    allocate(dtmp_m(self%nr, self%nt, 2))
    call self%derm(self%eqpsi_2d, dtmp_m, 'E')
    call self%eqdcart(dtmp_m, drm, dzm, self%dpcart)
    call self%eqdbish(self%dpcart, self%dpcart, self%dpbish)

    ! grad(B) in cartesian and bishop form
    call self%derm(self%B_psi, dtmp_m, 'E')
    call self%eqdcart(dtmp_m, drm, dzm, self%dbcart)
    call self%eqdbish(self%dbcart, self%dpcart, self%dbbish)

    ! grad(theta) in cartesian and bishop form
    call self%derm(self%eqth, dtmp_m, 'T')
    call self%eqdcart(dtmp_m, drm, dzm, self%dtcart)
    call self%eqdbish(self%dtcart, self%dpcart, self%dtbish)
    deallocate(dtmp_m, drm, dzm)
  end subroutine calculate_gradients

  !> Converts derivatives w.r.t. (psi_index,theta_index) to derivatives
  !> w.r.t. (R,Z). Note that `dfcart(1, :, :)` is not valid
  pure subroutine eqdcart(self, dfm, drm, dzm, dfcart)
    implicit none
    class(abstract_geo_type), intent(in) :: self
    !> dfm(:,:,1) is deriv w.r.t. psi_index (i.e. (df/dpsi)_theta * delta psi);
    !> dfm(:,:,2) is deriv w.r.t. theta_index
    real, dimension (:,:,:), intent(in)  :: dfm, drm, dzm
    !> dfcart(:,:,1) is deriv w.r.t. R;
    !> dfcart(:,:,2) is deriv w.r.t. Z
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm, 1), size(dfm, 2)) :: denom
    denom = drm(:,:,1) * dzm(:,:,2) - drm(:,:,2) * dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * dzm(:,:,2) - dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * drm(:,:,2) + drm(:,:,1) * dfm(:,:,2)
    ! Skips first point in first dimension because it's likely exactly on the
    ! magnetic axis, which means the Jacobian is zero.
    dfcart(2:, :, 1) = dfcart(2:, :, 1) / denom(2:, :)
    dfcart(2:, :, 2) = dfcart(2:, :, 2) / denom(2:, :)
    UNUSED_DUMMY(self)
  end subroutine eqdcart

  !> Convert gradients of a function f w.r.t. R,Z into bishop form.
  !> Note that `dbish(1, :, :)` is not valid
  pure subroutine eqdbish(self, dcart, dpcart, dbish)
    implicit none
    class(abstract_geo_type), intent(in) :: self
    !> dcart(:,:,1) is gradient of f w.r.t. R;
    !> dcart(:,:,2) is gradient of f w.r.t. Z
    real, dimension(:, :, :), intent (in) :: dcart, dpcart
    !> dbish(:,:,1) is set to (df/dR dpsi/dR + df/dZ dpsi/dZ)/|grad psi|;
    !> dbish(:,:,2) is set to (-df/dR dpsi/dZ + df/dZ dpsi/dR)/|grad psi|.
    !> Note that in the special case where f=psi:
    !> dbish(:,:,1) is |grad psi|; dbish(:,:,2) is 0
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart, 1) - 1, size(dcart, 2)) :: denom
    ! denom is |grad psi|
    denom = hypot(dpcart(2:,:,1), dpcart(2:,:,2))
    dbish(:,:,1) = dcart(:,:,1) * dpcart(:,:,1) + dcart(:,:,2) * dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * dpcart(:,:,2) + dcart(:,:,2) * dpcart(:,:,1)
    dbish(2:, :, 1) = dbish(2:, :, 1) / denom
    dbish(2:, :, 2) = dbish(2:, :, 2) / denom
    UNUSED_DUMMY(self)
  end subroutine eqdbish

  subroutine gradient(self, rgrid, theta, grad, char, rp, nth_used, ntm, use_bishop)
    use splines, only: inter_d_cspl
    implicit none
    class(abstract_geo_type), intent(in) :: self
    integer, intent (in) :: nth_used, ntm
    character(1), intent (in) :: char
    real, dimension (-ntm:), intent (in)  :: rgrid, theta
    real, dimension (-ntm:,:), intent (out) :: grad
    real, intent(in) :: rp
    logical, intent(in) :: use_bishop
    real, dimension(:, :, :), allocatable :: darr
    real, dimension (1) :: p_eqpsi, dp_deqpsi
    integer :: i
    if (use_bishop) then
       select case(char)
       case('B')
          allocate(darr, source = self%dbbish)
       case('P', 'R')
          allocate(darr, source = self%dpbish)
       case('T')
          allocate(darr, source = self%dtbish)
       end select
    else
       select case(char)
       case('B')
          allocate(darr, source = self%dbcart)
       case('P', 'R')
          allocate(darr, source = self%dpcart)
       case('T')
          allocate(darr, source = self%dtcart)
       end select
    end if

    do i = -nth_used, nth_used
       grad(i, 1) = self%eqitem(rgrid(i), theta(i), darr(:,:,1), 'R')
       grad(i, 2) = self%eqitem(rgrid(i), theta(i), darr(:,:,2), 'Z')
    end do

    if(char == 'T' .and. .not. self%has_full_theta_range) then
       where (theta(-nth_used:nth_used) < 0.0)
          grad(-nth_used:nth_used, 1) = -grad(-nth_used:nth_used, 1)
          grad(-nth_used:nth_used, 2) = -grad(-nth_used:nth_used, 2)
       end where
    end if

    !     to get grad(pressure), multiply grad(psi) by dpressure/dpsi
    ! Could use self%dpfun if took pbar rather than rp?
    if (char == 'R') then !Only actually used if bishop == 0
       call inter_d_cspl(self%eqpsi, self%pressure, [rp], p_eqpsi, dp_deqpsi)
       do i = -nth_used, nth_used
          grad(i, 1:2) = grad(i, 1:2) * dp_deqpsi(1) * 0.5 * self%beta_0
       end do
    end if
  end subroutine gradient

  !> Calculate the derivative of f w.r.t to the radial
  !> and poloidal indexes (i.e. calculate the finite
  !> differences without dividing by
  !> delta psi and delta theta).
  !> - dfm(:,:,1) is the psi derivative
  !> - dfm(:,:,2) is the theta derivative
  !> - char specifies the periodicity in theta
  !>   - 'E', 'O' mean continuous at theta = +/- pi
  !>   - 'T' means a jump of 2 pi at theta = +/- pi
  pure subroutine derm(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(abstract_geo_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nr, nt
    nr = self%nr ; nt = self%nt

    i = 1
    dfm(i,:,1) = -0.5 * (3 * f(i,:) - 4 * f(i+1,:) + f(i+2,:))

    i = nr
    dfm(i,:,1) = 0.5 * (3 * f(i,:) - 4 * f(i-1,:) + f(i-2,:))

    do i = 2, nr-1
       dfm(i,:,1) = 0.5 * (f(i+1,:) - f(i-1,:))
    end do

    do j = 2, nt-1
       dfm(:,j,2) = 0.5 * (f(:,j+1) - f(:,j-1))
    end do

    ! Handle periodic boundary
    if (self%has_full_theta_range) then
       ! Compute df/dtheta at -pi,+pi assuming
       !   (i) periodicity in theta,
       !  (ii) endpoints at -pi,+pi
       ! (iii) equispaced grid in theta, and
       !  (iv) 1 and nt correspond to -pi, +pi and are at the same place
       dfm(:,1,2) = 0.5 * (f(:,2) - f(:,nt-1))
       dfm(:,nt,2) = dfm(:, 1, 2)
       if (char == 'T') then
          dfm(:,1,2) = dfm(:,1,2) + pi
          dfm(:,nt,2) = dfm(:,nt,2) + pi
       end if
    else
       ! assume up-down symmetry for now:
       select case (char)
       case ('E')
          dfm(:,1,2) = 0.0
          dfm(:,nt,2) = 0.0
       case ('O')
          dfm(:,1,2) = f(:,2)
          dfm(:,nt,2) = -f(:,nt-1)
       case ('T')
          dfm(:,1,2) = f(:,2) !f(:,j+1) - f(:,j) -- assumes f(:,1) = 0
          dfm(:,nt,2) = pi - f(:,nt-1) !f(:, j) - f(:, j-1)  -- assumes f(:,nt) = pi
       end select
    end if
  end subroutine derm

  real function invR (self, r, theta)
    implicit none
    class(abstract_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%rpos(r, theta)
  end function invR

  real function rhofun_null (self, pbar) result(rhofun) !Can't be pure as deq variant isnt
    implicit none
    class(abstract_geo_type), intent(in) :: self
    real, intent (in) :: pbar
    rhofun = 0.
    UNUSED_DUMMY(self); UNUSED_DUMMY(pbar)
  end function rhofun_null

  !> FIXME : Add documentation
  subroutine hahm_burrell(self, a)
    use file_utils, only: open_output_file, close_output_file
    use mp, only: proc0
    implicit none
    class(abstract_geo_type), intent(in) :: self
    real, intent(in) :: a
    integer :: i, fort24_unit, nr
    real :: gradpsi, mag_B, rho_eq, rp1, rp2, rho1, rho2, drhodpsiq
    real, dimension(self%nr) :: gamma_shear, pbar, dp, d2p, pres, bs_coll, s__hat, bs_sh
    if (.not. proc0) return
    nr = self%nr
    gamma_shear = 0.
    pbar = (self%eqpsi - self%eqpsi(1)) / (self%eqpsi(nr) - self%eqpsi(1))

    do i = 1, nr
       dp(i)  = self%dpfun(pbar(i))
       pres(i) = self%pfun(pbar(i))
    end do

    d2p(1) = 0.0 ; d2p(nr) = 0.0
    do i = 2, nr - 1
       d2p(i) = (dp(i+1) - dp(i-1)) / (self%eqpsi(i+1) - self%eqpsi(i-1))
    end do

    do i = 2, nr- 1
       rp1 = self%eqpsi(i+1) ; rp2 = self%eqpsi(i-1)
       rho1 = 0.5 * self%diameter(rp1) ; rho2 = 0.5 * self%diameter(rp2)
       drhodpsiq = (rho1 - rho2) / (rp1 - rp2)
       gradpsi = self%eqitem(self%eqpsi(i), 0., self%dpbish(:,:,1), 'R')

       mag_B = hypot(self%btori(pbar(i)), gradpsi) * self%invR(self%eqpsi(i), 0.)
       !For geq and ideq could be mag_B = geom%eqitem(geom%eqpsi(i), 0., geom%B_psi, 'R')

       gamma_shear(i) = 0.01 * gradpsi**2 * (d2p(i) / pres(i) - a * (dp(i) / pres(i))**2) &
            * (2 * pres(i) / self%beta_0)**((1 - a) / 2) * (-pres(i) / (dp(i) / drhodpsiq)) &
            / mag_B

       s__hat(i) = (self%qfun(pbar(i+1)) - self%qfun(pbar(i-1))) / (rho2 - rho1)
       s__hat(i) = s__hat(i) * 0.5 * (rho1 + rho2) / self%qfun(pbar(i))
    end do

    ! Assume rho_*0 = 0.01, nu_*0 = 1.0
    do i = 2, nr - 1
       bs_coll(i) = 0.01 * ((2 * pres(i) / self%beta_0)**(1 - a))**2.5 &
             * self%qfun(0.) * 0.5 * self%diameter(self%eqpsi(i)) &
            * (-dp(i) / pres(i)) * self%R_psi(1,1)**1.5 / (2 * pres(i) / self%beta_0)**a
       bs_sh(i) = (-pres(i) / (dp(i) / drhodpsiq))*s__hat(i) / (self%R_psi(1,1) &
            * self%qfun(pbar(i)) * sqrt(pres(i)))
    end do

    call open_output_file(fort24_unit,".fort.24")
    do i = 2, nr - 1
       rho_eq = 0.5 * self%diameter(self%eqpsi(i))
       write(fort24_unit,'(i5,11(1x,e16.9))') i, pbar(i), rho_eq, pres(i), gamma_shear(i), bs_coll(i), bs_sh(i)
    end do
    call close_output_file(fort24_unit)
  end subroutine hahm_burrell

  !> Calculates fstar which is f interpolated at the location (r,theta).
  !> Here r is the normalised poloidal flux coordinate rp (= psi_pN + psi_0N)
  !> and theta_in is theta. f is a grid of values of f as a function of psi_p,theta
  !> First we find the radial and theta indices bounding our requested point and we
  !> then interpolate within the bound rectangle.
  real function eqitem_eqfile(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: r_in, theta_in
    real, dimension(:,:), intent(in) :: f
    character(1), intent(in) :: char
    integer :: istar, jstar
    real :: r, thet, eqpsi_sign, tps
    ! allow psi(r) to be a decreasing function
    eqpsi_sign = 1.
    if (self%eqpsi(2) < self%eqpsi(1)) eqpsi_sign = -1.
    associate (dim_1 => eqpsi_sign * self%eqpsi , mtheta => self%eqth(1, :)) !Non-contiguous copy - slow?

      !Copy input r to avoid any modification
      r = r_in ; thet = mod2pi(theta_in) ; tps = 1.

      if (.not. self%has_full_theta_range) then ! assume up-down symmetry
         if (char == 'Z' .and. abs(thet) >= 1.e-10) tps = sign(1.0, thet)
         thet = abs(thet)
      end if

      ! check for axis evaluation, no gradients available at axis, so do not get too close
      if (r < dim_1(2)) then
         write(*,*) r, theta_in, eqpsi_sign, dim_1(1), dim_1(2)
         call mp_abort('no evaluation too close to or inside axis allowed in eqitem', .true.)
      end if

      ! disallow evaluations on or outside the surface for now
      if (r >= dim_1(self%nr)) then
         write(*,*) r, theta_in, dim_1(self%nr), eqpsi_sign
         call mp_abort('No evaluation of eqitem allowed on or outside surface', .true.)
      end if

      istar = find_lower_bound(r, dim_1, 0)
      jstar = find_lower_bound(thet, mtheta, self%nt - 1)
      fstar = tps * bilinear_interpolate(f, dim_1, mtheta, istar, jstar, r, thet)
    end associate
  end function eqitem_eqfile

  subroutine setup_splines_eqfile(self)
    use splines, only: new_spline
    class(abstract_eqfile_geo_type), intent(in out) :: self
    self%diam_spl = new_spline(self%eqpsi, self%diam)
    self%rcenter_spl = new_spline(self%eqpsi, self%rc)
    self%btori_spl = new_spline(self%psi_bar, self%fp)
    self%qsf_spl = new_spline(self%psi_bar, self%qsf)
    self%pressure_spl = new_spline(self%psi_bar, self%pressure)
    self%beta_spl = new_spline(self%psi_bar, self%beta)
    call self%setup_special_splines
  end subroutine setup_splines_eqfile

  subroutine delete_splines_eqfile(self)
    use splines, only: delete_spline
    class(abstract_eqfile_geo_type), intent(in out) :: self
    call delete_spline(self%diam_spl)
    call delete_spline(self%rcenter_spl)
    call delete_spline(self%btori_spl)
    call delete_spline(self%qsf_spl)
    call delete_spline(self%pressure_spl)
    call delete_spline(self%beta_spl)
    call self%delete_special_splines
  end subroutine delete_splines_eqfile

  subroutine setup_special_splines_eqfile_null(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    UNUSED_DUMMY(self)
  end subroutine setup_special_splines_eqfile_null

  subroutine delete_special_splines_eqfile_null(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    UNUSED_DUMMY(self)
  end subroutine delete_special_splines_eqfile_null

  subroutine finalise_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    self%initialised = .false.
    call self%dealloc_arrays ; call self%delete_splines
  end subroutine finalise_eqfile

  subroutine finish_setup_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    call self%calculate_gradients ; call self%setup_splines
  end subroutine finish_setup_eqfile

  real function btori_eqfile (self, pbar) result(btori)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    btori = splint(pbar, self%btori_spl)
  end function btori_eqfile

  !> \(\frac{dI}{d\psi}\) at the given normalised \(\psi\)
  real function dbtori_eqfile (self, pbar) result(dbtori)
    use splines, only: dsplint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    dbtori = dsplint(pbar, self%btori_spl) / (self%psi_a - self%psi_0)
  end function dbtori_eqfile

  real function qfun_eqfile (self, pbar) result(qfun)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    qfun = splint(pbar, self%qsf_spl)
  end function qfun_eqfile

  real function pfun_eqfile (self, pbar) result(pfun)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    pfun = 0.5 * self%beta_0 * splint(pbar, self%pressure_spl)
  end function pfun_eqfile

  real function dpfun_eqfile (self, pbar) result(dpfun)
    use splines, only: dsplint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    ! p_N would be B**2/mu_0 => p = beta/2 in our units
    dpfun = 0.5 * self%beta_0 * dsplint(pbar, self%pressure_spl) / (self%psi_a - self%psi_0)
  end function dpfun_eqfile

  !> Return the diameter of the flux surface at a given major radius
  real function diameter_eqfile (self, rp) result(diameter)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: rp
    diameter = splint(rp, self%diam_spl)
  end function diameter_eqfile

  !> Return the major radius of the centre of a given flux surface
  real function rcenter_eqfile (self, rp) result(rcenter)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: rp
    rcenter = splint(rp, self%rcenter_spl)
  end function rcenter_eqfile

  real function psi_eqfile (self, r, theta) result(psi) !Can't mark pure as eqfile_cart variant is not pure
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    if (r==0.) then
       psi = self%psi_0
    else if (r==1. .and. theta == 0.0) then
       psi = self%psi_a
    else
       psi = r
    end if
  end function psi_eqfile

  real function betafun_eqfile (self, pbar) result(betafun)
    use splines, only: splint
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent(in) :: pbar
    if (pbar == 0.) then
       betafun = self%beta_0
    else
       betafun = splint(pbar, self%beta_spl)
    end if
  end function betafun_eqfile

  real function rpos_eqfile (self, r, theta) result(rpos)
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    rpos = self%eqitem(r, theta, self%R_psi, 'R')
  end function rpos_eqfile

  real function zpos_eqfile (self, r, theta) result(zpos)
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    zpos = self%eqitem(r, theta, self%Z_psi, 'Z')
  end function zpos_eqfile

  !> Calculate the derivative of f w.r.t. R, Z
  !> - dfm(:,:,1) is deriv w.r.t. R
  !> - dfm(:,:,2) is deriv w.r.t. Z
  pure subroutine derm_cart(self, f, dfm, char)
    use constants, only: pi
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, dimension(:,:), intent(in) :: f
    real, dimension(:,:,:), intent(out) :: dfm
    character(1), intent(in) :: char
    integer :: i, j, nw, nh
    nw = self%nw ; nh = self%nh
    i = 1
    dfm(i,:,1) = -0.5 * (3 * f(i,:) - 4 * f(i+1,:) + f(i+2,:))

    i = nw
    dfm(i,:,1) = 0.5 * (3 * f(i,:) - 4 * f(i-1,:) + f(i-2,:))

    j = 1
    dfm(:,j,2) = -0.5 * (3 * f(:,j) - 4 * f(:,j+1) + f(:,j+2))

    j = nh
    dfm(:,j,2) = 0.5 * (3 * f(:,j)-4 * f(:,j-1) + f(:,j-2))

    do i = 2, nw-1
       dfm(i,:,1) = 0.5 * (f(i+1,:) - f(i-1,:))
    end do

    do j = 2, nh-1
       dfm(:,j,2) = 0.5 * (f(:, j+1) - f(:, j-1))
    end do

    if (char == 'T') then !Could a similar approach be ported to derm_eqfile?
       where(dfm(:, :, 2) > pi * 0.5)
          dfm(:, :, 2) = dfm(:, :, 2) - pi
       end where
    end if
  end subroutine derm_cart

  !> fstar is f(R,Z) interpolated at the values (r,thetain). The parameter
  !> r is the distance to the magnetic axis.
  real function eqitem_cart(self, r_in, theta_in, f, char) result(fstar)
    use mp, only: mp_abort
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: r_in, theta_in
    real, dimension(:, :), intent(in) :: f
    character(1), intent(in) :: char
    integer :: istar, jstar
    real :: r, thet, r_pos, z_pos
    associate(dim_1 => self%eq_R, dim_2 => self%eq_Z)
      !Copy input r to avoid any modification
      r = r_in ; thet = mod2pi(theta_in)
      r_pos = self%Rpos(r, thet) ; z_pos = self%Zpos(r, thet)

      ! ensure point is on R mesh
      if (r_pos >= dim_1(self%nw) .or. r_pos <= dim_1(1)) then
         write(*,*) r, theta_in, dim_1(self%nw), r_pos
         call mp_abort('No evaluation of eqitem allowed outside or on edge of R domain', .true.)
      end if

      ! ensure point is on Z mesh
      if (z_pos >= dim_2(self%nh) .or. z_pos <= dim_2(1)) then
         write(*,*) r, theta_in, dim_2(1), dim_2(self%nh), z_pos
         call mp_abort('No evaluation of eqitem allowed outside or on edge of Z domain', .true.)
      end if

      istar = find_lower_bound(r_pos, dim_1, 0)
      jstar = find_lower_bound(z_pos, dim_2, 0)
      fstar = bilinear_interpolate(f, dim_1, dim_2, istar, jstar, r_pos, z_pos)
    end associate
  end function eqitem_cart

  pure integer function find_lower_bound(val, array, default) result(bound)
    real, intent(in) :: val
    real, dimension(:), intent(in) :: array
    integer, intent(in) :: default
    integer :: limit, ind
    limit = size(array) ; bound = default
    do ind = 2, limit
       if (val < array(ind)) then
          bound = ind - 1 ; exit
       end if
    end do
  end function find_lower_bound

  !> use opposite area stencil to interpolate to point bound by [istar,istar+1]
  !> and [jstar, jstar + 1].
  pure real function bilinear_interpolate(f, dim_1, dim_2, istar, jstar, &
       d1, d2) result(fstar)
    real, dimension(:, :), intent(in) :: f
    real, dimension(:), intent(in) :: dim_1, dim_2
    integer, intent(in) :: istar, jstar
    real, intent(in) :: d1, d2
    real :: sr, st, dr, dt
    dr = d1 - dim_1(istar) ; sr = dim_1(istar + 1) - d1
    dt = d2 - dim_2(jstar) ; st = dim_2(jstar + 1) - d2
    fstar = f(istar    , jstar    ) * sr * st &
         +f(istar + 1, jstar    ) * dr * st &
         +f(istar    , jstar + 1) * sr * dt &
         +f(istar + 1, jstar + 1) * dr * dt
    fstar = fstar / (abs(dim_1(istar + 1) - dim_1(istar)) * &
         (dim_2(jstar + 1) - dim_2(jstar)))
  end function bilinear_interpolate

  !> Find the value of the major radius on the plasma boundary
  !> at given geometric theta.
  real function bound_eqfile_cart(self, theta) result(bound)
    use splines, only: periodic_splint
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent(in) :: theta
    bound = periodic_splint(theta, self%rbound_spl)
  end function bound_eqfile_cart

  pure real function rpos_eqfile_cart (self, r, theta) result(rpos)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    rpos = self%R_mag + r * cos(theta)
  end function rpos_eqfile_cart

  pure real function zpos_eqfile_cart (self, r, theta) result(zpos)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    zpos = self%Z_mag + r * sin(theta)
  end function zpos_eqfile_cart

  real function psi_eqfile_cart (self, r, theta) result(psi)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    if (r==0.) then
       psi = self%psi_0
    else if (r==1. .and. theta == 0.0) then
       psi = self%psi_a
    else
       psi = self%eqitem(r, theta, self%eqpsi_2d, '?')
    end if
  end function psi_eqfile_cart

  real function diameter_eqfile_cart (self, rp) result(diameter)
    use constants, only: pi
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: rp
    diameter = self%rfun(rp, 0.) + self%rfun(rp, pi)
  end function diameter_eqfile_cart

  subroutine shared_setup_cart(self, nw_in, nh_in, spsi_bar, sefit_R, sefit_Z, f, p, q, &
       sefit_psi, nbbbs, rbbbs, zbbbs, bcentr, rleft, rwid, zoffset, zhei, big, calc_aminor)
    use splines, only: inter_cspl, fitp_surf2, fitp_surf1
    use constants, only: twopi
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in out) :: self
    integer, intent(in) :: nw_in, nh_in, nbbbs, big
    real, dimension(:), intent(in) :: spsi_bar, sefit_R, sefit_Z, f, p, q
    real, dimension(:), intent(in out) :: rbbbs, zbbbs
    real, intent(in) :: bcentr, rleft, rwid, zoffset, zhei
    real, dimension(:, :), intent(in) :: sefit_psi
    logical, intent(in) :: calc_aminor
    real, dimension(:), allocatable :: zp, temp, zx1, zxm, zy1, zyn, eq_R, eq_Z
    real :: zxy11, zxym1, zxy1n, zxymn, psi_N, f_N
    integer :: i, j, ierr
    self%nw = nw_in * big ; self%nr = self%nw
    self%nh = nh_in * big ; self%nt = self%nh
    call self%alloc_arrays()
    do i = 1, self%nw
       self%psi_bar(i) = float(i - 1) / float(self%nw - 1)
    end do
    call inter_cspl(spsi_bar, f, self%psi_bar, self%fp)
    call inter_cspl(spsi_bar, p, self%psi_bar, self%pressure)
    call inter_cspl(spsi_bar, q, self%psi_bar, self%qsf)

    allocate(eq_R(self%nw), eq_Z(self%nh))
    do i = 1, self%nw
       eq_R(i) = rleft + rwid * float(i - 1) / float(self%nw - 1)
    end do

    do j = 1, self%nh
       eq_Z(j) = zoffset + ((float(j - 1) / float(self%nh - 1) )) * zhei
    end do

    ! Setup 2D interpolation of psi(R, Z)
    allocate(zp(3 * nw_in * nh_in), temp(nw_in + 2 * nh_in))
    allocate(zx1(nh_in), zxm(nh_in), zy1(nw_in), zyn(nw_in))
    call fitp_surf1(nw_in, nh_in, sefit_R, sefit_Z, sefit_psi, &
         nw_in, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, &
         255, zp, temp, 1., ierr)

    do j = 1, self%nh
       do i = 1, self%nw
          self%eqpsi_2d(i, j) = fitp_surf2(eq_R(i), eq_Z(j), nw_in, nh_in, &
               sefit_R, sefit_Z, sefit_psi, nw_in, zp, 1.)
       end do
    end do
    deallocate(zp, temp)

    allocate(self%thetab(nbbbs), self%r_bound(nbbbs))
    self%thetab = atan2((zbbbs - self%Z_mag), (rbbbs - self%R_mag))
    self%r_bound = hypot((rbbbs - self%R_mag), (zbbbs - self%Z_mag))
    call sort(self%thetab, self%r_bound, zbbbs, rbbbs)

    ! Allow for duplicated points near +- pi:
    if (self%thetab(1) == self%thetab(2)) then
       self%thetab(1) = self%thetab(1) + twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    if (self%thetab(nbbbs-1) == self%thetab(nbbbs)) then
       self%thetab(nbbbs) = self%thetab(nbbbs) - twopi
       call sort(self%thetab, self%r_bound, zbbbs, rbbbs)
    end if

    ! It isn't likely that a duplicate point would exist near theta = 0,
    ! so I am not allowing this possibility for now.
    do i = 1, nbbbs-1
       if (self%thetab(i+1) == self%thetab(i)) then
          ! put in kluge for duplicate points, which happens near theta=0:
          self%thetab(i+1) = self%thetab(i+1) +  1.e-8
       end if
    end do

    if (calc_aminor) call a_minor(rbbbs, zbbbs, self%Z_mag, self%aminor)

    self%B_psi = 0.0 ; self%diam = 0.0 ; self%rc = 0.0
    self%r_bound = self%r_bound / self%aminor
    self%R_mag = self%R_mag / self%aminor
    self%Z_mag = self%Z_mag / self%aminor
    eq_R = eq_R / self%aminor ; eq_Z = eq_Z / self%aminor
    self%eq_R = eq_R ; self%eq_Z = eq_Z
    do j = 1, self%nh
       self%R_psi(:, j) = eq_R ; self%Z_psi(:, j) = eq_Z(j)
    end do

    self%B_T = abs(bcentr)
    psi_N = self%B_T * self%aminor**2
    self%psi_a = self%psi_a / psi_N
    self%psi_0 = self%psi_0 / psi_N
    self%eqpsi_2d = self%eqpsi_2d / psi_N
    f_N = self%B_T * self%aminor
    self%fp = self%fp / f_N
    ! MKS: beta = 2 mu_0 p / B**2
    self%beta = 4. * twopi * self%pressure * 1.e-7 / self%B_T**2
    self%beta_0 = self%beta(1)
    self%pressure = self%pressure / self%pressure(1)

    do j = 1, self%nh
       do i = 1, self%nw
          if (eq_Z(j) == self%Z_mag .and. eq_R(i) == self%R_mag) then
             self%eqth(i, j) = 0.  ! value should not matter
          else
             self%eqth(i, j) = atan2( (eq_Z(j) - self%Z_mag), (eq_R(i) - self%R_mag))
          end if
       end do
    end do
    self%eqpsi = self%psi_bar * (self%psi_a - self%psi_0) + self%psi_0
    self%has_full_theta_range = .true.
    if (.true.) then
       write (*,*) "Finished shared_setup... imported ",self%type_name," equilibrium"
       write (*,*) 'Some important quantities:'
       write (*,*) "aminor", self%aminor
       write (*,*) 'R_mag', self%R_mag
       write (*,*) 'B_T0', self%B_T
       write (*,*) 'beta', self%beta_0
    end if
  end subroutine shared_setup_cart

  !> FIXME : Add documentation | Designed for eeq/gs2deq only
  subroutine a_minor(r, z, Z_mag, a)
    use mp, only: mp_abort
    use splines, only: new_spline, splint, delete_spline, spline
    implicit none
    real, dimension(:), intent (in) :: r, z
    real, intent(in) :: Z_mag
    real, intent(out) :: a
    real :: r1, r2
    integer, parameter :: nz = 5, half_nz = int(nz/2.0), invalidPoint = -nz
    real, dimension(nz) :: rtmp, ztmp
    integer :: i, j, i1, n, k, index
    type (spline) :: spl
    !CMR, 28/10/08: add code to avoid duplicate points in 5 point spline
    !               to determine r2 on inboard mid-plane
    logical, parameter :: debug=.false.
    k = 0 ; n = size(r)
    if (debug) write(6,*) "aminor:"
    if (debug) write(6,fmt='("r=",10f8.4)') r
    if (debug) write(6,fmt='("z=",10f8.4)') z

    if (n < nz) call mp_abort('nbbbs < nz -- very strange.  Stopping. Look in geo_utils.fpp', .true.)

    j = 0
    do i = half_nz + 1, 1, -1
       j = j + 1
       ztmp(j) = z(i) ; rtmp(j) = r(i)
    end do

    if (r(1) == r(n) .and. z(1) == z(n)) then
       do i = n - 1, n - half_nz, -1
          j = j + 1
          ztmp(j) = z(i) ; rtmp(j) = r(i)
       end do
    else
       do i = n, n - half_nz + 1, -1
          j = j + 1
          ztmp(j) = z(i) ; rtmp(j) = r(i)
       end do
    end if

    if (debug) write(6,fmt='("rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("ztmp=",5f8.4)') ztmp
    if (debug) write(6,fmt='("Z_mag=",f8.4)') Z_mag

    spl = new_spline(ztmp, rtmp)
    r1 = splint(Z_mag, spl)
    call delete_spline(spl)

    ! find point near magnetic axis elevation on low field side
    ! Initialise i1 to some invalid value so we can check if an appropriate
    ! point is not found.
    i1 = invalidPoint

    ! Note there appears to be an assumption here about the order
    ! of z which may not be true
    do i = nz, n
       if (z(i) - Z_mag > 0.) then
          i1 = i - 1
          exit
       end if
    end do

    if (i1 == invalidPoint) then
       write(*, '("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")')
       call mp_abort("Could not find point near magnetic axis elevation on low field side in aminor of eeq.f90")
    end if

    !CMR, 28/10/08: modify this code to avoid duplicate points in 5 point spline
    !               to determine r2 on inboard mid-plane
    index = i1 - half_nz + k
    rtmp(1) = r(index)
    ztmp(1) = z(index)

    do i = 2, nz
       index = i1 - half_nz + i - 1 + k
       rtmp(i) = r(index)
       ztmp(i) = z(index)
       ! If we find a duplicate point then we increase the
       ! duplicate point counter in k, used to offset index and
       ! store the next value instead Note, here we assume we can't
       ! get three duplicate values. A better approach may be to
       ! instead use a while loop where we keep going until we've
       ! found enough points.
       if ((rtmp(i)-rtmp(i-1))**2+(ztmp(i)-ztmp(i-1))**2 < 1.0e-7) then
         k = k + 1
         if (debug) write(6,fmt='("a_minor: duplicate pt, set k=",i2)') k
         rtmp(i) = r(index + 1) ; ztmp(i) = z(index + 1)
       end if
    end do

    if (debug) write(6,fmt='("a_minor: rtmp=",5f8.4)') rtmp
    if (debug) write(6,fmt='("a_minor: ztmp=",5f8.4)') ztmp

    spl = new_spline(ztmp, rtmp)
    r2 = splint(Z_mag, spl)
    call delete_spline(spl)

    a = (r2 - r1) / 2.
    if (debug) write(6,*) "a_minor: r1,r2,a=",r1,r2,a
  end subroutine a_minor

  !> Put `theta` into the range \([-\pi, \pi]\)
  elemental real function mod2pi (theta)
    use constants, only: pi, twopi
    implicit none
    real, intent(in) :: theta
    mod2pi = theta
    if (theta > pi .or. theta < -pi) then
       mod2pi = modulo(mod2pi + pi, twopi) - pi
    end if
  end function mod2pi

  !> FIXME : Add documentation
  !>
  !> @note Could we use the sorting module here?
  pure subroutine sort(a, b, c, d)
    real, dimension(:), intent(in out) :: a, b, c, d
    real :: tmp
    integer :: i, j, jmax
    jmax = size(a)
    do j = 1, jmax
       do i = 1, jmax - j
          if (a(i+1) < a(i)) then
             tmp = a(i); a(i) = a(i+1); a(i+1) = tmp
             tmp = b(i); b(i) = b(i+1); b(i+1) = tmp
             tmp = c(i); c(i) = c(i+1); c(i+1) = tmp
             tmp = d(i); d(i) = d(i+1); d(i+1) = tmp
          end if
       end do
    end do
  end subroutine sort

  !> Returns the distance to the magnetic axis as a function of rp
  !> (the normalised poloidal flux variable) and theta.
  real function rfun(self, r, theta)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent(in) :: r, theta
    integer, parameter :: jmax = 5, kmax = 5
    real, parameter :: xerrsec = 1.e-9
    integer :: i, j, k, imax
    real :: fa, fb, bmult, thetroot, a, b, local_broot
    local_broot = self%bound(theta) ; thetroot = theta
    a = 0. ; b = local_broot ; bmult = 0.01
    fa = self%psi(a, thetroot) - r ; fb = self%psi(b, thetroot) - r
    ! need to bracket root, if we don't then try to search for a different b which does bracket
    if (fa * fb > 0.) then
       outer: do k = 1, kmax !Each k increases the distance we scan by (bmult)
          bmult = bmult * 2. ; imax = 10
          do j = 1, jmax !Each j increases the number of steps (imax)
             imax = imax * 2
             do i = 1, imax !Scan from local_broot to local_broot*(1 + bmult) in imax steps
                b = local_broot * (1 + bmult * float(i) / float(imax))
                fb = self%psi(b, thetroot) - r
                if (fa * fb <= 0.) exit outer !If a ; b brackets root skip ahead
             end do
             do i = 1, imax !Scan from local_broot to local_broot/(1 + bmult) in imax steps
                b = local_broot / (1 + bmult * float(i) / float(imax))
                fb = self%psi(b, thetroot) - r
                if (fa * fb <= 0.) exit outer  !If a ; b brackets root skip ahead
             end do
          end do
       end do outer !No promise that we have found a bracket here, but zbrent will report this.
    end if
    rfun = self%zbrent(a, b, r, thetroot, xerrsec)
  end function rfun

  !> FIXME : Add documentation
  real function zbrent(self, x1, x2, rootval, thetroot, tol)
    use mp, only: mp_abort
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: x1, x2, rootval, thetroot, tol
    real, parameter :: eps = 3.e-8, eps1 = 2.e-8
    integer, parameter :: itmax = 100
    real :: a, b, c, fa, fb, fc, d, e, tol1, xm, q, p, r, s
    integer :: iter
    zbrent = 0
    a = x1 ; b = x2
    fa = self%psi(a, thetroot) - rootval ; fb = self%psi(b, thetroot) - rootval
    if (fb * fa > 0.) then !If we're not bracketing the root
       ! check to see if fa or fb is close to zero already:
       if (abs(fa) < eps1) then
          zbrent = a ; return
       else if (abs(fb) < eps1) then
          zbrent = b ; return
       else
          write(*,*) a, b, fa, fb
          call mp_abort('root must be bracketed for zbrent.', .true.)
       end if
    end if

    c = 0. ; fc = fb
    do iter = 1, itmax
       if (fb * fc > 0.) then ! If root not bracketed by b and c swap a and c
          c = a ; fc = fa
          d = b - a ; e = d
       end if
       if (abs(fc) < abs(fb)) then !If c closer than b cyclic shuffle a<-b<-c<-a
          a = b ; b = c ; c = a
          fa = fb ; fb = fc ; fc = fa
       end if
       tol1 = 2. * eps * abs(b) + 0.5 * tol
       xm = 0.5 * (c - b) ! Bisect c and b
       if (abs(xm) <= tol1 .or. fb == 0.) then
          zbrent = b ; return
       end if
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then !e might not be set by this point?
          s = fb / fa
          if (a == c) then
             p = 2. * xm * s
             q = 1. - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.))
             q = (q - 1.) * (r - 1.) * (s - 1.)
          end if
          if (p > 0.) q = -q
          p = abs(p)
          if (2. * p  <  min(3. * xm * q - abs(tol1 * q), abs(e * q))) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a = b ; fa = fb
       if (abs(d)  >  tol1) then
          b = b + d
       else
          b = b + sign(tol1, xm)
       end if
       fb = self%psi(b, thetroot) - rootval
    end do
    call mp_abort('zbrent exceeding maximum iterations.',.true.)
  end function zbrent

  function new_flux_surface_type(geoType, Rmaj, R_geo, r, dr, aSurf, &
       sHorz, sVert, delm, deln, delmp, delnp, thm, thn, q, shat, nt, mMode, nMode, &
       n_mxh, c0_mxh, c_mxh, s_mxh, dc_mxh_dr, ds_mxh_dr) result(self)
    use optionals, only: get_option_with_default

    type(flux_surface_type) :: self
    integer, intent(in) :: geoType
    real, intent(in), optional :: Rmaj, R_geo, r, dr, aSurf, sHorz, sVert, delm, deln, delmp, delnp, thm, thn, q, shat, c0_mxh
    integer, intent(in), optional :: nt, mMode, nMode, n_mxh
    real, intent(in), optional, dimension(mxh_max_moments) :: c_mxh, dc_mxh_dr
    real, intent(in), optional, dimension(mxh_max_moments) :: s_mxh, ds_mxh_dr

    real, parameter :: invalid_value = -999.99
    integer, parameter :: invalid_int = -999

    self%geoType = geoType

    self%Rmaj = get_option_with_default(Rmaj, invalid_value)
    self%R_geo = get_option_with_default(R_geo, invalid_value)
    self%r = get_option_with_default(r, invalid_value)
    self%dr = get_option_with_default(dr, invalid_value)
    self%aSurf = get_option_with_default(aSurf, invalid_value)
    self%sHorz = get_option_with_default(sHorz, invalid_value)
    self%sVert = get_option_with_default(sVert, invalid_value)
    self%delm = get_option_with_default(delm, invalid_value)
    self%deln = get_option_with_default(deln, invalid_value)
    self%delmp = get_option_with_default(delmp, invalid_value)
    self%delnp = get_option_with_default(delnp, invalid_value)
    self%thm = get_option_with_default(thm, invalid_value)
    self%thn = get_option_with_default(thn, invalid_value)
    self%q = get_option_with_default(q, invalid_value)
    self%shat = get_option_with_default(shat, invalid_value)
    self%nt = get_option_with_default(nt, invalid_int)
    self%mMode = get_option_with_default(mMode, invalid_int)
    self%nMode = get_option_with_default(nMode, invalid_int)
    self%n_mxh = get_option_with_default(n_mxh, invalid_int)
    self%c0_mxh = get_option_with_default(c0_mxh, invalid_value)
    self%c_mxh = get_option_with_default(c_mxh, spread(invalid_value, 1, size(self%c_mxh)))
    self%s_mxh = get_option_with_default(s_mxh, spread(invalid_value, 1, size(self%s_mxh)))
    self%dc_mxh_dr = get_option_with_default(dc_mxh_dr, spread(invalid_value, 1, size(self%dc_mxh_dr)))
    self%ds_mxh_dr = get_option_with_default(ds_mxh_dr, spread(invalid_value, 1, size(self%ds_mxh_dr)))

  end function new_flux_surface_type
end module geo_utils
