!> A module to provide utilities used in geometry calculations.
!> Mainly to provide an abstract type used as a base for specific
!> geometry implementations.
module geo_utils
  use splines, only: spline, periodic_spline
  implicit none

  public :: abstract_geo_type, abstract_eqfile_geo_type, abstract_eqfile_cart_geo_type
  public :: geo_input_type, flux_surface_type, mod2pi, sort

  integer, parameter :: EQFILE_LENGTH = 800

  !> Helper type to wrap up some flux surface specific variables for
  !> analytic equilibrium type.
  type :: flux_surface_type
     !> Flux surface major radius, \(R_{0N}\), input variable: `Rmaj`
     real :: Rmaj
     !> Location of reference magnetic field, \(R_{geoN}\), input variable: `R_geo`
     real :: Rgeo
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
     !> Number of points in theta, input variable: `ntheta(nperiod - 1/2)`
     integer :: nt
     !> Geometry specification type, input variable: `geoType`
     integer :: geoType
     !> The first poloidal mode number, input variable: `mMode`
     integer :: mMode
     !> The second poloidal mode number, input variable: `nMode`
     integer :: nMode
  end type flux_surface_type

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
     logical :: initialised = .false.
     !> Typically the number of radial and theta grid points. Not used
     !> for cartesian like types (e.g. eeq)
     integer :: nr, nt
     real :: psi_0, psi_a, B_T, beta_0, R_mag, Z_mag, aminor
     !> 2D map of co-ordinates to R or Z
     real, allocatable, dimension (:,:)   :: R_psi, Z_psi
     ! The following quantities represent changes between adjacent grid
     ! points. The first two dimensions are minor radius and theta,
     ! respectively. The last dimension is the difference along either
     ! the (minor) radial or theta directions, respectively
     !> Major radius differences on theta, minor radius grids
     real, allocatable, dimension (:,:,:) :: drm
     !> Vertical differences on theta, minor radius grids
     real, allocatable, dimension (:,:,:) :: dzm
     !> Minor radius differences on theta, minor radius grids, note: `dpm(:,:,2) == 0.0`
     real, allocatable, dimension (:,:,:) :: dpm
     !> Theta differences on theta, minor radius grids, note: `dtm(:,:,1) == 0.0`
     real, allocatable, dimension (:,:,:) :: dtm
     !> Minor radius gradient in cylindrical coordinates \((\hat{e}_R,
     !> \hat{e}_Z, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dpcart
     !> Theta gradient in cylindrical coordinates \((\hat{e}_R, \hat{e}_Z, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dtcart
     !> Minor radius gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dpbish
     !> Theta gradient in Bishop coordinates \((\hat{e}_{r_\psi}, \hat{e}_l, \hat{e}_\zeta)\)
     real, allocatable, dimension (:,:,:) :: dtbish
     !> Equilibrium psi (not normalised) and pressure grids
     real, allocatable, dimension(:) :: eqpsi, pressure
   contains
     procedure(initialise_interface), deferred :: initialise
     procedure(noargs_interface), deferred :: finalise, calculate_gradients, finish_setup
     procedure(rtheta_interface), deferred :: rpos, zpos, invr, psi, rfun
     procedure(rp_interface), deferred :: rcenter, diameter
     procedure(pbar_interface), deferred :: btori, dbtori, qfun, pfun, dpfun, betafun
     !> Calculate the derivative of rp w.r.t. R and Z. Parameters
     !> are:
     !>   - rgrid: Submodule radial grid (r_circ for EFIT, rp otherwise)
     !>   - theta: grid of theta
     !>   - gradf(:,1): deriv w.r.t. R as fn of theta
     !>   - gradf(:,2): deriv w.r.t  Z as fn of theta
     !>   - char: if char = 'R', return the gradient of pressure instead
     !>   - rp: value of rp on the flux surface where we want the grad
     !>   - nth: number of theta points
     !>   - ntgrid: lower index of theta array is -ntgrid
     procedure(gradient_interface), deferred :: gradient
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
     procedure(gradient_interface), deferred :: bgradient
     procedure(eqitem_interface), deferred :: eqitem
     procedure(alloc_arrays_interface), deferred :: alloc_arrays, dealloc_arrays
     procedure :: eqdcart, eqdbish
  end type abstract_geo_type

  type, abstract, extends(abstract_geo_type) :: abstract_eqfile_geo_type
     character(len = EQFILE_LENGTH) :: filename
     real, allocatable, dimension(:) :: psi_bar, fp, qsf, beta, diam, rc
     type(spline) :: btori_spl, qsf_spl, pressure_spl, beta_spl, diam_spl, rcenter_spl
   contains
     procedure :: finalise => finalise_eqfile
     procedure :: finish_setup => finish_setup_eqfile
     procedure(splines_interface), deferred :: setup_splines, delete_splines
     procedure :: btori => btori_eqfile, dbtori => dbtori_eqfile
     procedure :: qfun => qfun_eqfile, betafun => betafun_eqfile
     procedure :: pfun => pfun_eqfile, dpfun => dpfun_eqfile
     procedure :: diameter => diameter_eqfile, rcenter => rcenter_eqfile
     procedure :: psi => psi_eqfile, rfun => psi_eqfile
     procedure :: invR => invR_eqfile, rpos => rpos_eqfile, zpos => zpos_eqfile
     procedure :: alloc_common_arrays, dealloc_common_arrays
     procedure :: alloc_arrays => alloc_arrays_eqfile
     procedure :: dealloc_arrays => dealloc_arrays_eqfile
     procedure(alloc_eqfile_arrays_interface), deferred :: alloc_special_arrays
     procedure(alloc_eqfile_arrays_interface), deferred :: dealloc_special_arrays
  end type abstract_eqfile_geo_type

  type, abstract, extends(abstract_eqfile_geo_type) :: abstract_eqfile_cart_geo_type
     integer :: nw, nh, nbbbs
     real :: eq_dR, eq_dZ
     real, allocatable, dimension(:) :: thetab, r_bound, eq_R, eq_Z
     real, allocatable, dimension (:,:) :: eq_psi
     type(periodic_spline) :: rbound_spl
   contains
     procedure :: bound => bound_eqfile_cart, rfun, zbrent
     procedure :: psi => psi_eqfile_cart, rpos => rpos_eqfile_cart, zpos => zpos_eqfile_cart
  end type abstract_eqfile_cart_geo_type

  ! Define the procedure interfaces.
  interface
     real function rtheta_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function rtheta_interface

     real function rp_interface(self, rp)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: rp
     end function rp_interface

     real function pbar_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function pbar_interface

     subroutine gradient_interface(self, rgrid, theta, grad, char, rp, nth_used, ntm)
       import
       class(abstract_geo_type), intent(in) :: self
       real, dimension(-ntm:), intent(in) :: rgrid, theta
       real, dimension(-ntm:, :), intent(out) :: grad
       character(len=1), intent(in) :: char
       real, intent(in) :: rp
       integer, intent(in) :: nth_used, ntm
     end subroutine gradient_interface

     ! It might make sense to convert this to a function?
     real function eqitem_interface(self, r_in, theta_in, f, char) result(fstar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r_in, theta_in
       real, dimension(:, :), intent(in) :: f
       character(len=1), intent(in) :: char
     end function eqitem_interface

     pure subroutine alloc_arrays_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine alloc_arrays_interface

     pure subroutine alloc_eqfile_arrays_interface(self)
       import
       class(abstract_eqfile_geo_type), intent(in out) :: self
     end subroutine alloc_eqfile_arrays_interface

     subroutine initialise_interface(self, inputs, psi_0_out, psi_a_out, rmaj, &
          B_T0, avgrmid, nthg)
       import
       class(abstract_geo_type), intent(in out) :: self
       type(geo_input_type), intent(in) :: inputs
       real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
       integer, intent(out) :: nthg
     end subroutine initialise_interface

     subroutine noargs_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine noargs_interface

     ! Additional interfaces for eqfile variant
     subroutine splines_interface(self)
       import
       class(abstract_eqfile_geo_type), intent(in out) :: self
     end subroutine splines_interface
  end interface

contains

  !> Converts derivatives w.r.t. (psi_index,theta_index) to derivatives
  !> w.r.t. (R,Z).
  !> Note that `dfcart(1, :, :)` is not valid
  pure subroutine eqdcart(self, dfm, dfcart)
    implicit none
    class(abstract_geo_type), intent(in) :: self
    !> dfm(:,:,1) is deriv w.r.t. psi_index (i.e. (df/dpsi)_theta * delta psi);
    !> dfm(:,:,2) is deriv w.r.t. theta_index
    real, dimension (:,:,:), intent(in)  :: dfm
    !> dfcart(:,:,1) is deriv w.r.t. R;
    !> dfcart(:,:,2) is deriv w.r.t. Z
    real, dimension (:,:,:), intent(out) :: dfcart
    real, dimension (size(dfm, 1), size(dfm, 2)) :: denom
    integer :: i, j
    denom = self%drm(:,:,1) * self%dzm(:,:,2) - self%drm(:,:,2) * self%dzm(:,:,1)
    dfcart(:,:,1) =   dfm(:,:,1) * self%dzm(:,:,2) - self%dzm(:,:,1) * dfm(:,:,2)
    dfcart(:,:,2) = - dfm(:,:,1) * self%drm(:,:,2) + self%drm(:,:,1) * dfm(:,:,2)

    ! Inner loop skips first point because it's likely exactly on the
    ! magnetic axis, which means the Jacobian is zero. Do we need
    ! need a better way to either calculate dfcart there, or to ensure
    ! we ignore it in other places?
    do j = 1, size(dfm, 2)
       do i = 2, size(dfm, 1)
          dfcart(i,j,:) = dfcart(i,j,:) / denom(i,j)
       end do
    end do
  end subroutine eqdcart

  !> Convert gradients of a function f w.r.t. R,Z into bishop form
  !> Note that `dbish(1, :, :)` is not valid
  pure subroutine eqdbish(self, dcart, dbish)
    implicit none
    class(abstract_geo_type), intent(in) :: self
    !> dcart(:,:,1) is gradient of f w.r.t. R;
    !> dcart(:,:,2) is gradient of f w.r.t. Z
    real, dimension(:, :, :), intent (in) :: dcart
    !> dbish(:,:,1) is set to (df/dR dpsi/dR + df/dZ dpsi/dZ)/|grad psi|;
    !> dbish(:,:,2) is set to (-df/dR dpsi/dZ + df/dZ dpsi/dR)/|grad psi|.
    !> Note that in the special case where f=psi:
    !> dbish(:,:,1) is |grad psi|;
    !> dbish(:,:,2) is 0
    real, dimension(:, :, :), intent(out) :: dbish
    real, dimension(size(dcart, 1), size(dcart, 2)) :: denom
    integer :: i, j
    ! denom is |grad psi|
    denom = sqrt(self%dpcart(:,:,1)**2 + self%dpcart(:,:,2)**2)
    dbish(:,:,1) = dcart(:,:,1) * self%dpcart(:,:,1) + dcart(:,:,2) * self%dpcart(:,:,2)
    dbish(:,:,2) =-dcart(:,:,1) * self%dpcart(:,:,2) + dcart(:,:,2) * self%dpcart(:,:,1)

    do j = 1, size(dcart, 2)
       do i = 2, size(dcart, 1)
          dbish(i,j,:) = dbish(i,j,:)/ denom(i,j)
       end do
    end do
  end subroutine eqdbish

  pure subroutine alloc_arrays_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    call self%alloc_common_arrays
    call self%alloc_special_arrays
  end subroutine alloc_arrays_eqfile

  pure subroutine dealloc_arrays_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    call self%dealloc_common_arrays
    call self%dealloc_special_arrays
  end subroutine dealloc_arrays_eqfile

  pure subroutine alloc_common_arrays(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    integer :: nr, nt
    nr = self%nr ; nt = self%nt
    allocate(self%R_psi(nr, nt), self%Z_psi(nr, nt), self%drm(nr, nt, 2), self%dzm(nr, nt, 2))
    allocate(self%dpm(nr, nt, 2), self%dtm(nr, nt, 2), self%dpcart(nr, nt, 2))
    allocate(self%dtcart(nr, nt, 2), self%dpbish(nr, nt, 2), self%dtbish(nr, nt, 2))
    allocate(self%eqpsi(nr), self%pressure(nr))
  end subroutine alloc_common_arrays

  pure subroutine dealloc_common_arrays(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    deallocate(self%R_psi, self%Z_psi, self%drm, self%dzm, self%dpm, self%dtm, self%dpcart)
    deallocate(self%dtcart, self%dpbish, self%dtbish, self%eqpsi, self%pressure)
  end subroutine dealloc_common_arrays

  subroutine finalise_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    self%initialised = .false.
    call self%dealloc_arrays
    call self%delete_splines
  end subroutine finalise_eqfile

  subroutine finish_setup_eqfile(self)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    call self%calculate_gradients
    call self%setup_splines
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
    dbtori = dsplint(pbar, self%btori_spl) / (self%psi_a-self%psi_0)
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
    dpfun = 0.5 * self%beta_0 * dsplint(pbar, self%pressure_spl) / (self%psi_a-self%psi_0)
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

  real function psi_eqfile (self, r, theta) result(psi)
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
       return
    end if
    betafun = splint(pbar, self%beta_spl)
  end function betafun_eqfile

  real function invR_eqfile (self, r, theta) result(invR)
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    invR = 1. / self%rpos(r, theta)
  end function invR_eqfile

  real function rpos_eqfile (self, r, theta) result(rpos)
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    rpos = self%eqitem(r, th, self%R_psi, 'R')
  end function rpos_eqfile

  real function zpos_eqfile (self, r, theta) result(zpos)
    implicit none
    class(abstract_eqfile_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    real :: th
    th = mod2pi(theta)
    zpos = self%eqitem(r, th, self%Z_psi, 'Z')
  end function zpos_eqfile

  !> Find the value of the major radius on the plasma boundary
  !> at given geometric theta.
  real function bound_eqfile_cart(self, theta) result(bound)
    use splines, only: periodic_splint
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent(in) :: theta
    bound = periodic_splint(theta, self%rbound_spl)
  end function bound_eqfile_cart

  real function rpos_eqfile_cart (self, r, theta) result(rpos)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent (in) :: r, theta
    rpos = self%R_mag + r * cos(theta)
  end function rpos_eqfile_cart

  real function zpos_eqfile_cart (self, r, theta) result(zpos)
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
       psi = self%eqitem(r, theta, self%eq_psi, '?')
    end if
  end function psi_eqfile_cart

  !> Put `theta` into the range \([-\pi, \pi]\)
  elemental real function mod2pi (theta)
    use constants, only: pi
    implicit none
    real, intent(in) :: theta
    real :: th
    logical :: out
    if (theta <= pi .and. theta >= -pi) then
       mod2pi = theta
    else
       th = theta
       out = .true.
       do while(out)
          if (th > pi) th = th - 2.*pi
          if (th <-pi) th = th + 2.*pi
          if (th <= pi .and. th >= -pi) out = .false.
       end do
       mod2pi = th
    end if
  end function mod2pi

  !> FIXME : Add documentation
  !>
  !> @note Could we use the sorting module here?
  subroutine sort(a, b, c, d)
    real, dimension(:), intent(in out) :: a, b, c, d
    real tmp
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

  !> Returns the distance to the magnetic axis
  !> as a function of rp (the normalised poloidal flux variable)
  !> and theta.
  real function rfun(self, r, theta)
    implicit none
    class(abstract_eqfile_cart_geo_type), intent(in) :: self
    real, intent(in) :: r, theta
    integer :: i, j, k, imax, jmax, kmax
    real :: fa, fb, fbroot, bmult, rootval, thetroot
    real :: a, b, xerrsec, local_broot

    local_broot = self%bound(theta)
    thetroot = theta
    rootval = r
    a = 0.
    b = local_broot

    fa = self%psi(a, thetroot) - rootval
    fbroot = self%psi(b, thetroot) - rootval
    fb = fbroot
    if (fa*fbroot <= 0.) goto 100

    ! need to bracket root.
    bmult = 0.01
    jmax = 5
    kmax = 5

    do k = 1, kmax
       bmult = bmult * 2.
       imax = 10
       do j = 1, jmax
          imax = imax * 2
          do i = 1, imax
             b = local_broot * (1 + bmult * float(i) / float(imax))
             fb = self%psi(b, thetroot) - rootval
             if (fa*fb <= 0.) goto 100
          end do
          do i = 1, imax
             b = local_broot / (1 + bmult * float(i) / float(imax))
             fb = self%psi(b, thetroot) - rootval
             if (fa*fb <= 0.) goto 100
          end do
       end do
    end do
    write(*,*) 'yuk: Did not find bracket of root in rfun'
100 continue
    xerrsec = 1.e-9
    rfun = self%zbrent(a, b, rootval, thetroot, xerrsec)
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

    a = x1
    b = x2
    c = 0.
    fa = self%psi(a, thetroot) - rootval
    fb = self%psi(b, thetroot) - rootval
    if (fb * fa > 0.) then
       ! check to see if fa or fb is close to zero already:
       if (abs(fa) < eps1) then
          zbrent = a
          return
       else if (abs(fb) < eps1) then
          zbrent = b
          return
       end if
       write(*,*) a, b, fa, fb
       call mp_abort('root must be bracketed for zbrent.', .true.)
    end if
    fc = fb

    do iter = 1, itmax
       if (fb * fc > 0.) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2. * eps * abs(b) + 0.5 * tol
       xm = 0.5 * (c - b)
       if (abs(xm) <= tol1 .or. fb == 0.) then
          zbrent = b
          return
       end if
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
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
       a = b
       fa = fb
       if (abs(d)  >  tol1) then
          b = b + d
       else
          b = b + sign(tol1, xm)
       end if
       fb = self%psi(b, thetroot) - rootval
    end do
    call mp_abort('zbrent exceeding maximum iterations.',.true.)
    zbrent = b
  end function zbrent

end module geo_utils
