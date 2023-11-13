!> A module to provide utilities used in geometry calculations.
!> Mainly to provide an abstract type used as a base for specific
!> geometry implementations.
module geo_utils
  implicit none

  public :: abstract_geo_type, abstract_eqfile_geo_type, abstract_eqfile_cart_geo_type
  public :: geo_input_type, flux_surface_type
  public :: mod2pi, sort, EQFILE_LENGTH

  integer, parameter :: EQFILE_LENGTH = 80

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
   contains
     procedure(initialise_interface), deferred :: initialise
     procedure(finalise_interface), deferred :: finalise
     procedure(calculate_gradients_interface), deferred, private :: calculate_gradients
     procedure(finish_setup_interface), deferred :: finish_setup
     procedure(rpos_interface), deferred :: rpos
     procedure(zpos_interface), deferred :: zpos
     procedure(invr_interface), deferred :: invr
     procedure(rcenter_interface), deferred :: rcenter
     procedure(diameter_interface), deferred :: diameter
     procedure(btori_interface), deferred :: btori
     procedure(dbtori_interface), deferred :: dbtori
     procedure(qfun_interface), deferred :: qfun
     procedure(pfun_interface), deferred :: pfun
     procedure(dpfun_interface), deferred :: dpfun
     procedure(betafun_interface), deferred :: betafun
     procedure(psi_interface), deferred :: psi
     procedure(rfun_interface), deferred :: rfun
     procedure(gradient_interface), deferred :: gradient
     procedure(bgradient_interface), deferred :: bgradient
     procedure(eqitem_interface), deferred, private :: eqitem
     procedure(alloc_arrays_interface), deferred, private :: alloc_arrays
     procedure(dealloc_arrays_interface), deferred, private :: dealloc_arrays
  end type abstract_geo_type

  type, abstract, extends(abstract_geo_type) :: abstract_eqfile_geo_type
   contains
     procedure :: finalise => finalise_eqfile
     procedure :: finish_setup => finish_setup_eqfile
     procedure(setup_splines_interface), deferred :: setup_splines
     procedure(delete_splines_interface), deferred :: delete_splines
  end type abstract_eqfile_geo_type

  type, abstract, extends(abstract_eqfile_geo_type) :: abstract_eqfile_cart_geo_type
   contains
     procedure(bound_interface), deferred :: bound
     procedure :: rfun
     procedure :: zbrent
  end type abstract_eqfile_cart_geo_type

  ! Define the procedure interfaces. Note we _could_ reuse
  ! several of these rather defining an interface for each routine.
  ! We don't do this at the moment for flexibility and clarity.
  interface
     real function rpos_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function rpos_interface

     real function zpos_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function zpos_interface

     real function invr_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function invr_interface

     real function rcenter_interface(self, rp)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: rp
     end function rcenter_interface

     real function diameter_interface(self, rp)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: rp
     end function diameter_interface

     real function btori_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function btori_interface

     real function dbtori_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function dbtori_interface

     real function qfun_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function qfun_interface

     real function pfun_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function pfun_interface

     real function dpfun_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function dpfun_interface

     real function betafun_interface(self, pbar)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: pbar
     end function betafun_interface

     real function psi_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function psi_interface

     real function rfun_interface(self, r, theta)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r, theta
     end function rfun_interface

     subroutine gradient_interface(self, rgrid, theta, grad, char, rp, nth_used, ntm)
       import
       class(abstract_geo_type), intent(in) :: self
       real, dimension(-ntm:), intent(in) :: rgrid, theta
       real, dimension(-ntm:, :), intent(out) :: grad
       character(len=1), intent(in) :: char
       real, intent(in) :: rp
       integer, intent(in) :: nth_used, ntm
     end subroutine gradient_interface

     subroutine bgradient_interface(self, rgrid, theta, grad, char, rp, nth_used, ntm)
       import
       class(abstract_geo_type), intent(in) :: self
       real, dimension(-ntm:), intent(in) :: rgrid, theta
       real, dimension(-ntm:, :), intent(out) :: grad
       character(len=1), intent(in) :: char
       real, intent(in) :: rp
       integer, intent(in) :: nth_used, ntm
     end subroutine bgradient_interface

     ! It might make sense to convert this to a function?
     subroutine eqitem_interface(self, r_in, theta_in, f, fstar, char)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r_in, theta_in
       real, dimension(:, :), intent(in) :: f
       real, intent(out) :: fstar
       character(len=1), intent(in) :: char
     end subroutine eqitem_interface

     pure subroutine alloc_arrays_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine alloc_arrays_interface

     pure subroutine dealloc_arrays_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine dealloc_arrays_interface

     subroutine initialise_interface(self, inputs, psi_0_out, psi_a_out, rmaj, &
          B_T0, avgrmid, nthg)
       import
       class(abstract_geo_type), intent(in out) :: self
       type(geo_input_type), intent(in) :: inputs
       real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
       integer, intent(out) :: nthg
     end subroutine initialise_interface

     subroutine finalise_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine finalise_interface

     subroutine calculate_gradients_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine calculate_gradients_interface

     subroutine finish_setup_interface(self)
       import
       class(abstract_geo_type), intent(in out) :: self
     end subroutine finish_setup_interface

     ! Additional interfaces for eqfile variant

     subroutine setup_splines_interface(self)
       import
       class(abstract_eqfile_geo_type), intent(in out) :: self
     end subroutine setup_splines_interface

     subroutine delete_splines_interface(self)
       import
       class(abstract_eqfile_geo_type), intent(in out) :: self
     end subroutine delete_splines_interface

     ! Additional interfaces for eqfile cartesian variant
     real function bound_interface(self, theta)
       import
       class(abstract_eqfile_cart_geo_type), intent(in) :: self
       real, intent(in) :: theta
     end function bound_interface

  end interface

contains

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
