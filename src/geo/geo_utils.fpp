!> A module to provide utilities used in geometry calculations.
!> Mainly to provide an abstract type used as a base for specific
!> geometry implementations.
module geo_utils
  implicit none

  public :: abstract_geo_type, abstract_eqfile_geo_type
  public :: mod2pi, EQFILE_LENGTH

  integer, parameter :: EQFILE_LENGTH = 80

  type, abstract :: abstract_geo_type
     logical :: initialised = .false.
   contains
     generic :: initialise => initialise_eqfile, initialise_analytic
     procedure(initialise_eqfile_interface), deferred :: initialise_eqfile
     procedure(initialise_analytic_interface), deferred :: initialise_analytic
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
     procedure(gradient_interface), deferred :: gradient
     procedure(bgradient_interface), deferred :: bgradient
     procedure(eqitem_interface), deferred, private :: eqitem
  end type abstract_geo_type

  type, abstract, extends(abstract_geo_type) :: abstract_eqfile_geo_type
   contains
     procedure :: initialise_analytic => initialise_analytic_null
  end type abstract_eqfile_geo_type

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

     real function psi_interface(self, r)
       import
       class(abstract_geo_type), intent(in) :: self
       real, intent(in) :: r
     end function psi_interface

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

     subroutine initialise_eqfile_interface(self, eqfile, psi_0_out, psi_a_out, rmaj, &
          B_T0, avgrmid, initeq, nthg)
       import
       class(abstract_geo_type), intent(in out) :: self
       character(len = EQFILE_LENGTH), intent(in) :: eqfile
       real, intent(out) :: psi_0_out, psi_a_out, rmaj, B_T0, avgrmid
       integer, intent(in) :: initeq
       integer, intent(out) :: nthg
     end subroutine initialise_eqfile_interface

     subroutine initialise_analytic_interface(self, geoType, Rmaj, Rgeo, r, dr, aSurf, &
          sHorz, sVert, mMode, nMode, delm, deln, delmp, delnp, thm, thn, q, &
          shat, nt_used, thShift)
       import
       class(abstract_geo_type), intent(in out) :: self
       real, intent(in) :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, delm, deln, &
            delmp, delnp, thm, thn, q, shat
       real, intent(in) :: thShift
       integer, intent(in) :: nt_used, geoType, mMode, nMode
     end subroutine initialise_analytic_interface

  end interface

contains

  subroutine initialise_analytic_null(self, geoType, Rmaj, Rgeo, r, dr, aSurf, &
       sHorz, sVert, mMode, nMode, delm, deln, delmp, delnp, thm, thn, q, &
       shat, nt_used, thShift)
    class(abstract_eqfile_geo_type), intent(in out) :: self
    real, intent(in) :: Rmaj, Rgeo, r, dr, aSurf, sHorz, sVert, delm, deln, &
         delmp, delnp, thm, thn, q, shat
    real, intent(in) :: thShift
    integer, intent(in) :: nt_used, geoType, mMode, nMode
    write(*,*) 'Not implemented'
  end subroutine initialise_analytic_null

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

end module geo_utils
