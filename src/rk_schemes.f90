!> A small helper module for storing and providing RK schemes
module rk_schemes
  implicit none

  private

  public :: rk_scheme, get_rk_scheme
  public :: rk_dopri

  !> Internal named type to represent scheme ids
  type rk_scheme_id
     integer, private :: id
  end type rk_scheme_id

  !> Public type storing the data for a particular RK scheme
  type :: rk_scheme
     integer :: number_of_stages
     integer :: order
     real, dimension(:), allocatable :: lower_order_coeffs
     real, dimension(:), allocatable :: high_order_coeffs
     real, dimension(:), allocatable :: time_points
     real, dimension(:, :), allocatable :: coeffs
     logical :: follow_high_order
  end type rk_scheme

  !> Dormand-Prince 45 method
  type(rk_scheme_id), parameter :: rk_dopri = rk_scheme_id(id=1)

contains

  !> Returns an rk_scheme instance describing the scheme represented by
  !> the passed variant type.
  pure type(rk_scheme) function get_rk_scheme(variant) result(scheme)
    implicit none
    type(rk_scheme_id), intent(in) :: variant

    select case(variant%id)
    case(rk_dopri%id)
       scheme%number_of_stages = 7
       scheme%order = 5
       scheme%lower_order_coeffs = [5179/57600.0, 0.0,7571/16695.0, 393/640.0, -92097/339200.0, 187/21000.0, 1/40.0]
       scheme%high_order_coeffs = [35/384.0, 0.0, 500/1113.0, 125/192.0, -2187/6784.0, 11/84.0, 0.0]
       scheme%time_points = [0.0, 1/5.0, 3/10.0, 4/5.0, 8/9.0, 1.0, 1.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/5.0]
       scheme%coeffs(:2, 3) = [3/40.0, 9/40.0]
       scheme%coeffs(:3, 4) = [44/45.0, -56/15.0, 32/9.0]
       scheme%coeffs(:4, 5) = [19372/6561.0, -25360/2187.0, 64448/6561.0, -212/729.0]
       scheme%coeffs(:5, 6) = [9017/3168.0, -355/33.0, 46732/5247.0, 49/176.0, -5103/18656.0]
       scheme%coeffs(:6, 7) = [35/384.0, 0.0, 500/1113.0, 125/192.0, -2187/6784.0, 11/84.0, 0.0]
    end select

  end function get_rk_scheme

end module rk_schemes
