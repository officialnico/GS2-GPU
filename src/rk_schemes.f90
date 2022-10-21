!> A small helper module for storing and providing RK schemes
module rk_schemes
  implicit none

  private

  public :: rk_scheme, get_rk_scheme, get_rk_scheme_by_id
  public :: get_rk_schemes_as_text_options
  public :: rk_dopri, rk_f45, rk_cashkarp
  public :: rk_heun, rk_ralston, rk_midpoint
  public :: rk_bs23, rk_bs45

  !> Internal named type to represent scheme ids
  type rk_scheme_id
     private
     integer :: id
     character(len=20) :: name
  end type rk_scheme_id

  !> Public type storing the data for a particular RK scheme
  type :: rk_scheme
     integer :: number_of_stages
     integer :: order
     character(len=20) :: name
     real, dimension(:), allocatable :: lower_order_coeffs
     real, dimension(:), allocatable :: high_order_coeffs
     real, dimension(:), allocatable :: time_points
     real, dimension(:, :), allocatable :: coeffs
     logical :: follow_high_order
   contains
     procedure :: validate => rk_scheme_validate
  end type rk_scheme

  !> Dormand-Prince 45 method
  type(rk_scheme_id), parameter :: rk_dopri = rk_scheme_id(id=1, name='dopri')

  !> Runge-Kutta-Fehlberg 45 method
  type(rk_scheme_id), parameter :: rk_f45 = rk_scheme_id(id=2, name='rkf45')

  !> Cash-Karp 45 method
  type(rk_scheme_id), parameter :: rk_cashkarp = rk_scheme_id(id=3, name='cashkarp')

  !> Heun-Euler 12 method
  type(rk_scheme_id), parameter :: rk_heun = rk_scheme_id(id=4, name='heun')

  !> Ralston-Euler 12 method
  type(rk_scheme_id), parameter :: rk_ralston = rk_scheme_id(id=5, name='ralston')

  !> Midpoint-Euler 12 method
  type(rk_scheme_id), parameter :: rk_midpoint = rk_scheme_id(id=6, name='midpoint')

  !> Bogacki-Shampine 23 method
  type(rk_scheme_id), parameter :: rk_bs23 = rk_scheme_id(id=7, name='bs23')

  !> Bogacki-Shampine 45 method
  type(rk_scheme_id), parameter :: rk_bs45 = rk_scheme_id(id=8, name='bs45')

contains

  !> Returns an array of text_option values representing the known rk schemes, mapping
  !> names to integer ids
  pure function get_rk_schemes_as_text_options() result(options)
    use text_options, only: text_option
    implicit none
    type(text_option), dimension(:), allocatable :: options

    options = [ &
         text_option(rk_dopri%name, rk_dopri%id), &
         text_option(rk_f45%name, rk_f45%id), &
         text_option(rk_cashkarp%name, rk_cashkarp%id), &
         text_option(rk_heun%name, rk_heun%id), &
         text_option(rk_ralston%name, rk_ralston%id), &
         text_option(rk_midpoint%name, rk_midpoint%id), &
         text_option(rk_bs23%name, rk_bs23%id), &
         text_option(rk_bs45%name, rk_bs45%id), &
         text_option('default', rk_heun%id) &
         ]
  end function get_rk_schemes_as_text_options

  !> Returns an rk_scheme instance describing the scheme with id matching
  !> the passed id value.
  type(rk_scheme) function get_rk_scheme_by_id(id) result(scheme)
    use mp, only: mp_abort
    implicit none
    integer, intent(in) :: id
    select case(id)
    case(rk_dopri%id)
       scheme = get_rk_scheme(rk_dopri)
    case(rk_f45%id)
       scheme = get_rk_scheme(rk_f45)
    case(rk_cashkarp%id)
       scheme = get_rk_scheme(rk_cashkarp)
    case(rk_heun%id)
       scheme = get_rk_scheme(rk_heun)
    case(rk_ralston%id)
       scheme = get_rk_scheme(rk_ralston)
    case(rk_midpoint%id)
       scheme = get_rk_scheme(rk_midpoint)
    case(rk_bs23%id)
       scheme = get_rk_scheme(rk_bs23)
    case(rk_bs45%id)
       scheme = get_rk_scheme(rk_bs45)
    case default
       call mp_abort("Invalid RK scheme id passed to get_rk_scheme_by_id", .true.)
    end select
  end function get_rk_scheme_by_id

  !> Returns an rk_scheme instance describing the scheme represented by
  !> the passed variant type.
  type(rk_scheme) function get_rk_scheme(variant) result(scheme)
    use mp, only: mp_abort, proc0
    implicit none
    type(rk_scheme_id), intent(in) :: variant
    integer :: validity

    scheme%name = variant%name

    select case(variant%id)
    case(rk_dopri%id)
       scheme%number_of_stages = 7
       scheme%order = 5
       scheme%lower_order_coeffs = [5179/57600.0, 0.0,7571/16695.0, 393/640.0, -92097/339200.0, 187/2100.0, 1/40.0]
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

    case(rk_f45%id)
       scheme%number_of_stages = 6
       scheme%order = 5
       scheme%lower_order_coeffs = [25/216.0, 0.0, 1408/2565.0, 2197/4104.0, -1/5.0, 0.0]
       scheme%high_order_coeffs = [16/135.0, 0.0, 6656/12825.0, 28561/56430.0, -9/50.0, 2/55.0]
       scheme%time_points = [0.0, 1/4.0, 3/8.0, 12/13.0, 1.0, 1/2.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/4.0]
       scheme%coeffs(:2, 3) = [3/32.0, 9/32.0]
       scheme%coeffs(:3, 4) = [1932/2197.0, -7200/2197.0, 7296/2197.0]
       scheme%coeffs(:4, 5) = [439/216.0, -8.0, 3680/513.0, -845/4104.0]
       scheme%coeffs(:5, 6) = [-8/27.0, 2.0, -3544/2565.0, 1859/4104.0, -11/40.0]

    case(rk_cashkarp%id)
       scheme%number_of_stages = 6
       scheme%order = 5
       scheme%lower_order_coeffs = [2825/27648.0, 0.0, 18575/48384.0, 13525/55296.0, 277/14336.0, 1/4.0]
       scheme%high_order_coeffs = [37/378.0, 0.0, 250/621.0, 125/594.0, 0.0, 512/1771.0]
       scheme%time_points = [0.0, 1/5.0, 3/10.0, 3/5.0, 1.0, 7/8.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/5.0]
       scheme%coeffs(:2, 3) = [3/40.0, 9/40.0]
       scheme%coeffs(:3, 4) = [3/10.0, -9/10.0, 6/5.0]
       scheme%coeffs(:4, 5) = [-11/54.0, 5/2.0, -70/27.0, 35/27.0]
       scheme%coeffs(:5, 6) = [1631/55296.0, 175/512.0, 575/13824.0, 44275/110592.0, 253/4096.0]

    case(rk_heun%id)
       scheme%number_of_stages = 2
       scheme%order = 2
       scheme%lower_order_coeffs = [1.0, 0.0]
       scheme%high_order_coeffs = [1/2.0, 1/2.0]
       scheme%time_points = [0.0, 1.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1.0]

    case(rk_ralston%id)
       scheme%number_of_stages = 2
       scheme%order = 2
       scheme%lower_order_coeffs = [1.0, 0.0]
       scheme%high_order_coeffs = [1/4.0, 3/4.0]
       scheme%time_points = [0.0, 2/3.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [2/3.0]

    case(rk_midpoint%id)
       scheme%number_of_stages = 2
       scheme%order = 2
       scheme%lower_order_coeffs = [1.0, 0.0]
       scheme%high_order_coeffs = [0.0, 1.0]
       scheme%time_points = [0.0, 1/2.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/2.0]

    case(rk_bs23%id)
       scheme%number_of_stages = 4
       scheme%order = 3
       scheme%lower_order_coeffs = [7/24.0, 1/4.0, 1/3.0, 1/8.0]
       scheme%high_order_coeffs = [2/9.0, 1/3.0, 4/9.0, 0.0]
       scheme%time_points = [0.0, 1/2.0, 3/4.0, 1.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/2.0]
       scheme%coeffs(:2, 3) = [0.0, 3/4.0]
       scheme%coeffs(:3, 4) = [2/9.0, 1/3.0, 4/9.0]

    case(rk_bs45%id)
       scheme%number_of_stages = 8
       scheme%order = 5
       scheme%lower_order_coeffs = [6059/80640.0, 0.0, 8559189/30983680.0, 26411/124800.0, &
            -927/89600.0, 443/1197.0, 7267/94080.0, 0.0]
       scheme%high_order_coeffs = [587/8064.0, 0.0, 4440339/15491840.0, 24353/124800.0, &
            387/44800.0, 2152/5985.0, 7267/94080.0, 0.0]
       scheme%time_points = [0.0, 1/6.0, 2/9.0, 3/7.0, 2/3.0, 3/4.0, 1.0, 1.0]
       scheme%follow_high_order = .true.
       allocate(scheme%coeffs(scheme%number_of_stages - 1, scheme%number_of_stages))
       scheme%coeffs = 0.0
       scheme%coeffs(:1, 2) = [1/6.0]
       scheme%coeffs(:2, 3) = [2/27.0, 4/27.0]
       scheme%coeffs(:3, 4) = [183/1372.0, -162/343.0, 1053/1372.0]
       scheme%coeffs(:4, 5) = [68/297.0, -4/11.0, 42/143.0, 1960/3861.0]
       scheme%coeffs(:5, 6) = [597/22528.0, 81/352.0, 63099/585728.0, 58653/366080.0, &
            4617/20480.0]
       scheme%coeffs(:6, 7) = [174197/959244.0, -30942/79937.0, 8152137/19744439.0, &
            666106/1039181.0, -29421/29068.0, 482048/414219.0]
       scheme%coeffs(:7, 8) = [587/8064.0, 0.0, 4440339/15491840.0, 24353/124800.0, &
            387/44800.0, 2152/5985.0, 7267/94080.0]
    case default
       call mp_abort("Invalid RK scheme passed to get_rk_scheme", .true.)
    end select

    validity = scheme%validate()
    if (proc0 .and. (validity /= 0)) then
       print*,'The RK scheme failed validity check with error code:',validity,'and scheme',scheme%name
       call mp_abort("RK scheme failed validity check", .true.)
    end if

  end function get_rk_scheme

  !> Perform a small check of the consistency. Returns an error
  !> code built from bitwise flags to give information about which
  !> checks, if any, failed. Success is repesented by zero.
  integer pure function rk_scheme_validate(self) result(valid)
    implicit none
    class(rk_scheme), intent(in) :: self
    integer :: i
    real, parameter :: tolerance = 10*epsilon(1.0)
    valid = 0
    if (.not.(abs(sum(self%lower_order_coeffs) - 1.0) < tolerance)) valid = valid + 1
    if (.not.(abs(sum(self%high_order_coeffs) - 1.0) < tolerance)) valid = valid + 2
    do i = 1, self%number_of_stages
       if (.not.(abs(sum(self%coeffs(:, i)) - self%time_points(i)) < tolerance)) valid = valid + (2**(i+1))
    end do
  end function rk_scheme_validate

end module rk_schemes
