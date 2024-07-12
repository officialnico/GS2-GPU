!> Set up ranges of kx and ky for linear runs.
module kt_grids_range
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN

  implicit none

  private

  public :: init_kt_grids_range, range_get_sizes, range_get_grids
  public :: check_kt_grids_range
  public :: wnml_kt_grids_range
  public :: read_parameters_range
  public :: finish_parameters_range

  public :: kt_grids_range_config_type
  public :: set_kt_grids_range_config
  public :: get_kt_grids_range_config

  integer :: naky, ntheta0, nn0, n0_min, n0_max
  real :: aky_min, aky_max, theta0_min, theta0_max
  real :: akx_min, akx_max, rhostar_range
  character(20) :: kyspacing_option
  integer :: kyspacingopt_switch

  logical :: parameters_read = .false.
  logical :: initialized = .false.

  !Note if we ever want to offer different spacing for theta0 we could
  !reuse these flags (rename to spacingopt_...).
  integer, parameter :: kyspacingopt_linear=1, kyspacingopt_exp=2

  !> Used to represent the input configuration of kt_grids_range
  type, extends(abstract_config_type) :: kt_grids_range_config_type
     ! namelist : kt_grids_range_parameters
     ! indexed : false
     !> Max kx for periodic finite kx ballooning space runs with
     !> \(\hat{s}=0\).
     real :: akx_max = 0.0
     !> Min kx for periodic finite kx ballooning space runs with
     !> \(\hat{s}=0\).
     real :: akx_min = 0.0
     !> Upper limit of \(k_y \rho\) range. Should set to something other
     !> than zero.
     real :: aky_max = 0.0
     !> Lower limit of \(k_y \rho\) range. Should typically set to
     !> something other than zero.
     real :: aky_min = 0.0
     !> Sets the type of spacing between ky grid points, available options are :
     !>
     !> -  'default' : Same as 'linear'
     !> -  'exponential' : Evenly spaced in log(ky).
     !> -  'linear' : Evenly spaced in ky.
     !>
     character(len = 20) :: kyspacing_option = 'default'
     !> Maximum toroidal mode number. Can use instead of `aky_max`.
     integer :: n0_max = 0
     !> Minimum toroidal mode number. Can use instead of `aky_min`.
     integer :: n0_min = 0
     !> The number of 'actual' ky modes.
     integer :: naky = 1
     !> Number of toroidal modes, only used if `n0_min`>0. Overrides `naky`
     !> in kt_grids_range_parameters. Note if the number of modes isn't compatible
     !> with the requested min and max toroidal mode numbers then we just run with
     !> one mode, determined by `n0_max`.
     integer :: nn0 = 1
     !> Number of \(\theta_0\) (kx) modes
     integer :: ntheta0 = 1
     !> Used to convert `n0_min`, `n0_max` range into `aky_min`,
     !> `aky_max`, if `n0_min`>0. If `n0_min` is set,
     !> `aky_min=n0_min*drhodpsi*rhostar_range` and
     !> `aky_max=n0_max*drhodpsi*rhostar_range` where `drhodpsi` is
     !> calculated as a part of the geometry setup.
     real :: rhostar_range = 1.0e-4
     !> Upper limit of `theta_0` range
     real :: theta0_max = 0.0
     !> Lower limit of `theta_0` range
     real :: theta0_min = 0.0
   contains
     procedure, public :: read => read_kt_grids_range_config
     procedure, public :: write => write_kt_grids_range_config
     procedure, public :: reset => reset_kt_grids_range_config
     procedure, public :: broadcast => broadcast_kt_grids_range_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_range_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_range_config
  end type kt_grids_range_config_type

  type(kt_grids_range_config_type) :: kt_grids_range_config

contains

  !> FIXME : Add documentation
  subroutine read_parameters_range(kt_grids_range_config_in)
    use file_utils, only: input_unit, input_unit_exist, error_unit
    use text_options, only: text_option, get_option_value
    implicit none
    type(kt_grids_range_config_type), intent(in), optional :: kt_grids_range_config_in
    integer :: ierr
    logical :: exist
    type (text_option), dimension(3), parameter :: kyspacingopts = &
         (/ text_option('default', kyspacingopt_linear), &
         text_option('linear', kyspacingopt_linear), &
         text_option('exponential', kyspacingopt_exp) /)

    if (parameters_read) return
    parameters_read = .true.

    if (present(kt_grids_range_config_in)) kt_grids_range_config = kt_grids_range_config_in

    call kt_grids_range_config%init(name = 'kt_grids_range_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    akx_max = kt_grids_range_config%akx_max
    akx_min = kt_grids_range_config%akx_min
    aky_max = kt_grids_range_config%aky_max
    aky_min = kt_grids_range_config%aky_min
    kyspacing_option = kt_grids_range_config%kyspacing_option
    n0_max = kt_grids_range_config%n0_max
    n0_min = kt_grids_range_config%n0_min
    naky = kt_grids_range_config%naky
    nn0 = kt_grids_range_config%nn0
    ntheta0 = kt_grids_range_config%ntheta0
    ntheta0 = kt_grids_range_config%ntheta0
    rhostar_range = kt_grids_range_config%rhostar_range
    theta0_max = kt_grids_range_config%theta0_max
    theta0_min = kt_grids_range_config%theta0_min

    exist = kt_grids_range_config%exist

    ierr = error_unit()
    call get_option_value(kyspacing_option, kyspacingopts, kyspacingopt_switch,&
         ierr, "kyspacing_option in kt_grids_range_parameters",.true.)
  end subroutine read_parameters_range

  !> FIXME : Add documentation
  subroutine finish_parameters_range
    implicit none
    parameters_read = .false.
    initialized = .false.
    call kt_grids_range_config%reset()
  end subroutine finish_parameters_range

  !> FIXME : Add documentation
  subroutine init_kt_grids_range(kt_grids_range_config_in)
    !CMR, 14/10/2013:
    ! New namelist variables nn0, n0_min, n0_max, rhostar_range to set ky grid
    !                                             using toroidal mode numbers.
    ! Toroidal modenumbers are used if n0_min> 0 prescribed in input file.
    use theta_grid, only: drhodpsi
    use file_utils, only: error_unit
    implicit none
    type(kt_grids_range_config_type), intent(in), optional :: kt_grids_range_config_in
    integer :: ierr
    !! Temporary value of n0 as floating point
    real :: n0_tmp
    !! Temporary variable for swapping n0_{min,max}
    integer :: n0_swap

    if (initialized) return
    initialized = .false.

    call read_parameters_range(kt_grids_range_config_in)

    ierr = error_unit()

    !Override kyspacing_option in certain cases
    select case (kyspacingopt_switch)
    case (kyspacingopt_exp)
       if(aky_min.le.0) then
          write(ierr,'("Cannot use kyspacing_option=",A," with aky_min<=0.0 --> setting to",A)') &
               "'exponential'","'linear'"
          kyspacingopt_switch=kyspacingopt_linear
       endif
    end select

    if (n0_min .gt. 0) then
       !CMR if n0_min>0 then override aky inputs and use nn0, n0_min, n0_max to determine aky range

       !Important to only do the following check and fix if nn0 > 0 as the second argument of `mod`
       !must be non-zero. Failing to guard against this can lead to different outcomes with different
       !compilers
       if(nn0 > 1) then

          !If toroidal mode number range would lead to non-integer mode numbers then
          !we adjust the range to try to fix this.
          if(mod(n0_max-n0_min,nn0-1).ne.0) then
             !Give a warning message that we're changing things
             write(ierr,'("Warning: toroidal mode number range setup would lead to non-integer")')
             write(ierr,'("         mode numbers --> Attempting to adjust range slightly.")')

             !n0_max should be n0_min + I*(nn0-1) or n0_min + (I+1)*(nn0-1)
             !where I is int(n0_max-n0_min/(nn0-1)) and we add 1 if the
             !remainder (n0_max-n0_min)/(nn0-1) - I is > 0.5
             !Note it's nn0-1 rather than nn0 as this is number of intervals

             !First calculate the floating step size
             n0_tmp = (n0_max-n0_min*1.0)/(nn0-1)

             !Now construct the new upper limit, int(n0_tmp) = I
             !nint(n0_tmp-int(n0_tmp)) should be either 0 or 1 depending
             !on if the remainder is < or > 0.5
             n0_max = n0_min + (int(n0_tmp)+nint(n0_tmp-int(n0_tmp)))*(nn0-1)

             !Double check it's all OK now, should always be fine but just
             !putting this here in case of logic error or strange cases.
             if(mod(n0_max-n0_min,nn0-1).ne.0) then
                write(ierr,'("Error: Attempt to fix toroidal mode number range failed. Forcing nn0=1")')
                nn0=1
                n0_max = n0_min
             endif
          endif
       end if

       !If n0_min>n0_max swap values
       if(n0_min.gt.n0_max) then
          write(ierr,'("Warning: Swapping max and min n0 values")')
          n0_swap = n0_min
          n0_min = n0_max
          n0_max = n0_swap
       endif

       !If n0_min == n0_max ensure nn0=1
       if(n0_min.eq.n0_max) then
          if(nn0.ne.1) write(ierr,'("Warning: Forcing nn0=1 as n0_min==n0_max.")')
          nn0 = 1
       endif

       !If there's only one mode then force n0_max=n0_min
       if(nn0 .eq. 1) then
          n0_max = n0_min
       endif

       !Set the upper and lower aky limits
       aky_max=n0_max*drhodpsi*rhostar_range
       aky_min=n0_min*drhodpsi*rhostar_range

       !Set the number of aky values
       naky=nn0
    endif

  end subroutine init_kt_grids_range

  !> FIXME : Add documentation
  subroutine wnml_kt_grids_range(unit)
    implicit none
    integer, intent(in) :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_range_parameters"
    write (unit, fmt="(' naky = ',i3)") naky
    write (unit, fmt="(' aky_min = ',e17.10)") aky_min
    write (unit, fmt="(' aky_max = ',e17.10)") aky_max
    write (unit, fmt="(' nn0 = ',i3)") nn0
    write (unit, fmt="(' n0_min = ',i10)") n0_min
    write (unit, fmt="(' n0_max = ',i10)") n0_max
    write (unit, fmt="(' rhostar_range = ',e17.10)") rhostar_range
    write (unit, fmt="(' ntheta0 = ',i3)") ntheta0
    write (unit, fmt="(' theta0_min = ',e17.10)") theta0_min
    write (unit, fmt="(' theta0_max = ',e17.10)") theta0_max
    write (unit, fmt="(' akx_min = ',e17.10)") akx_min
    write (unit, fmt="(' akx_max = ',e17.10)") akx_max
    select case(kyspacingopt_switch)
    case (kyspacingopt_linear)
       write (unit, fmt="(' kyspacing_option = ',A)") "linear"
    case (kyspacingopt_exp)
       write (unit, fmt="(' kyspacing_option = ',A)") "exponential"
    end select
    write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_range

  !> FIXME : Add documentation
  subroutine range_get_sizes (naky_x, ntheta0_x, nx, ny)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx, ny
    naky_x = naky  ;  ntheta0_x = ntheta0
    nx = 0         ;  ny = 0
  end subroutine range_get_sizes

  !> FIXME : Add documentation
  subroutine range_get_grids (aky, theta0, akx, ikx)
    ! BD: Could add some logic here to set theta0 if akx is given?  When do we need what?
    use theta_grid, only: shat
    use mp, only: mp_abort
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    !> Discrete kx wavenumber grid indices
    integer, dimension (:), intent (out) :: ikx

    real :: dkx, dky, dtheta0
    integer :: i, j

    if ( size(aky) /= naky) then
       call mp_abort('range_get_grids: size(aky) /= naky',.true.)
    endif

    if ( size(akx) /= ntheta0) then
       call mp_abort('range_get_grids: size(akx) /= ntheta0',.true.)
    endif

    dky = 0.0
    if (naky > 1)then
       select case (kyspacingopt_switch)
       case (kyspacingopt_linear)
          dky = (aky_max - aky_min)/real(naky - 1)
          aky = (/ (aky_min + dky*real(i), i = 0,naky-1) /)
       case (kyspacingopt_exp)
          dky = (log(aky_max) - log(aky_min))/real(naky - 1)
          aky = (/ (exp(log(aky_min) + dky*real(i)), i = 0,naky-1) /)
       end select
    else
       aky = (/ (aky_min, i = 0,naky-1) /)
    endif

    ! set default theta0 to 0
    theta0=0.0

    !
    ! BD: Assumption here differs from convention that abs(shat) <= 1.e-5 triggers periodic bc
    !
    if (shat /= 0.0d0) then  ! ie assumes boundary_option .eq. 'linked'
       dtheta0 = 0.0
       if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)

       do j = 1, naky
          theta0(:,j) &
               = (/ (theta0_min + dtheta0*real(i), i=0,ntheta0-1) /)
       end do

       !<DD>Adding support for ky=0, kx/=0
       if(aky(1)==0)then
          if(naky>1)then
             akx = theta0(:,2) * shat * aky(2)
          else
             dkx = 0.0
             if (ntheta0 > 1) dkx = (akx_max - akx_min)/real(ntheta0 - 1)
             akx = (/ (akx_min + dkx*real(i), i = 0,ntheta0-1) /)
          end if
       else
          !This is the original behaviour
          akx = theta0(:,1) * shat * aky(1)
       endif
    else

       !CMR, 22/9/2010:  ie here assume boundary_option .eq. 'periodic'
       !new code for periodic finite kx ballooning space runs with shat=0
       dkx = 0.0
       if (ntheta0 > 1) dkx = (akx_max - akx_min)/real(ntheta0 - 1)
       akx = (/ (akx_min + dkx*real(i), i = 0,ntheta0-1) /)
    endif

    do j = 1, ntheta0
       ikx(j) = j - 1
    end do
  end subroutine range_get_grids

  !> FIXME : Add documentation
  subroutine check_kt_grids_range(report_unit)
    use constants, only: twopi
    use theta_grid, only: shat
    implicit none
    integer, intent(in) :: report_unit
    real :: dtheta0
    integer :: i, j
    real, dimension(:), allocatable:: aky, akx
    real, dimension(:,:), allocatable:: theta0
    integer, dimension(:), allocatable :: ikx

    write (report_unit, *)
    write (report_unit, fmt="('A range of k_perps will be evolved.')")
    if (n0_min .gt.0) write (report_unit, fmt="('ky set using toroidal mode numbers with n0_min=',i8/T34,'rhostar_range=',1pe12.4)") n0_min,rhostar_range
    write (report_unit, *)
    write (report_unit, fmt="('There are ',i3,' values of ky rho and ',i3,' values of theta_0/kx rho:')") naky, ntheta0
    write (report_unit, *)

    !<DD>Calculate the kt grids
    allocate(aky(naky),theta0(ntheta0,naky),akx(ntheta0), ikx(ntheta0))
    call range_get_grids(aky, theta0, akx, ikx)

    !Report grid values
    do j = 1, naky
       do i = 1, ntheta0
          write (report_unit, fmt="('ky rho = ',e11.4,' theta0 = ',e11.4,' kx rho = ',e11.4)") &
               aky(j),theta0(i,j),akx(i)
       end do
    end do
    deallocate(aky,theta0,akx)

    ! CMR, add some !!!error checking!!! for ballooning space runs for shat /= 0
    ! using flow shear: check that the constraints on theta0 grid are satisfied!

    if (shat /= 0) then
       !It would be nice to only write this information if g_exb*gexbfac/=0 but currently
       !dependencies prevent this.
       dtheta0 = 0.0    ;  if (ntheta0 > 1) dtheta0 = (theta0_max - theta0_min)/real(ntheta0 - 1)
       if (abs(mod(twopi-theta0_max+theta0_min,twopi)-dtheta0) > 1.0e-3*dtheta0) then
          write (report_unit, *)
          write (report_unit, fmt="('IF using perp ExB flow shear in BALLOONING SPACE there is an ERROR that will corrupt results.')")
          write (report_unit, fmt="('check_kt_grids_range: inappropriate theta0 grid')")
          write (report_unit, fmt="('In ballooning space with sheared flow, 2pi-theta0_max+theta0_min =',e11.4,' must be set equal to dtheta = ',e11.4)") twopi-theta0_max+theta0_min, dtheta0
       endif
    endif

  end subroutine check_kt_grids_range

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_range_config(kt_grids_range_config_in)
    use mp, only: mp_abort
    type(kt_grids_range_config_type), intent(in), optional :: kt_grids_range_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_range_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_range_config_in)) then
       kt_grids_range_config = kt_grids_range_config_in
    end if
  end subroutine set_kt_grids_range_config

#include "kt_grids_range_auto_gen.inc"
end module kt_grids_range
