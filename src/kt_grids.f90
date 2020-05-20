!> Set up values of kx and ky for linear runs that use a single k_perp mode.
module kt_grids_single
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: init_kt_grids_single, single_get_sizes, single_get_grids
  public :: check_kt_grids_single, wnml_kt_grids_single 
  public :: read_parameters_single
  public :: finish_parameters_single

  public :: kt_grids_single_config_type
  public :: set_kt_grids_single_config
  
  real :: akx, aky, theta0, rhostar_single
  integer :: n0
  logical :: parameters_read = .false.
  logical :: initialized = .false.
  
  !> Used to represent the input configuration of kt_grids_single
  type, extends(abstract_config_type) :: kt_grids_single_config_type
     ! namelist : kt_grids_single_parameters
     ! indexed : false     
     !> \(k_x \rho\) for the reference species (but recommended to set `theta0` instead).
     real :: akx = 0.0
     !> \(k_y \rho\) for the reference species.
     real :: aky = 0.4
     !> if `n0`>0 use toroidal mode number to override `aky` and set
     !> `aky=n0*drhodpsi*rhostar_single` where `drhodpsi` is calculated
     !> as a part of the geometry setup.
     integer :: n0 = 0
     !> Used in conjunction with `n0`:
     !> `aky=n0*drhodpsi*rhostar_single` (if `n0` is set) where
     !> `drhodpsi` is calculated as a part of the geometry setup.
     real :: rhostar_single = 1.0e-4
     !> \(\theta_0\) is the ballooning angle, sets the point in
     !> \(\theta\) where the radial wavenumber is zero
     real :: theta0 = 0.0
   contains
     procedure, public :: read => read_kt_grids_single_config
     procedure, public :: write => write_kt_grids_single_config
     procedure, public :: reset => reset_kt_grids_single_config
     procedure, public :: broadcast => broadcast_kt_grids_single_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_single_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_single_config
  end type kt_grids_single_config_type
  
  type(kt_grids_single_config_type) :: kt_grids_single_config
  
contains

  !> FIXME : Add documentation
  subroutine read_parameters_single(kt_grids_single_config_in)
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in
    logical :: exist

    if (parameters_read) return
    parameters_read = .true.

    if (present(kt_grids_single_config_in)) kt_grids_single_config = kt_grids_single_config_in

    call kt_grids_single_config%init(name = 'kt_grids_single_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    akx = kt_grids_single_config%akx
    aky = kt_grids_single_config%aky
    n0 = kt_grids_single_config%n0
    rhostar_single = kt_grids_single_config%rhostar_single
    theta0 = kt_grids_single_config%theta0

    exist = kt_grids_single_config%exist
  end subroutine read_parameters_single

  !> FIXME : Add documentation
  subroutine finish_parameters_single
    implicit none
    parameters_read = .false.
    initialized = .false.
    call kt_grids_single_config%reset()
  end subroutine finish_parameters_single

  !> FIXME : Add documentation  
  subroutine init_kt_grids_single(kt_grids_single_config_in)
!CMR, 14/10/2013: 
! New namelist variables n0, rhostar_single to set aky using toroidal mode number.
! Toroidal modenumber used if n0> 0 prescribed in input file. 
    use theta_grid, only: drhodpsi
    implicit none
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in

    if(initialized) return
    initialized = .true.
    
    call read_parameters_single(kt_grids_single_config_in)

    if (n0 .gt. 0) then
!CMR if n0>0 then override aky inputs and use n0 to determine aky
       aky=n0*drhodpsi*rhostar_single
    endif

  end subroutine init_kt_grids_single

  !> FIXME : Add documentation  
  subroutine wnml_kt_grids_single(unit)
    implicit none
    integer, intent(in) :: unit
    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_single_parameters"
    write (unit, fmt="(' aky = ',e17.10)") aky
    write (unit, fmt="(' theta0 = ',e17.10)") theta0
    write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_single

  !> FIXME : Add documentation  
  subroutine single_get_sizes (naky, ntheta0, nx, ny)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny

    naky = 1  ;  ntheta0 = 1
    nx = 0    ;  ny = 0

  end subroutine single_get_sizes

  !> FIXME : Add documentation  
  subroutine single_get_grids (aky_out, theta0_out, akx_out)
    implicit none
    real, dimension (:), intent (out) :: aky_out, akx_out
    real, dimension (:,:), intent (out) :: theta0_out

    aky_out = aky
    theta0_out = theta0
    akx_out = akx
  end subroutine single_get_grids

  !> FIXME : Add documentation  
  subroutine check_kt_grids_single(report_unit)
    implicit none
    integer, intent(in) :: report_unit

    write (report_unit, *) 
    write (report_unit, fmt="('A single k_perp will be evolved, with: ')")
    if (n0 .gt.0) write (report_unit, fmt="('ky set using toroidal mode number, n0=',i8/T24,'rhostar_single=',1pe12.4)") n0, rhostar_single
    write (report_unit, *) 
    write (report_unit, fmt="('ky rho = ',f10.4)") aky
    write (report_unit, fmt="('theta_0 = ',f10.4)") theta0
    if (akx /= 0.) then
       write (report_unit, *) 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, fmt="('The value of akx in the kt_grids_single_parameters namelist is ignored.')") 
       write (report_unit, fmt="('You have set akx to a non-zero value.')") 
       write (report_unit, fmt="('THIS IS PROBABLY AN ERROR.')") 
       write (report_unit, fmt="('################# WARNING #######################')")
       write (report_unit, *) 
    end if
  end subroutine check_kt_grids_single

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_single_config(kt_grids_single_config_in)
    use mp, only: mp_abort
    type(kt_grids_single_config_type), intent(in), optional :: kt_grids_single_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_single_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_single_config_in)) then
       kt_grids_single_config = kt_grids_single_config_in
    end if
  end subroutine set_kt_grids_single_config
  
  !---------------------------------------
  ! Following is for the config_type
  !--------------------------------------- 
  
  !> Reads in the kt_grids_single_parameters namelist and populates the member variables
  subroutine read_kt_grids_single_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_single_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: n0
    real :: akx, aky, rhostar_single, theta0

    namelist /kt_grids_single_parameters/ akx, aky, n0, rhostar_single, theta0

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    akx = self%akx
    aky = self%aky
    n0 = self%n0
    rhostar_single = self%rhostar_single
    theta0 = self%theta0

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_single_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_single_parameters)
    endif

    ! Now copy from local variables into type members
    self%akx = akx
    self%aky = aky
    self%n0 = n0
    self%rhostar_single = rhostar_single
    self%theta0 = theta0

    self%exist = exist
  end subroutine read_kt_grids_single_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_single_config(self, unit)
    implicit none
    class(kt_grids_single_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("akx", self%akx, unit_internal)
    call self%write_key_val("aky", self%aky, unit_internal)
    call self%write_key_val("n0", self%n0, unit_internal)
    call self%write_key_val("rhostar_single", self%rhostar_single, unit_internal)
    call self%write_key_val("theta0", self%theta0, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_single_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_single_config(self)
    class(kt_grids_single_config_type), intent(in out) :: self
    type(kt_grids_single_config_type) :: empty
    select type (self)
    type is (kt_grids_single_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_single_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_single_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_single_config_type), intent(in out) :: self
    call broadcast(self%akx)
    call broadcast(self%aky)
    call broadcast(self%n0)
    call broadcast(self%rhostar_single)
    call broadcast(self%theta0)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_single_config

  !> Gets the default name for this namelist
  function get_default_name_kt_grids_single_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_single_config
    get_default_name_kt_grids_single_config = "kt_grids_single_parameters"
  end function get_default_name_kt_grids_single_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_single_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_single_config
    get_default_requires_index_kt_grids_single_config = .false.
  end function get_default_requires_index_kt_grids_single_config
end module kt_grids_single

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
  subroutine range_get_grids (aky, theta0, akx)
    ! BD: Could add some logic here to set theta0 if akx is given?  When do we need what?
    use theta_grid, only: shat
    use mp, only: mp_abort
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0

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

    write (report_unit, *) 
    write (report_unit, fmt="('A range of k_perps will be evolved.')")
    if (n0_min .gt.0) write (report_unit, fmt="('ky set using toroidal mode numbers with n0_min=',i8/T34,'rhostar_range=',1pe12.4)") n0_min,rhostar_range
    write (report_unit, *) 
    write (report_unit, fmt="('There are ',i3,' values of ky rho and ',i3,' values of theta_0/kx rho:')") naky, ntheta0
    write (report_unit, *) 

    !<DD>Calculate the kt grids
    allocate(aky(naky),theta0(ntheta0,naky),akx(ntheta0))
    call range_get_grids(aky, theta0, akx)

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
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the kt_grids_range_parameters namelist and populates the member variables
  subroutine read_kt_grids_range_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_range_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    character(len = 20) :: kyspacing_option
    integer :: n0_max, n0_min, naky, nn0, ntheta0
    real :: akx_max, akx_min, aky_max, aky_min, rhostar_range, theta0_max, theta0_min

    namelist /kt_grids_range_parameters/ akx_max, akx_min, aky_max, aky_min, kyspacing_option, n0_max, n0_min, naky, nn0, ntheta0, rhostar_range, theta0_max, theta0_min

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    akx_max = self%akx_max
    akx_min = self%akx_min
    aky_max = self%aky_max
    aky_min = self%aky_min
    kyspacing_option = self%kyspacing_option
    n0_max = self%n0_max
    n0_min = self%n0_min
    naky = self%naky
    nn0 = self%nn0
    ntheta0 = self%ntheta0
    rhostar_range = self%rhostar_range
    theta0_max = self%theta0_max
    theta0_min = self%theta0_min

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_range_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_range_parameters)
    endif

    ! Now copy from local variables into type members
    self%akx_max = akx_max
    self%akx_min = akx_min
    self%aky_max = aky_max
    self%aky_min = aky_min
    self%kyspacing_option = kyspacing_option
    self%n0_max = n0_max
    self%n0_min = n0_min
    self%naky = naky
    self%nn0 = nn0
    self%ntheta0 = ntheta0
    self%rhostar_range = rhostar_range
    self%theta0_max = theta0_max
    self%theta0_min = theta0_min

    self%exist = exist
  end subroutine read_kt_grids_range_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_range_config(self, unit)
    implicit none
    class(kt_grids_range_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("akx_max", self%akx_max, unit_internal)
    call self%write_key_val("akx_min", self%akx_min, unit_internal)
    call self%write_key_val("aky_max", self%aky_max, unit_internal)
    call self%write_key_val("aky_min", self%aky_min, unit_internal)
    call self%write_key_val("kyspacing_option", self%kyspacing_option, unit_internal)
    call self%write_key_val("n0_max", self%n0_max, unit_internal)
    call self%write_key_val("n0_min", self%n0_min, unit_internal)
    call self%write_key_val("naky", self%naky, unit_internal)
    call self%write_key_val("nn0", self%nn0, unit_internal)
    call self%write_key_val("ntheta0", self%ntheta0, unit_internal)
    call self%write_key_val("rhostar_range", self%rhostar_range, unit_internal)
    call self%write_key_val("theta0_max", self%theta0_max, unit_internal)
    call self%write_key_val("theta0_min", self%theta0_min, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_range_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_range_config(self)
    class(kt_grids_range_config_type), intent(in out) :: self
    type(kt_grids_range_config_type) :: empty
    select type (self)
    type is (kt_grids_range_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_range_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_range_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_range_config_type), intent(in out) :: self
    call broadcast(self%akx_max)
    call broadcast(self%akx_min)
    call broadcast(self%aky_max)
    call broadcast(self%aky_min)
    call broadcast(self%kyspacing_option)
    call broadcast(self%n0_max)
    call broadcast(self%n0_min)
    call broadcast(self%naky)
    call broadcast(self%nn0)
    call broadcast(self%ntheta0)
    call broadcast(self%rhostar_range)
    call broadcast(self%theta0_max)
    call broadcast(self%theta0_min)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_range_config

  !> Gets the default name for this namelist
  function get_default_name_kt_grids_range_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_range_config
    get_default_name_kt_grids_range_config = "kt_grids_range_parameters"
  end function get_default_name_kt_grids_range_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_range_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_range_config
    get_default_requires_index_kt_grids_range_config = .false.
  end function get_default_requires_index_kt_grids_range_config
end module kt_grids_range

!> Set up sets of (kx, ky) values for linear runs.
module kt_grids_specified
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: init_kt_grids_specified, specified_get_sizes, specified_get_grids
  public :: check_kt_grids_specified, wnml_kt_grids_specified, finish_parameters_specified

  public :: kt_grids_specified_config_type 
  public :: kt_grids_specified_element_config_type
  public :: set_kt_grids_specified_config
  
  integer :: naky, ntheta0, nx, ny
  logical :: initialized = .false.
  
  !> Used to represent the input configuration of kt_grids_specified
  type, extends(abstract_config_type) :: kt_grids_specified_config_type
     ! namelist : kt_grids_specified_parameters
     ! indexed : false     
     !> The number of `ky` values evolved. The total number of modes evolved is given by
     !> the maximum of `naky` and `ntheta0`. One must also set up the appropriate number of
     !> [[kt_grids_specified_element]]`_<i>` namelists.
     integer :: naky = 1
     !> The number of `theta0` values evolved. The total number of modes evolved is given by
     !> the maximum of `naky` and `ntheta0`. One must also set up the appropriate number of
     !> [[kt_grids_specified_element]]`_<i>` namelists.
     integer :: ntheta0 = 1
     !> Mostly ignored. Does not set the number of modes used.
     !>
     !> @todo Consider removing
     integer :: nx = 0
     !> Mostly ignored. Does not set the number of modes used.
     !>
     !> @todo Consider removing
     integer :: ny = 0
   contains
     procedure, public :: read => read_kt_grids_specified_config
     procedure, public :: write => write_kt_grids_specified_config
     procedure, public :: reset => reset_kt_grids_specified_config
     procedure, public :: broadcast => broadcast_kt_grids_specified_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_specified_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_specified_config
  end type kt_grids_specified_config_type
  
  type(kt_grids_specified_config_type) :: kt_grids_specified_config

  !> Used to represent the input configuration of a specific specified
  !> wavenumber pair
  type, extends(abstract_config_type) :: kt_grids_specified_element_config_type
     ! namelist : kt_grids_specified_element
     ! indexed : true     
     !> Sets the \(k_y \rho\) value for this wavenumber element.
     real :: aky = 0.4
     !> Sets the \(k_x \rho\) value for this wavenumber element (but should set `theta0` instead).
     real :: akx = 0.0
     !> Sets the \(\theta_0\) value for this wavenumber element.
     real :: theta0 = 0.0
   contains
     procedure, public :: read => read_kt_grids_specified_element_config
     procedure, public :: write => write_kt_grids_specified_element_config
     procedure, public :: reset => reset_kt_grids_specified_element_config
     procedure, public :: broadcast => broadcast_kt_grids_specified_element_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_specified_element_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_specified_element_config
  end type kt_grids_specified_element_config_type
  
  type(kt_grids_specified_element_config_type), dimension(:), allocatable :: kt_grids_specified_element_config
  
contains

  !> FIXME : Add documentation
  subroutine init_kt_grids_specified(kt_grids_specified_config_in)
    use file_utils, only: input_unit, input_unit_exist
    implicit none
    type(kt_grids_specified_config_type), intent(in), optional :: kt_grids_specified_config_in
    logical :: exist

    if (initialized) return
    initialized = .true.
    
    if (present(kt_grids_specified_config_in)) kt_grids_specified_config = kt_grids_specified_config_in

    call kt_grids_specified_config%init(name = 'kt_grids_specified_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    naky = kt_grids_specified_config%naky
    ntheta0 = kt_grids_specified_config%ntheta0
    nx = kt_grids_specified_config%nx
    ny = kt_grids_specified_config%ny

    exist = kt_grids_specified_config%exist

  end subroutine init_kt_grids_specified

  !> FIXME : Add documentation  
  subroutine finish_parameters_specified
    implicit none
    initialized = .false.
    call kt_grids_specified_config%reset()
    if(allocated(kt_grids_specified_element_config)) deallocate(kt_grids_specified_element_config)
  end subroutine finish_parameters_specified
  
  !> FIXME : Add documentation  
  subroutine wnml_kt_grids_specified (unit, aky, theta0)
    implicit none
    integer, intent(in) :: unit
    real, dimension(:), intent(in) :: aky
    real, dimension(:,:), intent(in) :: theta0
    integer :: i
    character (200) :: line

    call kt_grids_specified_config%write(unit)
    do i=1,max(naky,ntheta0)
       write (unit, *)
       write (line, *) i
       write (unit, fmt="(' &',a)") &
            & trim("kt_grids_specified_element_"//trim(adjustl(line)))
       write (unit, fmt="(' aky = ',e13.6,' theta0 = ',e13.6,'  /')") aky(i), theta0(i,1)
       write (unit, fmt="(' /')")
    end do
  end subroutine wnml_kt_grids_specified

  !> FIXME : Add documentation  
  subroutine specified_get_sizes (naky_x, ntheta0_x, nx_x, ny_x)
    implicit none
    integer, intent (out) :: naky_x, ntheta0_x, nx_x, ny_x

    naky_x = naky  ;   ntheta0_x = ntheta0
    nx_x = nx      ;   ny_x = ny

 end subroutine specified_get_sizes

 !> FIXME : Add documentation
 !>
 !> @todo : Consider moving the content of this routine to init_kt_grids_specified
 subroutine specified_get_grids (aky, theta0, akx, kt_grids_specified_element_config_in)
   use mp, only: proc0
   implicit none
   real, dimension (:), intent (out) :: akx, aky
   real, dimension (:,:), intent (out) :: theta0
   type(kt_grids_specified_element_config_type), intent(in), dimension(:), optional :: kt_grids_specified_element_config_in    
   real :: aky_dummy, theta0_dummy, akx_dummy
   integer :: i, naky, ntheta0, nmodes

   ! Note we shadow the module naky/ntheta variables here
   naky = size(aky)  ;  ntheta0 = size(akx)
   nmodes = max(naky,ntheta0)
   
   if(.not.allocated(kt_grids_specified_element_config)) allocate(kt_grids_specified_element_config(nmodes))
   
   if (size(kt_grids_specified_element_config) .ne. nmodes) then
      if(proc0) print*,"Warning: inconsistent number of config elements."
   end if
   
   do i = 1, nmodes
      call read_element (i, aky_dummy, theta0_dummy, akx_dummy, kt_grids_specified_element_config_in)
      if (i <= naky) aky(i) = aky_dummy
      if (i <= ntheta0) theta0(i,:) = theta0_dummy
      if (i <= ntheta0) akx(i) = akx_dummy
    end do
    
  end subroutine specified_get_grids
  
  !> FIXME : Add documentation  
  subroutine read_element (i, aky_dummy, theta0_dummy, akx_dummy, kt_grids_specified_element_config_in)
    use file_utils, only: get_indexed_namelist_unit
    implicit none
    integer, intent (in) :: i
    real, intent (out) :: aky_dummy, theta0_dummy, akx_dummy
    type(kt_grids_specified_element_config_type), intent(in), dimension(:), optional :: kt_grids_specified_element_config_in    
    real :: akx, aky, theta0
    logical :: exist

    if (present(kt_grids_specified_element_config_in)) kt_grids_specified_element_config(i) = kt_grids_specified_element_config_in(i)

    call kt_grids_specified_element_config(i)%init(name = 'kt_grids_specified_element', &
         requires_index = .true., index = i)

    ! Copy out internal values into module level parameters
    akx = kt_grids_specified_element_config(i)%akx
    aky = kt_grids_specified_element_config(i)%aky
    theta0 = kt_grids_specified_element_config(i)%theta0
    
    exist = kt_grids_specified_element_config(i)%exist
    
    aky_dummy = aky
    theta0_dummy = theta0
    akx_dummy = akx
  end subroutine read_element

  subroutine check_kt_grids_specified (report_unit, aky, theta0)
    implicit none
    integer, intent(in) :: report_unit
    real, dimension (:), intent (in) :: aky
    real, dimension (:), intent (in) :: theta0
    integer :: i

    write (report_unit, *) 
    write (report_unit, fmt="('A set of ',i3,' k_perps will be evolved.')") max(naky,ntheta0)
    write (report_unit, *) 
    do i=1, max(naky,ntheta0)
       write (report_unit, fmt="('ky rho = ',e11.4,' theta0 = ',e11.4)") aky(i), theta0(i)
    end do
  end subroutine check_kt_grids_specified

  !> Set the module level config types
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_specified_config(kt_grids_specified_config_in, kt_grids_specified_element_config_in)
    use mp, only: mp_abort
    type(kt_grids_specified_config_type), intent(in), optional :: kt_grids_specified_config_in
    type(kt_grids_specified_element_config_type), intent(in), dimension(:), optional :: kt_grids_specified_element_config_in    

    if (initialized) then
       call mp_abort("Trying to set kt_grids_specified_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_specified_config_in)) then
       kt_grids_specified_config = kt_grids_specified_config_in
    end if
    if (present(kt_grids_specified_element_config_in)) then
       kt_grids_specified_element_config = kt_grids_specified_element_config_in
    end if    
  end subroutine set_kt_grids_specified_config   
  
  !---------------------------------------
  ! Following is for the config_type
  !--------------------------------------- 
  
  !> Reads in the kt_grids_specified_parameters namelist and populates the member variables
  subroutine read_kt_grids_specified_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_specified_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    integer :: naky, ntheta0, nx, ny

    namelist /kt_grids_specified_parameters/ naky, ntheta0, nx, ny

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    naky = self%naky
    ntheta0 = self%ntheta0
    nx = self%nx
    ny = self%ny

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_specified_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_specified_parameters)
    endif

    ! Now copy from local variables into type members
    self%naky = naky
    self%ntheta0 = ntheta0
    self%nx = nx
    self%ny = ny

    self%exist = exist
  end subroutine read_kt_grids_specified_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_specified_config(self, unit)
    implicit none
    class(kt_grids_specified_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("naky", self%naky, unit_internal)
    call self%write_key_val("ntheta0", self%ntheta0, unit_internal)
    call self%write_key_val("nx", self%nx, unit_internal)
    call self%write_key_val("ny", self%ny, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_specified_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_specified_config(self)
    class(kt_grids_specified_config_type), intent(in out) :: self
    type(kt_grids_specified_config_type) :: empty
    select type (self)
    type is (kt_grids_specified_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_specified_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_specified_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_specified_config_type), intent(in out) :: self
    call broadcast(self%naky)
    call broadcast(self%ntheta0)
    call broadcast(self%nx)
    call broadcast(self%ny)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_specified_config
  
  !> Gets the default name for this namelist
  function get_default_name_kt_grids_specified_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_specified_config
    get_default_name_kt_grids_specified_config = "kt_grids_specified_parameters"
  end function get_default_name_kt_grids_specified_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_specified_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_specified_config
    get_default_requires_index_kt_grids_specified_config = .false.
  end function get_default_requires_index_kt_grids_specified_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the kt_grids_specified_element namelist and populates the member variables
  subroutine read_kt_grids_specified_element_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_specified_element_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    real :: akx, aky, theta0

    namelist /kt_grids_specified_element/ akx, aky, theta0

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    akx = self%akx
    aky = self%aky
    theta0 = self%theta0

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_specified_element, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_specified_element)
    endif

    ! Now copy from local variables into type members
    self%akx = akx
    self%aky = aky
    self%theta0 = theta0

    self%exist = exist
  end subroutine read_kt_grids_specified_element_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_specified_element_config(self, unit)
    implicit none
    class(kt_grids_specified_element_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("akx", self%akx, unit_internal)
    call self%write_key_val("aky", self%aky, unit_internal)
    call self%write_key_val("theta0", self%theta0, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_specified_element_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_specified_element_config(self)
    class(kt_grids_specified_element_config_type), intent(in out) :: self
    type(kt_grids_specified_element_config_type) :: empty
    select type (self)
    type is (kt_grids_specified_element_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_specified_element_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_specified_element_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_specified_element_config_type), intent(in out) :: self
    call broadcast(self%akx)
    call broadcast(self%aky)
    call broadcast(self%theta0)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_specified_element_config

  !> Gets the default name for this namelist
  function get_default_name_kt_grids_specified_element_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_specified_element_config
    get_default_name_kt_grids_specified_element_config = "kt_grids_specified_element"
  end function get_default_name_kt_grids_specified_element_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_specified_element_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_specified_element_config
    get_default_requires_index_kt_grids_specified_element_config = .true.
  end function get_default_requires_index_kt_grids_specified_element_config
end module kt_grids_specified

!> Set the perpendicular box size and resolution for linear or nonlinear runs.
module kt_grids_box
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  
  implicit none

  private

  public :: init_kt_grids_box, box_get_sizes, box_get_grids
  public :: read_parameters_box
  public :: finish_parameters_box
  public :: box_set_overrides
  public :: check_kt_grids_box, wnml_kt_grids_box
  public :: x0, y0, jtwist !RN> Caution: these are not broadcasted!
  public :: gryfx

  public :: kt_grids_box_config_type
  public :: set_kt_grids_box_config
  
  real :: ly, y0, x0, rtwist, rhostar_box
  integer :: naky_private, ntheta0_private, nx_private, ny_private
  integer :: nkpolar_private, n0
  integer :: jtwist
  integer :: naky, ntheta0, nx, ny, nkpolar
  logical :: gryfx
  logical :: parameters_read = .false.
  logical :: initialized = .false.
  
  !> Used to represent the input configuration of kt_grids_box
  type, extends(abstract_config_type) :: kt_grids_box_config_type
     ! namelist : kt_grids_box_parameters
     ! indexed : false
     !> FIXME: Add documentation
     !>
     !> @note We should probably remove/move this flag to some gryfx
     !> specific location.
     logical :: gryfx = .false.
     !> For finite magnetic shear determines the box size in the x
     !> direction according to \(L_x = L_y
     !> \textrm{jtwist}/2\pi\hat{s}\).  This also affects the number
     !> of connections at each `ky` when linked boundary conditions
     !> are selected in the [[dist_fn_knobs]] namelist.
     !>
     !> This gets a smart default Initialised to `jtwist =
     !> max(int(2*pi*shat+0.5,1))` so that \(L_x\approx L_y\). If the
     !> magnetic shear is less than around 0.16 then `jtwist` will
     !> default to the minimum allowed value of 1.
     integer :: jtwist = 1
     !> Sets the box size in the y direction. If set to 0 (the default) then
     !> we set `ly=2*pi*y0`.
     real :: ly = 0.0
     !> If set greater than zero (the default) then this sets the
     !> toroidal mode number of the first non-zero `ky` by overriding
     !> the value given for `y0` through `y0 =
     !> 1.0/(n0*rhostar_box*drhodpsi)` where `drhodpsi` is determined
     !> during geometry setup.
     integer :: n0 = 0
     !> The actual number of ky modes. For nonlinear runs it is
     !> recommended to use `ny` instead. If set to 0 (the default)
     !> this will be automatically set to `1 + (ny-1)/3`.
     integer :: naky = 0
     !> FIXME: Add documentation
     !>
     !> Not used anywhere so should be removed.
     !> If set to 0 (the default) then will be set to `int(real(naky-1.)*sqrt(2.))`
     !>
     !> @todo Consider for removal.
     integer :: nkpolar = 0
     !> The actual number of theta0 modes. For nonlinear runs it is
     !> recommended to use `nx` instead. If set to 0 (the default)
     !> this will be automatically set to `1 + 2*(nx-1)/3`.
     integer :: ntheta0 = 0
     !> The number of kx points in inputs to the fft routines, and
     !> hence the number of radial points used in real space
     !> calculations. This differs from the actual number of kx points
     !> simulated in the rest of the code due to the need to prevent
     !> aliasing. The number of kx modes actually simulated
     !> (`ntheta0`) is, by default, equal to `1 + 2*(nx - 1)/3`.
     integer :: nx = 0
     !> The number of ky points in inputs to the fft routines, and
     !> hence the number of binormal points used in real space
     !> calculations. This differs from the actual number of ky points
     !> simulated in the rest of the code due to the need to prevent
     !> aliasing. The number of ky modes actually simulated
     !> (`naky`) is, by default, equal to `1 + (ny - 1)/3`.
     integer :: ny = 0
     !> The rhostar (`rhoref/Lref`) to use. Only used if `n0` also set
     !> greater than zero.  If `rhostar_box` and `n0` are greater than
     !> zero then `y0=1.0/(n0*rhostar_box*drhodpsi)`, which
     !> effectively sets the minimum non-zero `ky` used in the
     !> simulation.
     real :: rhostar_box = 0.0
     !> Expert usage only -- more documentation required.
     !>
     !> Used to control the kx spacing in simulations with effectively
     !> zero shear (\(\< 10^{-5}\)) where linked boundaries are not
     !> appropriate so periodic boundaries are used. Also only used if
     !> `x0` has been set to zero (the default). If `rtwist` is set to
     !> 0.0 (the default) then it is set to the value of
     !> `jtwist`. Effectively ends up setting the box size in the x
     !> direction, as \(L_x = L_y \textrm{rtwist}\) if `rtwist > 0 `
     !> and \(L_x = L_y / \textrm{rtwist}\) if `rtwist < 0`. See
     !> [[x0]] for an alternative.
     real :: rtwist = 0.0
     !> Controls the box length in the x direction (measured in
     !> the reference Larmour radius) if the magnetic shear is small
     !> (\(\< 10^{-5}\)). The box size in the x direction is
     !> given by \(L_x = 2\pi x_0\). See [[rtwist]] for an alternative.
     real :: x0 = 0.
     !> Controls the box length in the y direction (measured in the
     !> reference Larmour radius). The box size in the y direction is
     !> given by \(L_y = 2\pi y_0\). Note if `y0` is set negative
     !> then, it is replaced with `-1/y0` and Effectively sets the
     !> minimum wavelength captured by the box.
     real :: y0 = 2.0
   contains
     procedure, public :: read => read_kt_grids_box_config
     procedure, public :: write => write_kt_grids_box_config
     procedure, public :: reset => reset_kt_grids_box_config
     procedure, public :: broadcast => broadcast_kt_grids_box_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_box_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_box_config
  end type kt_grids_box_config_type
  
  type(kt_grids_box_config_type) :: kt_grids_box_config
  
contains

  !> FIXME : Add documentation
  subroutine read_parameters_box(kt_grids_box_config_in)
    use file_utils, only: input_unit, input_unit_exist
    use theta_grid, only: init_theta_grid, shat
    use constants, only: pi
    implicit none
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in    
    integer :: in_file
    logical :: exist

    if (parameters_read) return
    parameters_read = .true. 

    call init_theta_grid

    if (present(kt_grids_box_config_in)) kt_grids_box_config = kt_grids_box_config_in

    ! Set the smart defaults
    if (.not.kt_grids_box_config%is_initialised()) then
       kt_grids_box_config%jtwist = max(int(2.0*pi*shat + 0.5),1)  ! default jtwist -- MAB
    end if

    call kt_grids_box_config%init(name = 'kt_grids_box_parameters', requires_index = .false.)

    ! Copy out internal values into module level parameters
    gryfx = kt_grids_box_config%gryfx
    jtwist = kt_grids_box_config%jtwist
    ly = kt_grids_box_config%ly
    n0 = kt_grids_box_config%n0
    naky = kt_grids_box_config%naky
    nkpolar = kt_grids_box_config%nkpolar
    ntheta0 = kt_grids_box_config%ntheta0
    nx = kt_grids_box_config%nx
    ny = kt_grids_box_config%ny
    rhostar_box = kt_grids_box_config%rhostar_box
    rtwist = kt_grids_box_config%rtwist
    x0 = kt_grids_box_config%x0
    y0 = kt_grids_box_config%y0
    
    exist = kt_grids_box_config%exist
  end subroutine read_parameters_box

  !> FIXME : Add documentation  
  subroutine finish_parameters_box
    implicit none
    parameters_read = .false.
    initialized = .false.
    call kt_grids_box_config%reset()
  end subroutine finish_parameters_box

  !> FIXME : Add documentation  
  subroutine init_kt_grids_box(kt_grids_box_config_in)
!CMR, 14/10/2013: 
! New namelist variables: n0, rhostar_box. 
! If n0 and rhostar_box defined, set ky(1) using toroidal mode number.

    use file_utils, only: error_unit
    use theta_grid, only: shat, drhodpsi
    use constants, only: pi
    use mp, only: mp_abort, proc0
    implicit none
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in

    if (initialized) return
    initialized = .true.
    
    call read_parameters_box(kt_grids_box_config_in)

    if (ny==0 .and. naky==0) call mp_abort("ERROR: ny==0 .and. naky==0", .true.) 
    if (nx==0 .and. ntheta0==0) call mp_abort("ERROR: nx==0 .and. ntheta0==0", .true.) 

    if (rhostar_box .gt. 0.0 .and. n0 .gt. 0) y0=1.0/(n0*rhostar_box*drhodpsi)

    if (gryfx) then
      !WHEN RUNNING IN GRYFX, ONLY EVOLVE ky=0 MODES.
      naky = 1
      ny = 1

      !NEED TO ACCOUNT FOR SQRT(2) DIFFERENCE BETWEEN rho_GS2 and rho_GryfX, so
      !change x0 and y0
      y0 = y0/sqrt(2.)  
      !this still needs to be renormalized even though only ky=0
      !is running because sometimes x0 is set from y0 and jtwist later.
      x0 = x0/sqrt(2.) 
      !this is in case x0 is not set from y0 and jtwist later.
    end if


    if (y0 < 0) y0 = -1./y0

    !EGH This line does not affect any existing gs2 runs but is
    !here because gryfx uses jtwist < 0 to signal using the
    ! default
    if (jtwist < 0) jtwist = max(int(2.0*pi*shat + 0.5),1)

    if (ly == 0.) ly = 2.0*pi*y0
    if (naky == 0) naky = (ny-1)/3 + 1
    if (ntheta0 == 0) ntheta0 = 2*((nx-1)/3) + 1
    if (rtwist == 0.) rtwist = real(jtwist)
    if (nkpolar == 0) nkpolar = int(real(naky-1.)*sqrt(2.))

    ! Now we make sure that we set ny and nx for given 
    ! choices of naky and ntheta0. If e.g. both ny and naky
    ! are set and they are not consistent with each other
    ! raise an error.

    if (ny == 0) then 
      if (proc0) write (error_unit(), *) "INFO: ny set from naky"
      ny = (naky - 1)*  3 + 1
    else if (naky /= (ny-1)/3 + 1) then
      if (proc0) write (error_unit(), *) "ERROR: naky and ny both set and inconsistent... set one or the other"
      call mp_abort("ERROR: naky and ny both set and inconsistent... set one or the other")
    end if
    if (nx == 0) then 
      if (proc0) write (error_unit(), *) "INFO: nx set from ntheta0"
      if (mod(ntheta0, 2) /= 1) then
        call mp_abort("ERROR: ntheta0 must be an odd number in box mode", .true.)
      end if
      nx = ((ntheta0 - 1) /  2) * 3 + 1
    else if (ntheta0 /= 2*((nx-1)/3) + 1) then
      if (proc0) write (error_unit(), *) "ERROR: ntheta0 and nx both set and inconsistent... set one or the other"
      call mp_abort("")
    end if

    
    nkpolar_private = nkpolar
    naky_private = naky
    ntheta0_private = ntheta0
    nx_private = nx
    ny_private = ny
  end subroutine init_kt_grids_box

  !> FIXME : Add documentation  
  subroutine wnml_kt_grids_box (unit)
    implicit none
    integer, intent(in) :: unit

    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_box_parameters"
    write (unit, fmt="(' nx = ',i4)") nx_private
    write (unit, fmt="(' ny = ',i4)") ny_private
    write (unit, fmt="(' Ly = ',e17.10)") ly
    if (rtwist /= 0.) then
       write (unit, fmt="(' rtwist = ',e17.10)") rtwist
    else
       write (unit, fmt="(' jtwist = ',i4)") jtwist
    end if
    write (unit, fmt="(' /')")
  end subroutine wnml_kt_grids_box

  !> FIXME : Add documentation  
  subroutine box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
    implicit none
    integer, intent (out) :: naky, ntheta0, nx, ny, nkpolar
    naky = naky_private
    ntheta0 = ntheta0_private
    nx = nx_private
    ny = ny_private
    nkpolar = nkpolar_private
  end subroutine box_get_sizes

  !> FIXME : Add documentation  
  subroutine box_set_overrides(grids_ov)
    use overrides, only: kt_grids_overrides_type
    type(kt_grids_overrides_type), intent(in) :: grids_ov
    if (grids_ov%override_naky) naky = grids_ov%naky
    if (grids_ov%override_ny) ny = grids_ov%ny
    if (grids_ov%override_ntheta0) ntheta0 = grids_ov%ntheta0
    if (grids_ov%override_nx) nx = grids_ov%nx
    if (grids_ov%override_x0) x0 = grids_ov%x0
    if (grids_ov%override_y0) y0 = grids_ov%y0
    if (grids_ov%override_gryfx) gryfx = grids_ov%gryfx
  end subroutine box_set_overrides

  !> FIXME : Add documentation  
  subroutine box_get_grids (aky, theta0, akx, ikx, iky)
    use theta_grid, only: shat
    use constants, only: pi
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    integer, dimension (:), intent (out) :: ikx, iky

    real :: dkx, dky, ratio
    integer :: i, naky, ntheta0

    naky = size(aky)    
    ntheta0 = size(akx)

    dky = 2.0*pi/ly

    if(abs(shat) <=  1.e-5) then   ! non-quantized b/c assumed to be periodic instead linked boundary conditions       

       if (x0 == 0.) then          
          
          if (rtwist > 0) then 
             ratio = rtwist
          else 
             ratio = 1. / abs(rtwist)
          end if
          
          dkx = dky / ratio
          
       else

          if (x0 > 0.) then
             dkx = 1./x0
          else
             dkx = -x0
          end if

       end if

    else
       if (jtwist /= 0) then
          dkx = dky * 2.0*pi*abs(shat)/real(jtwist)
       else
          dkx = dky
       end if
    endif

    do i = 1, naky
       iky(i) = i-1
       aky(i) = real(i-1)*dky
    end do

    do i = 1, (ntheta0+1)/2
       ikx(i) = i-1
       akx(i) = real(i-1)*dkx
    end do

    do i = (ntheta0+1)/2+1, ntheta0
       ikx(i) = i-ntheta0-1
       akx(i) = real(i-ntheta0-1)*dkx
    end do

    if (shat /= 0.) then
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = akx(i)/(aky(2:)*shat)
       end do
    else
       do i = 1, ntheta0
          theta0(i,1) = 0.0
          theta0(i,2:) = - akx(i)/aky(2:)   ! meaningless, so be careful
       end do
    end if

  end subroutine box_get_grids

  !> FIXME : Add documentation  
  subroutine check_kt_grids_box(report_unit)
    use theta_grid, only: shat_real => shat
    use constants, only: pi
    implicit none
    integer, intent(in) :: report_unit
    real :: lx, shat
    integer :: nx, ny, naky, ntheta0

    naky=naky_private
    ntheta0=ntheta0_private
    nx=nx_private
    ny=ny_private
    shat=shat_real

    if (y0 /= 2.) then
       if (abs(2.0*pi*y0 - ly) > 1.0e-7) then
          write (report_unit, *) 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, fmt="('You cannot specify both ly and y0.')")
          write (report_unit, fmt="(' ly=',e12.4,'  2.0*pi*y0=',2e12.4)") ly
          write (report_unit, fmt="('THIS IS AN ERROR.')") 
          write (report_unit, fmt="('################# WARNING #######################')")
          write (report_unit, *) 
       end if
    end if

    write (report_unit, *) 
    write (report_unit, fmt="('A rectangular simulation domain has been selected.')")
    if (rhostar_box .gt. 0.0 .and. n0 .gt. 0) write (report_unit, fmt="('The flux-tube size corresponds to toroidal mode number, n0=',i8/T44,' at rhostar_box=',1pe12.4)") n0,rhostar_box
    write (report_unit, *) 
    write (report_unit, fmt="('The domain is ',f10.4,' rho in the y direction.')") ly

    
    if (abs(shat) <= 1.e-5) then
       if (x0 == 0.) then
          if (rtwist > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  ly*rtwist
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -ly/rtwist
          end if
       else
          if (x0 > 0) then
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)")  2.*pi*x0
          else
             write (report_unit, fmt="('At theta=0, the domain has Lx = ',f10.5)") -2.*pi/x0
          end if
       end if
    else
       lx = ly * jtwist / (2.*pi*abs(shat))
       write (report_unit, fmt="('At theta=0, the domain is ',f10.4,' rho in the x direction.')") lx
       write (report_unit,*) ly, rtwist, jtwist, pi, shat
    end if
    
    write (report_unit, *) 
    write (report_unit, fmt="('The nonlinear terms will be evaluated on a grid with ',&
         & i4,' points in x and ',i4,' points in y.')") nx, ny
    write (report_unit, *) 
    naky = (ny-1)/3+1
    ntheta0 = 2*((nx-1)/3)+1
    write (report_unit, fmt="('After de-aliasing, there will be ',i4,'  ky >= 0 modes and ',i4,' kx modes.')") naky, ntheta0
    write (report_unit, fmt="('The modes with ky < 0 are determined by the reality condition.')")
  end subroutine check_kt_grids_box

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_box_config(kt_grids_box_config_in)
    use mp, only: mp_abort
    type(kt_grids_box_config_type), intent(in), optional :: kt_grids_box_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_box_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_box_config_in)) then
       kt_grids_box_config = kt_grids_box_config_in
    end if
  end subroutine set_kt_grids_box_config
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the kt_grids_box_parameters namelist and populates the member variables
  subroutine read_kt_grids_box_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_box_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    integer :: jtwist, n0, naky, nkpolar, ntheta0, nx, ny
    logical :: gryfx
    real :: ly, rhostar_box, rtwist, x0, y0

    namelist /kt_grids_box_parameters/ gryfx, jtwist, ly, n0, naky, nkpolar, ntheta0, nx, ny, rhostar_box, rtwist, x0, y0

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    gryfx = self%gryfx
    jtwist = self%jtwist
    ly = self%ly
    n0 = self%n0
    naky = self%naky
    nkpolar = self%nkpolar
    ntheta0 = self%ntheta0
    nx = self%nx
    ny = self%ny
    rhostar_box = self%rhostar_box
    rtwist = self%rtwist
    x0 = self%x0
    y0 = self%y0

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_box_parameters, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_box_parameters)
    endif

    ! Now copy from local variables into type members
    self%gryfx = gryfx
    self%jtwist = jtwist
    self%ly = ly
    self%n0 = n0
    self%naky = naky
    self%nkpolar = nkpolar
    self%ntheta0 = ntheta0
    self%nx = nx
    self%ny = ny
    self%rhostar_box = rhostar_box
    self%rtwist = rtwist
    self%x0 = x0
    self%y0 = y0

    self%exist = exist
  end subroutine read_kt_grids_box_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_box_config(self, unit)
    implicit none
    class(kt_grids_box_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("gryfx", self%gryfx, unit_internal)
    call self%write_key_val("jtwist", self%jtwist, unit_internal)
    call self%write_key_val("ly", self%ly, unit_internal)
    call self%write_key_val("n0", self%n0, unit_internal)
    call self%write_key_val("naky", self%naky, unit_internal)
    call self%write_key_val("nkpolar", self%nkpolar, unit_internal)
    call self%write_key_val("ntheta0", self%ntheta0, unit_internal)
    call self%write_key_val("nx", self%nx, unit_internal)
    call self%write_key_val("ny", self%ny, unit_internal)
    call self%write_key_val("rhostar_box", self%rhostar_box, unit_internal)
    call self%write_key_val("rtwist", self%rtwist, unit_internal)
    call self%write_key_val("x0", self%x0, unit_internal)
    call self%write_key_val("y0", self%y0, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_box_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_box_config(self)
    class(kt_grids_box_config_type), intent(in out) :: self
    type(kt_grids_box_config_type) :: empty
    select type (self)
    type is (kt_grids_box_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_box_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_box_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_box_config_type), intent(in out) :: self
    call broadcast(self%gryfx)
    call broadcast(self%jtwist)
    call broadcast(self%ly)
    call broadcast(self%n0)
    call broadcast(self%naky)
    call broadcast(self%nkpolar)
    call broadcast(self%ntheta0)
    call broadcast(self%nx)
    call broadcast(self%ny)
    call broadcast(self%rhostar_box)
    call broadcast(self%rtwist)
    call broadcast(self%x0)
    call broadcast(self%y0)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_box_config

  !> Gets the default name for this namelist
  function get_default_name_kt_grids_box_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_box_config
    get_default_name_kt_grids_box_config = "kt_grids_box_parameters"
  end function get_default_name_kt_grids_box_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_box_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_box_config
    get_default_requires_index_kt_grids_box_config = .false.
  end function get_default_requires_index_kt_grids_box_config
end module kt_grids_box

!> Set up the perpendicular wavenumbers by calling the appropriate sub-modules. 
module kt_grids
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use kt_grids_box, only: jtwist
  
  implicit none

  private

  public :: init_kt_grids, box, finish_kt_grids, check_kt_grids, wnml_kt
  public :: init_kt_grids_parameters, finish_kt_grids_parameters
  public :: set_overrides
  public :: aky, theta0, akx
  public :: naky, ntheta0, nx, ny, reality
  public :: nkpolar 
  public :: ikx, iky, jtwist_out
  public :: gridopt_switch, grid_option
  public :: gridopt_single, gridopt_range, gridopt_specified, gridopt_box
  public :: kwork_filter, kperp2
  public :: gryfx

  public :: kt_grids_config_type
  public :: set_kt_grids_config
  
  logical, dimension(:,:), allocatable :: kwork_filter
  real, dimension (:,:,:), allocatable :: kperp2
  real, dimension (:,:), allocatable :: theta0
  real, dimension (:), allocatable :: aky, akx
  integer, dimension(:), allocatable :: ikx, iky
  integer :: naky, ntheta0, nx, ny, nkpolar, jtwist_out
  character(20) :: grid_option

  ! internal variables
  integer :: gridopt_switch
  integer, parameter :: gridopt_single = 1, gridopt_range = 2, &
       gridopt_specified = 3, gridopt_box = 4
  logical :: reality = .false.
  logical :: box = .false.
  logical :: initialized = .false.
  logical :: kp2init=.false.
  logical :: nml_exist
  logical :: parameters_read = .false.
  
  !> Used to represent the input configuration of kt_grids
  type, extends(abstract_config_type) :: kt_grids_config_type
     ! namelist : kt_grids_knobs
     ! indexed : false     
     !> Controls the type of perpendicular wavenumber grid to use.
     !> Can be one of
     !>
     !> - 'single' Evolve a single Fourier mode. Set up
     !> with [[kt_grids_single_parameters]].
     !> - 'default' Same as 'single'
     !> - 'range' Evolve a range of equally spaced Fourier
     !> modes. Set up with [[kt_grids_range_parameters]].
     !> - 'specified' Evolve an arbitrary set of Fourier modes. Set
     !> up with [[kt_grids_specified_parameters]] and corresponding
     !> [[kt_grids_specified_elements]].
     !> - 'box' Simulation domain is logically rectangular. Set up with
     !> [[kt_grids_box_parameters]]. Required for nonlinear simulations.
     !> - 'nonlinear' Same as 'box.'
     !> - 'xbox' Experimental.
     !>
     character(len = 20) :: grid_option = "default"
   contains
     procedure, public :: read => read_kt_grids_config
     procedure, public :: write => write_kt_grids_config
     procedure, public :: reset => reset_kt_grids_config
     procedure, public :: broadcast => broadcast_kt_grids_config
     procedure, public, nopass :: get_default_name => get_default_name_kt_grids_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_kt_grids_config
  end type kt_grids_config_type
  
  type(kt_grids_config_type) :: kt_grids_config
  
contains

  !> FIXME : Add documentation
  subroutine init_kt_grids_parameters(kt_grids_config_in)
    use theta_grid, only: init_theta_grid
    use kt_grids_single, only: read_parameters_single
    use kt_grids_range, only: read_parameters_range
    use kt_grids_box, only: read_parameters_box
    implicit none
    type(kt_grids_config_type), intent(in), optional :: kt_grids_config_in
    
    if (parameters_read) return
    parameters_read = .true.
    
    call init_theta_grid
    nkpolar = 0   ! will be set to non-zero value only in box case; only used for an MHD diagnostic
    
    call read_parameters_internal(kt_grids_config_in)

    ! Read all namelists in case grid_option is overriden
    ! kt_grids_specified parameters cannot be overriden at this time
    ! so kt_grids_specified parameters are read in get_sizes
    !
    ! @note Should pass in optional configs here
    call read_parameters_single
    call read_parameters_range
    call read_parameters_box

  end subroutine init_kt_grids_parameters

  !> FIXME : Add documentation  
  subroutine finish_kt_grids_parameters
    use kt_grids_single, only: finish_parameters_single
    use kt_grids_range, only: finish_parameters_range
    use kt_grids_specified, only: finish_parameters_specified
    use kt_grids_box, only: finish_parameters_box
    implicit none

    call finish_parameters_single
    call finish_parameters_specified    
    call finish_parameters_range
    call finish_parameters_box

    call kt_grids_config%reset()
    parameters_read = .false.
  end subroutine finish_kt_grids_parameters

  !> FIXME : Add documentation  
  subroutine set_overrides(grids_ov)
    use overrides, only: kt_grids_overrides_type
    use kt_grids_box, only: box_set_overrides
    use mp, only: mp_abort, proc0
    implicit none
    type(kt_grids_overrides_type), intent(in) :: grids_ov
    
    if (proc0) then
      select case (gridopt_switch)          
      case (gridopt_box)
        call box_set_overrides(grids_ov)
      case default
        call mp_abort("Overrides currently only implemented for kt_grids_box", .true.)
      end select
    end if
  end subroutine set_overrides

  !> FIXME : Add documentation  
  function gryfx()
    use mp, only: broadcast
    use kt_grids_box, only: grbx=>gryfx
    implicit none
    logical :: gryfx
    gryfx = grbx
    !call broadcast(gryfx)
  end function gryfx

  !> FIXME : Add documentation
  !>
  !> @todo : Provide pass-through of other kt_grid module configs
  subroutine init_kt_grids(kt_grids_config_in)
    use mp, only: proc0, broadcast
    implicit none
    type(kt_grids_config_type), intent(in), optional :: kt_grids_config_in    

    if (initialized) return
    initialized = .true.

    call init_kt_grids_parameters(kt_grids_config_in)

    if (proc0) then
       call get_sizes
       jtwist_out = jtwist
    end if

    call broadcast (reality)
    call broadcast (box)
    call broadcast (naky)
    call broadcast (nkpolar)
    call broadcast (ntheta0)
    call broadcast (ny)
    call broadcast (nx)
    call broadcast (gridopt_switch)
    call allocate_arrays

    if (proc0) call get_grids
    call broadcast (ikx)     ! MR
    call broadcast (aky)
    call broadcast (akx)
    call broadcast (jtwist_out)
    call broadcast (theta0)
    allocate(kwork_filter(ntheta0,naky))
    kwork_filter=.false.
    call init_kperp2
  end subroutine init_kt_grids

  !> FIXME : Add documentation  
  subroutine read_parameters_internal(kt_grids_config_in)
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    implicit none
    type(kt_grids_config_type), intent(in), optional :: kt_grids_config_in
    type (text_option), dimension (6), parameter :: gridopts = &
         (/ text_option('default', gridopt_single), &
            text_option('single', gridopt_single), &
            text_option('range', gridopt_range), &
            text_option('specified', gridopt_specified), &
            text_option('box', gridopt_box), &
            text_option('nonlinear', gridopt_box) /)
    integer :: ierr

    if (present(kt_grids_config_in)) kt_grids_config = kt_grids_config_in

    call kt_grids_config%init(name = 'kt_grids_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    grid_option = kt_grids_config%grid_option
    
    nml_exist = kt_grids_config%exist
    
    ierr = error_unit()
    call get_option_value (grid_option, gridopts, gridopt_switch, &
         ierr, "grid_option in kt_grids_knobs",.true.)

  end subroutine read_parameters_internal

  !> FIXME : Add documentation  
  subroutine wnml_kt(unit)
    use kt_grids_single, only: wnml_kt_grids_single
    use kt_grids_range, only: wnml_kt_grids_range
    use kt_grids_specified, only: wnml_kt_grids_specified
    use kt_grids_box, only: wnml_kt_grids_box
    implicit none
    integer, intent(in) :: unit

    if (.not. nml_exist) return
    write (unit, *)
    write (unit, fmt="(' &',a)") "kt_grids_knobs"
    select case (gridopt_switch)          
    case (gridopt_single)
       write (unit, fmt="(' grid_option = ',a)") '"single"'
    case (gridopt_range)
       write (unit, fmt="(' grid_option = ',a)") '"range"'
    case (gridopt_specified)
       write (unit, fmt="(' grid_option = ',a)") '"specified"'
    case (gridopt_box)
       write (unit, fmt="(' grid_option = ',a)") '"box"'
    end select

    select case (gridopt_switch)          
    case (gridopt_single)
       call wnml_kt_grids_single (unit)
    case (gridopt_range)
       call wnml_kt_grids_range (unit)
    case (gridopt_specified)
       call wnml_kt_grids_specified (unit, aky, theta0)
    case (gridopt_box)
       call wnml_kt_grids_box (unit)
    end select
  end subroutine wnml_kt

  !> FIXME : Add documentation  
  subroutine allocate_arrays
    implicit none
    allocate (akx(ntheta0))
    allocate (aky(naky))
    allocate (theta0(ntheta0,naky))
    allocate (ikx(ntheta0))
    allocate (iky(naky))
  end subroutine allocate_arrays

  !> FIXME : Add documentation  
  subroutine get_sizes
    use kt_grids_single, only: init_kt_grids_single, single_get_sizes
    use kt_grids_range, only: init_kt_grids_range, range_get_sizes
    use kt_grids_specified, only: init_kt_grids_specified, specified_get_sizes
    use kt_grids_box, only: init_kt_grids_box, box_get_sizes
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call init_kt_grids_single
       call single_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_range)
       call init_kt_grids_range
       call range_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_specified)
       call init_kt_grids_specified
       call specified_get_sizes (naky, ntheta0, nx, ny)
    case (gridopt_box)
       call init_kt_grids_box
       call box_get_sizes (naky, ntheta0, nx, ny, nkpolar)
       reality = .true.
       box = .true.
    end select
  end subroutine get_sizes

  !> FIXME : Add documentation  
  subroutine get_grids
    use kt_grids_single, only: single_get_grids
    use kt_grids_range, only: range_get_grids
    use kt_grids_specified, only: specified_get_grids
    use kt_grids_box, only: box_get_grids
    implicit none
    select case (gridopt_switch)
    case (gridopt_single)
       call single_get_grids (aky, theta0, akx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx, ikx, iky)
    end select
  end subroutine get_grids

  !> FIXME : Add documentation  
  subroutine init_kperp2
    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat
    implicit none
    integer :: ik, it

    if (kp2init) return
    kp2init = .true.

    allocate (kperp2(-ntgrid:ntgrid,ntheta0,naky))
    do ik = 1, naky
       if (aky(ik) == 0.0) then
         do it = 1, ntheta0
             kperp2(:,it,ik) = akx(it)*akx(it)*gds22/(shat*shat)
          end do
       else
          do it = 1, ntheta0
             kperp2(:,it,ik) = aky(ik)*aky(ik) &
                  *(gds2 + 2.0*theta0(it,ik)*gds21 &
                  + theta0(it,ik)*theta0(it,ik)*gds22)
          end do
       end if
    end do
  end subroutine init_kperp2

  !> FIXME : Add documentation  
  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (akx, aky, theta0, ikx, iky)
    if (allocated(kwork_filter)) deallocate(kwork_filter)
    if (allocated(kperp2)) deallocate(kperp2)
    reality = .false. ; box = .false.
    kp2init = .false.
    initialized = .false.
    call finish_kt_grids_parameters
  end subroutine finish_kt_grids

  !> FIXME : Add documentation  
  subroutine check_kt_grids(report_unit)
!CMR, 22/9/2010:
! add routine to GS2 to perform ingen checking and reporting
    use kt_grids_single, only: check_kt_grids_single
    use kt_grids_range,  only: check_kt_grids_range
    use kt_grids_box,    only: check_kt_grids_box
    use kt_grids_specified, only: check_kt_grids_specified
    implicit none
    integer, intent(in) :: report_unit

    select case (gridopt_switch) 
    case (gridopt_single) 
       call check_kt_grids_single (report_unit)
    case (gridopt_range)
       call check_kt_grids_range (report_unit)
    case (gridopt_specified) 
       call check_kt_grids_specified (report_unit, aky, theta0(:,1))
    case (gridopt_box)
       call check_kt_grids_box (report_unit)
    end select
  end subroutine check_kt_grids

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_kt_grids_config(kt_grids_config_in)
    use mp, only: mp_abort
    type(kt_grids_config_type), intent(in), optional :: kt_grids_config_in
    if (initialized) then
       call mp_abort("Trying to set kt_grids_config when already initialized.", to_screen = .true.)
    end if
    if (present(kt_grids_config_in)) then
       kt_grids_config = kt_grids_config_in
    end if
  end subroutine set_kt_grids_config

  !---------------------------------------
  ! Following is for the config_type
  !--------------------------------------- 
  
  !> Reads in the kt_grids_knobs namelist and populates the member variables
  subroutine read_kt_grids_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(kt_grids_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    character(len = 20) :: grid_option

    namelist /kt_grids_knobs/ grid_option

    ! Only proc0 reads from file
    if (.not. proc0) return

    ! First set local variables from current values
    grid_option = self%grid_option

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = kt_grids_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)       
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = kt_grids_knobs)
    endif

    ! Now copy from local variables into type members
    self%grid_option = grid_option

    self%exist = exist
  end subroutine read_kt_grids_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_kt_grids_config(self, unit)
    implicit none
    class(kt_grids_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
    integer :: unit_internal

    unit_internal = 6 ! @todo -- get stdout from file_utils
    if (present(unit)) then
       unit_internal = unit
    endif

    call self%write_namelist_header(unit_internal)
    call self%write_key_val("grid_option", self%grid_option, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_kt_grids_config

  !> Resets the config object to the initial empty state
  subroutine reset_kt_grids_config(self)
    class(kt_grids_config_type), intent(in out) :: self
    type(kt_grids_config_type) :: empty
    select type (self)
    type is (kt_grids_config_type)
       self = empty
    end select
  end subroutine reset_kt_grids_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_kt_grids_config(self)
    use mp, only: broadcast
    implicit none
    class(kt_grids_config_type), intent(in out) :: self
    call broadcast(self%grid_option)

    call broadcast(self%exist)
  end subroutine broadcast_kt_grids_config

  !> Gets the default name for this namelist
  function get_default_name_kt_grids_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_config
    get_default_name_kt_grids_config = "kt_grids_knobs"
  end function get_default_name_kt_grids_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_kt_grids_config()
    implicit none
    logical :: get_default_requires_index_kt_grids_config
    get_default_requires_index_kt_grids_config = .false.
  end function get_default_requires_index_kt_grids_config
end module kt_grids

