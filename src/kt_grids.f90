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
  public :: ikx, jtwist_out
  public :: gridopt_switch, grid_option
  public :: gridopt_single, gridopt_range, gridopt_specified, gridopt_box
  public :: kwork_filter, kperp2, inv_kperp2

  public :: kt_grids_config_type
  public :: set_kt_grids_config
  public :: get_kt_grids_config

  logical, dimension(:,:), allocatable :: kwork_filter
  real, dimension (:,:,:), allocatable :: kperp2, inv_kperp2
  !> The \(\theta_0(k_x, k_y)\) grid
  real, dimension (:,:), allocatable :: theta0
  !> The \(k_x \rho\) grid
  real, dimension (:), allocatable :: akx
  !> The \(k_y \rho\) grid
  real, dimension (:), allocatable :: aky
  !> Discrete kx grid index
  integer, dimension(:), allocatable :: ikx
  !> Number of \(k_y \rho\) points
  integer :: naky
  !> Number of \(\theta_0\) points
  integer :: ntheta0
  !> Number of (real space) \(x\) points
  integer :: nx
  !> Number of (real space) \(y\) points
  integer :: ny
  !> Related to the box size in x. See
  !> [[kt_grids_box_parameters:jtwist]]
  integer :: jtwist_out
  !> The type of perpendicular wavenumber grid used. See
  !> [[kt_grids_knobs:grid_option]]
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
     !> [[kt_grids_specified_element]].
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
    write(unit, fmt="(' /')")

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
       call box_get_sizes (naky, ntheta0, nx, ny)
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
       call single_get_grids (aky, theta0, akx, ikx)
    case (gridopt_range)
       call range_get_grids (aky, theta0, akx, ikx)
    case (gridopt_specified)
       call specified_get_grids (aky, theta0, akx, ikx)
    case (gridopt_box)
       call box_get_grids (aky, theta0, akx, ikx)
    end select
  end subroutine get_grids

  !> FIXME : Add documentation  
  subroutine init_kperp2
    use theta_grid, only: ntgrid, gds2, gds21, gds22, shat
    implicit none
    integer :: ik, it, ig

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

    allocate (inv_kperp2(-ntgrid:ntgrid,ntheta0,naky))

    do ik = 1, naky
       do it = 1, ntheta0
          do ig = -ntgrid, ntgrid
             if (abs(kperp2(ig, it, ik)) > epsilon(0.0)) then
                inv_kperp2(ig, it, ik) = 1.0/kperp2(ig, it, ik)
             else
                inv_kperp2(ig, it, ik) = 0.0
             end if
          end do
       end do
    end do
  end subroutine init_kperp2

  !> FIXME : Add documentation  
  subroutine finish_kt_grids

    implicit none

    if (allocated(aky)) deallocate (akx, aky, theta0, ikx)
    if (allocated(kwork_filter)) deallocate(kwork_filter)
    if (allocated(kperp2)) deallocate(kperp2)
    if (allocated(inv_kperp2)) deallocate(inv_kperp2)
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

  !> Get the module level config instance
  function get_kt_grids_config()
    type(kt_grids_config_type) :: get_kt_grids_config
    get_kt_grids_config = kt_grids_config
  end function get_kt_grids_config

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
    integer :: in_file

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
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = kt_grids_knobs)

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
