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
  public :: get_kt_grids_specified_config
  public :: get_kt_grids_specified_element_config

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
  subroutine specified_get_grids (aky, theta0, akx, ikx, kt_grids_specified_element_config_in)
    use mp, only: proc0
    implicit none
    real, dimension (:), intent (out) :: akx, aky
    real, dimension (:,:), intent (out) :: theta0
    !> Discrete kx wavenumber grid indices
    integer, dimension (:), intent (out) :: ikx
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
       if (i <= ntheta0) ikx(i) = i - 1
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

  !> Get the module level config instance
  function get_kt_grids_specified_config()
    type(kt_grids_specified_config_type) :: get_kt_grids_specified_config
    get_kt_grids_specified_config = kt_grids_specified_config
  end function get_kt_grids_specified_config

  !> Get the module level config instance
  function get_kt_grids_specified_element_config()
    type(kt_grids_specified_element_config_type), allocatable :: get_kt_grids_specified_element_config(:)
    if (allocated(kt_grids_specified_element_config)) then
       get_kt_grids_specified_element_config = kt_grids_specified_element_config
    else
       allocate(get_kt_grids_specified_element_config(0))
    end if
  end function get_kt_grids_specified_element_config

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
    integer :: in_file

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
    in_file = input_unit_exist(self%get_name(), exist)
    if (exist) read(in_file, nml = kt_grids_specified_parameters)

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
    integer :: in_file

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
    call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index, exist)
    if (exist) read(in_file, nml = kt_grids_specified_element)
    close(unit = in_file)

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
