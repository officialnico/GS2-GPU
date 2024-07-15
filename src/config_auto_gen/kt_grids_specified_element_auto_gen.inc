
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
  integer, intent(in), optional:: unit
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

!> Get the array of module level config instances. If it isn't allocated,
!> then return a zero-length array
pure function get_kt_grids_specified_element_config()
  type(kt_grids_specified_element_config_type), allocatable, dimension(:) :: get_kt_grids_specified_element_config
  if (allocated(kt_grids_specified_element_config)) then
    get_kt_grids_specified_element_config = kt_grids_specified_element_config
  else
    allocate(get_kt_grids_specified_element_config(0))
  end if
end function get_kt_grids_specified_element_config
