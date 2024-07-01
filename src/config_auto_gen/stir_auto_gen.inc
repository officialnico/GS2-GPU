
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the stir namelist and populates the member variables
subroutine read_stir_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(stir_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: kx, ky, kz
  logical :: travel
  real :: a, b

  namelist /stir/ a, b, kx, ky, kz, travel

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  a = self%a
  b = self%b
  kx = self%kx
  ky = self%ky
  kz = self%kz
  travel = self%travel

  ! Now read in the main namelist
  call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index, exist)
  if (exist) read(in_file, nml = stir)
  close(unit = in_file)

  ! Now copy from local variables into type members
  self%a = a
  self%b = b
  self%kx = kx
  self%ky = ky
  self%kz = kz
  self%travel = travel

  self%exist = exist
end subroutine read_stir_config

!> Writes out a namelist representing the current state of the config object
subroutine write_stir_config(self, unit)
  implicit none
  class(stir_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("a", self%a, unit_internal)
  call self%write_key_val("b", self%b, unit_internal)
  call self%write_key_val("kx", self%kx, unit_internal)
  call self%write_key_val("ky", self%ky, unit_internal)
  call self%write_key_val("kz", self%kz, unit_internal)
  call self%write_key_val("travel", self%travel, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_stir_config

!> Resets the config object to the initial empty state
subroutine reset_stir_config(self)
  class(stir_config_type), intent(in out) :: self
  type(stir_config_type) :: empty
  select type (self)
  type is (stir_config_type)
     self = empty
  end select
end subroutine reset_stir_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_stir_config(self)
  use mp, only: broadcast
  implicit none
  class(stir_config_type), intent(in out) :: self
  call broadcast(self%a)
  call broadcast(self%b)
  call broadcast(self%kx)
  call broadcast(self%ky)
  call broadcast(self%kz)
  call broadcast(self%travel)

  call broadcast(self%exist)
end subroutine broadcast_stir_config

!> Gets the default name for this namelist
function get_default_name_stir_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_stir_config
  get_default_name_stir_config = "stir"
end function get_default_name_stir_config

!> Gets the default requires index for this namelist
function get_default_requires_index_stir_config()
  implicit none
  logical :: get_default_requires_index_stir_config
  get_default_requires_index_stir_config = .true.
end function get_default_requires_index_stir_config

!> Get the array of module level config instances. If it isn't allocated,
!> then return a zero-length array
pure function get_stir_config()
  type(stir_config_type), allocatable, dimension(:) :: get_stir_config
  if (allocated(stir_config)) then
    get_stir_config = stir_config
  else
    allocate(get_stir_config(0))
  end if
end function get_stir_config
