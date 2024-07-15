
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the ingen_knobs namelist and populates the member variables
subroutine read_ingen_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(ingen_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: ncut, npmax
  logical :: scan, stdin

  namelist /ingen_knobs/ ncut, npmax, scan, stdin

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  ncut = self%ncut
  npmax = self%npmax
  scan = self%scan
  stdin = self%stdin

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = ingen_knobs)

  ! Now copy from local variables into type members
  self%ncut = ncut
  self%npmax = npmax
  self%scan = scan
  self%stdin = stdin

  self%exist = exist
end subroutine read_ingen_config

!> Writes out a namelist representing the current state of the config object
subroutine write_ingen_config(self, unit)
  implicit none
  class(ingen_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("ncut", self%ncut, unit_internal)
  call self%write_key_val("npmax", self%npmax, unit_internal)
  call self%write_key_val("scan", self%scan, unit_internal)
  call self%write_key_val("stdin", self%stdin, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_ingen_config

!> Resets the config object to the initial empty state
subroutine reset_ingen_config(self)
  class(ingen_config_type), intent(in out) :: self
  type(ingen_config_type) :: empty
  select type (self)
  type is (ingen_config_type)
     self = empty
  end select
end subroutine reset_ingen_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_ingen_config(self)
  use mp, only: broadcast
  implicit none
  class(ingen_config_type), intent(in out) :: self
  call broadcast(self%ncut)
  call broadcast(self%npmax)
  call broadcast(self%scan)
  call broadcast(self%stdin)

  call broadcast(self%exist)
end subroutine broadcast_ingen_config

!> Gets the default name for this namelist
function get_default_name_ingen_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_ingen_config
  get_default_name_ingen_config = "ingen_knobs"
end function get_default_name_ingen_config

!> Gets the default requires index for this namelist
function get_default_requires_index_ingen_config()
  implicit none
  logical :: get_default_requires_index_ingen_config
  get_default_requires_index_ingen_config = .false.
end function get_default_requires_index_ingen_config

!> Get the module level config instance
pure function get_ingen_config()
  type(ingen_config_type) :: get_ingen_config
  get_ingen_config = ingen_config
end function get_ingen_config
