
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
  integer, intent(in), optional:: unit
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

!> Get the module level config instance
pure function get_kt_grids_specified_config()
  type(kt_grids_specified_config_type) :: get_kt_grids_specified_config
  get_kt_grids_specified_config = kt_grids_specified_config
end function get_kt_grids_specified_config
