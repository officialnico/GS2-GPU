
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the theta_grid_file_knobs namelist and populates the member variables
subroutine read_theta_grid_file_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(theta_grid_file_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = run_name_size) :: gridout_file
  logical :: no_geo_info

  namelist /theta_grid_file_knobs/ gridout_file, no_geo_info

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  gridout_file = self%gridout_file
  no_geo_info = self%no_geo_info

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = theta_grid_file_knobs)

  ! Now copy from local variables into type members
  self%gridout_file = gridout_file
  self%no_geo_info = no_geo_info

  self%exist = exist
end subroutine read_theta_grid_file_config

!> Writes out a namelist representing the current state of the config object
subroutine write_theta_grid_file_config(self, unit)
  implicit none
  class(theta_grid_file_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("gridout_file", self%gridout_file, unit_internal)
  call self%write_key_val("no_geo_info", self%no_geo_info, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_theta_grid_file_config

!> Resets the config object to the initial empty state
subroutine reset_theta_grid_file_config(self)
  class(theta_grid_file_config_type), intent(in out) :: self
  type(theta_grid_file_config_type) :: empty
  select type (self)
  type is (theta_grid_file_config_type)
     self = empty
  end select
end subroutine reset_theta_grid_file_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_theta_grid_file_config(self)
  use mp, only: broadcast
  implicit none
  class(theta_grid_file_config_type), intent(in out) :: self
  call broadcast(self%gridout_file)
  call broadcast(self%no_geo_info)

  call broadcast(self%exist)
end subroutine broadcast_theta_grid_file_config

!> Gets the default name for this namelist
function get_default_name_theta_grid_file_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_file_config
  get_default_name_theta_grid_file_config = "theta_grid_file_knobs"
end function get_default_name_theta_grid_file_config

!> Gets the default requires index for this namelist
function get_default_requires_index_theta_grid_file_config()
  implicit none
  logical :: get_default_requires_index_theta_grid_file_config
  get_default_requires_index_theta_grid_file_config = .false.
end function get_default_requires_index_theta_grid_file_config

!> Get the module level config instance
pure function get_theta_grid_file_config()
  type(theta_grid_file_config_type) :: get_theta_grid_file_config
  get_theta_grid_file_config = theta_grid_file_config
end function get_theta_grid_file_config
