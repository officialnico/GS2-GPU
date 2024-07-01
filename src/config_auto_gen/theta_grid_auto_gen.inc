
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the theta_grid_knobs namelist and populates the member variables
subroutine read_theta_grid_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(theta_grid_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: equilibrium_option
  logical :: gb_to_cv
  real :: cvdriftknob, gbdriftknob

  namelist /theta_grid_knobs/ cvdriftknob, equilibrium_option, gb_to_cv, gbdriftknob

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  cvdriftknob = self%cvdriftknob
  equilibrium_option = self%equilibrium_option
  gb_to_cv = self%gb_to_cv
  gbdriftknob = self%gbdriftknob

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = theta_grid_knobs)

  ! Now copy from local variables into type members
  self%cvdriftknob = cvdriftknob
  self%equilibrium_option = equilibrium_option
  self%gb_to_cv = gb_to_cv
  self%gbdriftknob = gbdriftknob

  self%exist = exist
end subroutine read_theta_grid_config

!> Writes out a namelist representing the current state of the config object
subroutine write_theta_grid_config(self, unit)
  implicit none
  class(theta_grid_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("cvdriftknob", self%cvdriftknob, unit_internal)
  call self%write_key_val("equilibrium_option", self%equilibrium_option, unit_internal)
  call self%write_key_val("gb_to_cv", self%gb_to_cv, unit_internal)
  call self%write_key_val("gbdriftknob", self%gbdriftknob, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_theta_grid_config

!> Resets the config object to the initial empty state
subroutine reset_theta_grid_config(self)
  class(theta_grid_config_type), intent(in out) :: self
  type(theta_grid_config_type) :: empty
  select type (self)
  type is (theta_grid_config_type)
     self = empty
  end select
end subroutine reset_theta_grid_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_theta_grid_config(self)
  use mp, only: broadcast
  implicit none
  class(theta_grid_config_type), intent(in out) :: self
  call broadcast(self%cvdriftknob)
  call broadcast(self%equilibrium_option)
  call broadcast(self%gb_to_cv)
  call broadcast(self%gbdriftknob)

  call broadcast(self%exist)
end subroutine broadcast_theta_grid_config

!> Gets the default name for this namelist
function get_default_name_theta_grid_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_config
  get_default_name_theta_grid_config = "theta_grid_knobs"
end function get_default_name_theta_grid_config

!> Gets the default requires index for this namelist
function get_default_requires_index_theta_grid_config()
  implicit none
  logical :: get_default_requires_index_theta_grid_config
  get_default_requires_index_theta_grid_config = .false.
end function get_default_requires_index_theta_grid_config

!> Get the module level config instance
pure function get_theta_grid_config()
  type(theta_grid_config_type) :: get_theta_grid_config
  get_theta_grid_config = theta_grid_config
end function get_theta_grid_config
