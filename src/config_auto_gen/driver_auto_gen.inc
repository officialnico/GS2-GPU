
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the driver namelist and populates the member variables
subroutine read_driver_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(driver_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  complex :: w_antenna
  integer :: nk_stir
  logical :: ant_off, restarting, write_antenna
  real :: amplitude, t0, w_dot

  namelist /driver/ amplitude, ant_off, nk_stir, restarting, t0, w_antenna, w_dot, write_antenna

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  amplitude = self%amplitude
  ant_off = self%ant_off
  nk_stir = self%nk_stir
  restarting = self%restarting
  t0 = self%t0
  w_antenna = self%w_antenna
  w_dot = self%w_dot
  write_antenna = self%write_antenna

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = driver)

  ! Now copy from local variables into type members
  self%amplitude = amplitude
  self%ant_off = ant_off
  self%nk_stir = nk_stir
  self%restarting = restarting
  self%t0 = t0
  self%w_antenna = w_antenna
  self%w_dot = w_dot
  self%write_antenna = write_antenna

  self%exist = exist
end subroutine read_driver_config

!> Writes out a namelist representing the current state of the config object
subroutine write_driver_config(self, unit)
  implicit none
  class(driver_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("amplitude", self%amplitude, unit_internal)
  call self%write_key_val("ant_off", self%ant_off, unit_internal)
  call self%write_key_val("nk_stir", self%nk_stir, unit_internal)
  call self%write_key_val("restarting", self%restarting, unit_internal)
  call self%write_key_val("t0", self%t0, unit_internal)
  call self%write_key_val("w_antenna", self%w_antenna, unit_internal)
  call self%write_key_val("w_dot", self%w_dot, unit_internal)
  call self%write_key_val("write_antenna", self%write_antenna, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_driver_config

!> Resets the config object to the initial empty state
subroutine reset_driver_config(self)
  class(driver_config_type), intent(in out) :: self
  type(driver_config_type) :: empty
  select type (self)
  type is (driver_config_type)
     self = empty
  end select
end subroutine reset_driver_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_driver_config(self)
  use mp, only: broadcast
  implicit none
  class(driver_config_type), intent(in out) :: self
  call broadcast(self%amplitude)
  call broadcast(self%ant_off)
  call broadcast(self%nk_stir)
  call broadcast(self%restarting)
  call broadcast(self%t0)
  call broadcast(self%w_antenna)
  call broadcast(self%w_dot)
  call broadcast(self%write_antenna)

  call broadcast(self%exist)
end subroutine broadcast_driver_config

!> Gets the default name for this namelist
function get_default_name_driver_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_driver_config
  get_default_name_driver_config = "driver"
end function get_default_name_driver_config

!> Gets the default requires index for this namelist
function get_default_requires_index_driver_config()
  implicit none
  logical :: get_default_requires_index_driver_config
  get_default_requires_index_driver_config = .false.
end function get_default_requires_index_driver_config

!> Get the module level config instance
pure function get_driver_config()
  type(driver_config_type) :: get_driver_config
  get_driver_config = driver_config
end function get_driver_config
