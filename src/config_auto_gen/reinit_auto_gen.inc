
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the reinit_knobs namelist and populates the member variables
subroutine read_reinit_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(reinit_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  logical :: abort_rapid_time_step_change, in_memory
  real :: delt_adj, delt_cushion, delt_minimum, dt0

  namelist /reinit_knobs/ abort_rapid_time_step_change, delt_adj, delt_cushion, delt_minimum, dt0, in_memory

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  abort_rapid_time_step_change = self%abort_rapid_time_step_change
  delt_adj = self%delt_adj
  delt_cushion = self%delt_cushion
  delt_minimum = self%delt_minimum
  dt0 = self%dt0
  in_memory = self%in_memory

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = reinit_knobs)

  ! Now copy from local variables into type members
  self%abort_rapid_time_step_change = abort_rapid_time_step_change
  self%delt_adj = delt_adj
  self%delt_cushion = delt_cushion
  self%delt_minimum = delt_minimum
  self%dt0 = dt0
  self%in_memory = in_memory

  self%exist = exist
end subroutine read_reinit_config

!> Writes out a namelist representing the current state of the config object
subroutine write_reinit_config(self, unit)
  implicit none
  class(reinit_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("abort_rapid_time_step_change", self%abort_rapid_time_step_change, unit_internal)
  call self%write_key_val("delt_adj", self%delt_adj, unit_internal)
  call self%write_key_val("delt_cushion", self%delt_cushion, unit_internal)
  call self%write_key_val("delt_minimum", self%delt_minimum, unit_internal)
  call self%write_key_val("dt0", self%dt0, unit_internal)
  call self%write_key_val("in_memory", self%in_memory, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_reinit_config

!> Resets the config object to the initial empty state
subroutine reset_reinit_config(self)
  class(reinit_config_type), intent(in out) :: self
  type(reinit_config_type) :: empty
  select type (self)
  type is (reinit_config_type)
     self = empty
  end select
end subroutine reset_reinit_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_reinit_config(self)
  use mp, only: broadcast
  implicit none
  class(reinit_config_type), intent(in out) :: self
  call broadcast(self%abort_rapid_time_step_change)
  call broadcast(self%delt_adj)
  call broadcast(self%delt_cushion)
  call broadcast(self%delt_minimum)
  call broadcast(self%dt0)
  call broadcast(self%in_memory)

  call broadcast(self%exist)
end subroutine broadcast_reinit_config

!> Gets the default name for this namelist
function get_default_name_reinit_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_reinit_config
  get_default_name_reinit_config = "reinit_knobs"
end function get_default_name_reinit_config

!> Gets the default requires index for this namelist
function get_default_requires_index_reinit_config()
  implicit none
  logical :: get_default_requires_index_reinit_config
  get_default_requires_index_reinit_config = .false.
end function get_default_requires_index_reinit_config

!> Get the module level config instance
pure function get_reinit_config()
  type(reinit_config_type) :: get_reinit_config
  get_reinit_config = reinit_config
end function get_reinit_config
