
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the fields_knobs namelist and populates the member variables
subroutine read_fields_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(fields_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: field_option
  character(len = 256) :: response_dir, response_file
  integer :: minnrow
  logical :: do_smart_update, dump_response, field_local_allreduce, field_local_allreduce_sub, field_local_nonblocking_collectives, field_local_tuneminnrow, field_subgath, force_maxwell_reinit, read_response
  logical :: remove_zonal_flows_switch

  namelist /fields_knobs/ do_smart_update, dump_response, field_local_allreduce, field_local_allreduce_sub, field_local_nonblocking_collectives, field_local_tuneminnrow, field_option, field_subgath, force_maxwell_reinit, &
    minnrow, read_response, remove_zonal_flows_switch, response_dir, response_file

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  do_smart_update = self%do_smart_update
  dump_response = self%dump_response
  field_local_allreduce = self%field_local_allreduce
  field_local_allreduce_sub = self%field_local_allreduce_sub
  field_local_nonblocking_collectives = self%field_local_nonblocking_collectives
  field_local_tuneminnrow = self%field_local_tuneminnrow
  field_option = self%field_option
  field_subgath = self%field_subgath
  force_maxwell_reinit = self%force_maxwell_reinit
  minnrow = self%minnrow
  read_response = self%read_response
  remove_zonal_flows_switch = self%remove_zonal_flows_switch
  response_dir = self%response_dir
  response_file = self%response_file

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = fields_knobs)

  ! Now copy from local variables into type members
  self%do_smart_update = do_smart_update
  self%dump_response = dump_response
  self%field_local_allreduce = field_local_allreduce
  self%field_local_allreduce_sub = field_local_allreduce_sub
  self%field_local_nonblocking_collectives = field_local_nonblocking_collectives
  self%field_local_tuneminnrow = field_local_tuneminnrow
  self%field_option = field_option
  self%field_subgath = field_subgath
  self%force_maxwell_reinit = force_maxwell_reinit
  self%minnrow = minnrow
  self%read_response = read_response
  self%remove_zonal_flows_switch = remove_zonal_flows_switch
  self%response_dir = response_dir
  self%response_file = response_file

  self%exist = exist
end subroutine read_fields_config

!> Writes out a namelist representing the current state of the config object
subroutine write_fields_config(self, unit)
  implicit none
  class(fields_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("do_smart_update", self%do_smart_update, unit_internal)
  call self%write_key_val("dump_response", self%dump_response, unit_internal)
  call self%write_key_val("field_local_allreduce", self%field_local_allreduce, unit_internal)
  call self%write_key_val("field_local_allreduce_sub", self%field_local_allreduce_sub, unit_internal)
  call self%write_key_val("field_local_nonblocking_collectives", self%field_local_nonblocking_collectives, unit_internal)
  call self%write_key_val("field_local_tuneminnrow", self%field_local_tuneminnrow, unit_internal)
  call self%write_key_val("field_option", self%field_option, unit_internal)
  call self%write_key_val("field_subgath", self%field_subgath, unit_internal)
  call self%write_key_val("force_maxwell_reinit", self%force_maxwell_reinit, unit_internal)
  call self%write_key_val("minnrow", self%minnrow, unit_internal)
  call self%write_key_val("read_response", self%read_response, unit_internal)
  call self%write_key_val("remove_zonal_flows_switch", self%remove_zonal_flows_switch, unit_internal)
  call self%write_key_val("response_dir", self%response_dir, unit_internal)
  call self%write_key_val("response_file", self%response_file, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_fields_config

!> Resets the config object to the initial empty state
subroutine reset_fields_config(self)
  class(fields_config_type), intent(in out) :: self
  type(fields_config_type) :: empty
  select type (self)
  type is (fields_config_type)
     self = empty
  end select
end subroutine reset_fields_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_fields_config(self)
  use mp, only: broadcast
  implicit none
  class(fields_config_type), intent(in out) :: self
  call broadcast(self%do_smart_update)
  call broadcast(self%dump_response)
  call broadcast(self%field_local_allreduce)
  call broadcast(self%field_local_allreduce_sub)
  call broadcast(self%field_local_nonblocking_collectives)
  call broadcast(self%field_local_tuneminnrow)
  call broadcast(self%field_option)
  call broadcast(self%field_subgath)
  call broadcast(self%force_maxwell_reinit)
  call broadcast(self%minnrow)
  call broadcast(self%read_response)
  call broadcast(self%remove_zonal_flows_switch)
  call broadcast(self%response_dir)
  call broadcast(self%response_file)

  call broadcast(self%exist)
end subroutine broadcast_fields_config

!> Gets the default name for this namelist
function get_default_name_fields_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_fields_config
  get_default_name_fields_config = "fields_knobs"
end function get_default_name_fields_config

!> Gets the default requires index for this namelist
function get_default_requires_index_fields_config()
  implicit none
  logical :: get_default_requires_index_fields_config
  get_default_requires_index_fields_config = .false.
end function get_default_requires_index_fields_config

!> Get the module level config instance
pure function get_fields_config()
  type(fields_config_type) :: get_fields_config
  get_fields_config = fields_config
end function get_fields_config
