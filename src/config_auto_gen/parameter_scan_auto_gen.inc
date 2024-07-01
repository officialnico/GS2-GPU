
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the parameter_scan_knobs namelist and populates the member variables
subroutine read_parameter_scan_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(parameter_scan_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: inc_con, scan_par, scan_type, target_par
  integer :: nstep_inc, nstep_init, scan_spec
  logical :: scan_restarted
  real :: delta_t_inc, delta_t_init, par_end, par_inc, par_start, target_val

  namelist /parameter_scan_knobs/ delta_t_inc, delta_t_init, inc_con, nstep_inc, nstep_init, par_end, par_inc, par_start, scan_par, scan_restarted, scan_spec, scan_type, target_par, target_val

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  delta_t_inc = self%delta_t_inc
  delta_t_init = self%delta_t_init
  inc_con = self%inc_con
  nstep_inc = self%nstep_inc
  nstep_init = self%nstep_init
  par_end = self%par_end
  par_inc = self%par_inc
  par_start = self%par_start
  scan_par = self%scan_par
  scan_restarted = self%scan_restarted
  scan_spec = self%scan_spec
  scan_type = self%scan_type
  target_par = self%target_par
  target_val = self%target_val

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = parameter_scan_knobs)

  ! Now copy from local variables into type members
  self%delta_t_inc = delta_t_inc
  self%delta_t_init = delta_t_init
  self%inc_con = inc_con
  self%nstep_inc = nstep_inc
  self%nstep_init = nstep_init
  self%par_end = par_end
  self%par_inc = par_inc
  self%par_start = par_start
  self%scan_par = scan_par
  self%scan_restarted = scan_restarted
  self%scan_spec = scan_spec
  self%scan_type = scan_type
  self%target_par = target_par
  self%target_val = target_val

  self%exist = exist
end subroutine read_parameter_scan_config

!> Writes out a namelist representing the current state of the config object
subroutine write_parameter_scan_config(self, unit)
  implicit none
  class(parameter_scan_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("delta_t_inc", self%delta_t_inc, unit_internal)
  call self%write_key_val("delta_t_init", self%delta_t_init, unit_internal)
  call self%write_key_val("inc_con", self%inc_con, unit_internal)
  call self%write_key_val("nstep_inc", self%nstep_inc, unit_internal)
  call self%write_key_val("nstep_init", self%nstep_init, unit_internal)
  call self%write_key_val("par_end", self%par_end, unit_internal)
  call self%write_key_val("par_inc", self%par_inc, unit_internal)
  call self%write_key_val("par_start", self%par_start, unit_internal)
  call self%write_key_val("scan_par", self%scan_par, unit_internal)
  call self%write_key_val("scan_restarted", self%scan_restarted, unit_internal)
  call self%write_key_val("scan_spec", self%scan_spec, unit_internal)
  call self%write_key_val("scan_type", self%scan_type, unit_internal)
  call self%write_key_val("target_par", self%target_par, unit_internal)
  call self%write_key_val("target_val", self%target_val, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_parameter_scan_config

!> Resets the config object to the initial empty state
subroutine reset_parameter_scan_config(self)
  class(parameter_scan_config_type), intent(in out) :: self
  type(parameter_scan_config_type) :: empty
  select type (self)
  type is (parameter_scan_config_type)
     self = empty
  end select
end subroutine reset_parameter_scan_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_parameter_scan_config(self)
  use mp, only: broadcast
  implicit none
  class(parameter_scan_config_type), intent(in out) :: self
  call broadcast(self%delta_t_inc)
  call broadcast(self%delta_t_init)
  call broadcast(self%inc_con)
  call broadcast(self%nstep_inc)
  call broadcast(self%nstep_init)
  call broadcast(self%par_end)
  call broadcast(self%par_inc)
  call broadcast(self%par_start)
  call broadcast(self%scan_par)
  call broadcast(self%scan_restarted)
  call broadcast(self%scan_spec)
  call broadcast(self%scan_type)
  call broadcast(self%target_par)
  call broadcast(self%target_val)

  call broadcast(self%exist)
end subroutine broadcast_parameter_scan_config

!> Gets the default name for this namelist
function get_default_name_parameter_scan_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_parameter_scan_config
  get_default_name_parameter_scan_config = "parameter_scan_knobs"
end function get_default_name_parameter_scan_config

!> Gets the default requires index for this namelist
function get_default_requires_index_parameter_scan_config()
  implicit none
  logical :: get_default_requires_index_parameter_scan_config
  get_default_requires_index_parameter_scan_config = .false.
end function get_default_requires_index_parameter_scan_config

!> Get the module level config instance
pure function get_parameter_scan_config()
  type(parameter_scan_config_type) :: get_parameter_scan_config
  get_parameter_scan_config = parameter_scan_config
end function get_parameter_scan_config
