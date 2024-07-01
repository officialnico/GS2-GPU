
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the split_nonlinear_terms_knobs namelist and populates the member variables
subroutine read_split_nonlinear_terms_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(split_nonlinear_terms_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: rk_method, split_method
  logical :: advance_nonadiabatic_dfn, show_statistics
  real :: absolute_tolerance, convergence_tolerance, relative_tolerance, time_step_safety_factor

  namelist /split_nonlinear_terms_knobs/ absolute_tolerance, advance_nonadiabatic_dfn, convergence_tolerance, relative_tolerance, rk_method, show_statistics, split_method, time_step_safety_factor

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  absolute_tolerance = self%absolute_tolerance
  advance_nonadiabatic_dfn = self%advance_nonadiabatic_dfn
  convergence_tolerance = self%convergence_tolerance
  relative_tolerance = self%relative_tolerance
  rk_method = self%rk_method
  show_statistics = self%show_statistics
  split_method = self%split_method
  time_step_safety_factor = self%time_step_safety_factor

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = split_nonlinear_terms_knobs)

  ! Now copy from local variables into type members
  self%absolute_tolerance = absolute_tolerance
  self%advance_nonadiabatic_dfn = advance_nonadiabatic_dfn
  self%convergence_tolerance = convergence_tolerance
  self%relative_tolerance = relative_tolerance
  self%rk_method = rk_method
  self%show_statistics = show_statistics
  self%split_method = split_method
  self%time_step_safety_factor = time_step_safety_factor

  self%exist = exist
end subroutine read_split_nonlinear_terms_config

!> Writes out a namelist representing the current state of the config object
subroutine write_split_nonlinear_terms_config(self, unit)
  implicit none
  class(split_nonlinear_terms_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("absolute_tolerance", self%absolute_tolerance, unit_internal)
  call self%write_key_val("advance_nonadiabatic_dfn", self%advance_nonadiabatic_dfn, unit_internal)
  call self%write_key_val("convergence_tolerance", self%convergence_tolerance, unit_internal)
  call self%write_key_val("relative_tolerance", self%relative_tolerance, unit_internal)
  call self%write_key_val("rk_method", self%rk_method, unit_internal)
  call self%write_key_val("show_statistics", self%show_statistics, unit_internal)
  call self%write_key_val("split_method", self%split_method, unit_internal)
  call self%write_key_val("time_step_safety_factor", self%time_step_safety_factor, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_split_nonlinear_terms_config

!> Resets the config object to the initial empty state
subroutine reset_split_nonlinear_terms_config(self)
  class(split_nonlinear_terms_config_type), intent(in out) :: self
  type(split_nonlinear_terms_config_type) :: empty
  select type (self)
  type is (split_nonlinear_terms_config_type)
     self = empty
  end select
end subroutine reset_split_nonlinear_terms_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_split_nonlinear_terms_config(self)
  use mp, only: broadcast
  implicit none
  class(split_nonlinear_terms_config_type), intent(in out) :: self
  call broadcast(self%absolute_tolerance)
  call broadcast(self%advance_nonadiabatic_dfn)
  call broadcast(self%convergence_tolerance)
  call broadcast(self%relative_tolerance)
  call broadcast(self%rk_method)
  call broadcast(self%show_statistics)
  call broadcast(self%split_method)
  call broadcast(self%time_step_safety_factor)

  call broadcast(self%exist)
end subroutine broadcast_split_nonlinear_terms_config

!> Gets the default name for this namelist
function get_default_name_split_nonlinear_terms_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_split_nonlinear_terms_config
  get_default_name_split_nonlinear_terms_config = "split_nonlinear_terms_knobs"
end function get_default_name_split_nonlinear_terms_config

!> Gets the default requires index for this namelist
function get_default_requires_index_split_nonlinear_terms_config()
  implicit none
  logical :: get_default_requires_index_split_nonlinear_terms_config
  get_default_requires_index_split_nonlinear_terms_config = .false.
end function get_default_requires_index_split_nonlinear_terms_config

!> Get the module level config instance
pure function get_split_nonlinear_terms_config()
  type(split_nonlinear_terms_config_type) :: get_split_nonlinear_terms_config
  get_split_nonlinear_terms_config = split_nonlinear_terms_config
end function get_split_nonlinear_terms_config
