
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the nonlinear_terms_knobs namelist and populates the member variables
subroutine read_nonlinear_terms_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(nonlinear_terms_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: nonlinear_mode
  integer :: istep_error_start
  logical :: include_apar, include_bpar, include_phi, nl_forbid_force_zero, split_nonlinear, use_2d_cfl, use_cfl_limit, use_order_based_error, zip
  real :: cfl, error_target

  namelist /nonlinear_terms_knobs/ cfl, error_target, include_apar, include_bpar, include_phi, istep_error_start, nl_forbid_force_zero, nonlinear_mode, split_nonlinear, use_2d_cfl, use_cfl_limit, use_order_based_error, zip

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  cfl = self%cfl
  error_target = self%error_target
  include_apar = self%include_apar
  include_bpar = self%include_bpar
  include_phi = self%include_phi
  istep_error_start = self%istep_error_start
  nl_forbid_force_zero = self%nl_forbid_force_zero
  nonlinear_mode = self%nonlinear_mode
  split_nonlinear = self%split_nonlinear
  use_2d_cfl = self%use_2d_cfl
  use_cfl_limit = self%use_cfl_limit
  use_order_based_error = self%use_order_based_error
  zip = self%zip

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = nonlinear_terms_knobs)

  ! Now copy from local variables into type members
  self%cfl = cfl
  self%error_target = error_target
  self%include_apar = include_apar
  self%include_bpar = include_bpar
  self%include_phi = include_phi
  self%istep_error_start = istep_error_start
  self%nl_forbid_force_zero = nl_forbid_force_zero
  self%nonlinear_mode = nonlinear_mode
  self%split_nonlinear = split_nonlinear
  self%use_2d_cfl = use_2d_cfl
  self%use_cfl_limit = use_cfl_limit
  self%use_order_based_error = use_order_based_error
  self%zip = zip

  self%exist = exist
end subroutine read_nonlinear_terms_config

!> Writes out a namelist representing the current state of the config object
subroutine write_nonlinear_terms_config(self, unit)
  implicit none
  class(nonlinear_terms_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("cfl", self%cfl, unit_internal)
  call self%write_key_val("error_target", self%error_target, unit_internal)
  call self%write_key_val("include_apar", self%include_apar, unit_internal)
  call self%write_key_val("include_bpar", self%include_bpar, unit_internal)
  call self%write_key_val("include_phi", self%include_phi, unit_internal)
  call self%write_key_val("istep_error_start", self%istep_error_start, unit_internal)
  call self%write_key_val("nl_forbid_force_zero", self%nl_forbid_force_zero, unit_internal)
  call self%write_key_val("nonlinear_mode", self%nonlinear_mode, unit_internal)
  call self%write_key_val("split_nonlinear", self%split_nonlinear, unit_internal)
  call self%write_key_val("use_2d_cfl", self%use_2d_cfl, unit_internal)
  call self%write_key_val("use_cfl_limit", self%use_cfl_limit, unit_internal)
  call self%write_key_val("use_order_based_error", self%use_order_based_error, unit_internal)
  call self%write_key_val("zip", self%zip, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_nonlinear_terms_config

!> Resets the config object to the initial empty state
subroutine reset_nonlinear_terms_config(self)
  class(nonlinear_terms_config_type), intent(in out) :: self
  type(nonlinear_terms_config_type) :: empty
  select type (self)
  type is (nonlinear_terms_config_type)
     self = empty
  end select
end subroutine reset_nonlinear_terms_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_nonlinear_terms_config(self)
  use mp, only: broadcast
  implicit none
  class(nonlinear_terms_config_type), intent(in out) :: self
  call broadcast(self%cfl)
  call broadcast(self%error_target)
  call broadcast(self%include_apar)
  call broadcast(self%include_bpar)
  call broadcast(self%include_phi)
  call broadcast(self%istep_error_start)
  call broadcast(self%nl_forbid_force_zero)
  call broadcast(self%nonlinear_mode)
  call broadcast(self%split_nonlinear)
  call broadcast(self%use_2d_cfl)
  call broadcast(self%use_cfl_limit)
  call broadcast(self%use_order_based_error)
  call broadcast(self%zip)

  call broadcast(self%exist)
end subroutine broadcast_nonlinear_terms_config

!> Gets the default name for this namelist
function get_default_name_nonlinear_terms_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_nonlinear_terms_config
  get_default_name_nonlinear_terms_config = "nonlinear_terms_knobs"
end function get_default_name_nonlinear_terms_config

!> Gets the default requires index for this namelist
function get_default_requires_index_nonlinear_terms_config()
  implicit none
  logical :: get_default_requires_index_nonlinear_terms_config
  get_default_requires_index_nonlinear_terms_config = .false.
end function get_default_requires_index_nonlinear_terms_config

!> Get the module level config instance
pure function get_nonlinear_terms_config()
  type(nonlinear_terms_config_type) :: get_nonlinear_terms_config
  get_nonlinear_terms_config = nonlinear_terms_config
end function get_nonlinear_terms_config
