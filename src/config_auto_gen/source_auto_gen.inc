
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the source_knobs namelist and populates the member variables
subroutine read_source_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(source_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: source_option
  real :: gamma0, omega0, phi_ext, source0, t0

  namelist /source_knobs/ gamma0, omega0, phi_ext, source0, source_option, t0

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  gamma0 = self%gamma0
  omega0 = self%omega0
  phi_ext = self%phi_ext
  source0 = self%source0
  source_option = self%source_option
  t0 = self%t0

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = source_knobs)

  ! Now copy from local variables into type members
  self%gamma0 = gamma0
  self%omega0 = omega0
  self%phi_ext = phi_ext
  self%source0 = source0
  self%source_option = source_option
  self%t0 = t0

  self%exist = exist
end subroutine read_source_config

!> Writes out a namelist representing the current state of the config object
subroutine write_source_config(self, unit)
  implicit none
  class(source_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("gamma0", self%gamma0, unit_internal)
  call self%write_key_val("omega0", self%omega0, unit_internal)
  call self%write_key_val("phi_ext", self%phi_ext, unit_internal)
  call self%write_key_val("source0", self%source0, unit_internal)
  call self%write_key_val("source_option", self%source_option, unit_internal)
  call self%write_key_val("t0", self%t0, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_source_config

!> Resets the config object to the initial empty state
subroutine reset_source_config(self)
  class(source_config_type), intent(in out) :: self
  type(source_config_type) :: empty
  select type (self)
  type is (source_config_type)
     self = empty
  end select
end subroutine reset_source_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_source_config(self)
  use mp, only: broadcast
  implicit none
  class(source_config_type), intent(in out) :: self
  call broadcast(self%gamma0)
  call broadcast(self%omega0)
  call broadcast(self%phi_ext)
  call broadcast(self%source0)
  call broadcast(self%source_option)
  call broadcast(self%t0)

  call broadcast(self%exist)
end subroutine broadcast_source_config

!> Gets the default name for this namelist
function get_default_name_source_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_source_config
  get_default_name_source_config = "source_knobs"
end function get_default_name_source_config

!> Gets the default requires index for this namelist
function get_default_requires_index_source_config()
  implicit none
  logical :: get_default_requires_index_source_config
  get_default_requires_index_source_config = .false.
end function get_default_requires_index_source_config

!> Get the module level config instance
pure function get_source_config()
  type(source_config_type) :: get_source_config
  get_source_config = source_config
end function get_source_config
