
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the kt_grids_range_parameters namelist and populates the member variables
subroutine read_kt_grids_range_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(kt_grids_range_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: kyspacing_option
  integer :: n0_max, n0_min, naky, nn0, ntheta0
  real :: akx_max, akx_min, aky_max, aky_min, rhostar_range, theta0_max, theta0_min

  namelist /kt_grids_range_parameters/ akx_max, akx_min, aky_max, aky_min, kyspacing_option, n0_max, n0_min, naky, nn0, ntheta0, rhostar_range, theta0_max, theta0_min

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  akx_max = self%akx_max
  akx_min = self%akx_min
  aky_max = self%aky_max
  aky_min = self%aky_min
  kyspacing_option = self%kyspacing_option
  n0_max = self%n0_max
  n0_min = self%n0_min
  naky = self%naky
  nn0 = self%nn0
  ntheta0 = self%ntheta0
  rhostar_range = self%rhostar_range
  theta0_max = self%theta0_max
  theta0_min = self%theta0_min

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = kt_grids_range_parameters)

  ! Now copy from local variables into type members
  self%akx_max = akx_max
  self%akx_min = akx_min
  self%aky_max = aky_max
  self%aky_min = aky_min
  self%kyspacing_option = kyspacing_option
  self%n0_max = n0_max
  self%n0_min = n0_min
  self%naky = naky
  self%nn0 = nn0
  self%ntheta0 = ntheta0
  self%rhostar_range = rhostar_range
  self%theta0_max = theta0_max
  self%theta0_min = theta0_min

  self%exist = exist
end subroutine read_kt_grids_range_config

!> Writes out a namelist representing the current state of the config object
subroutine write_kt_grids_range_config(self, unit)
  implicit none
  class(kt_grids_range_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("akx_max", self%akx_max, unit_internal)
  call self%write_key_val("akx_min", self%akx_min, unit_internal)
  call self%write_key_val("aky_max", self%aky_max, unit_internal)
  call self%write_key_val("aky_min", self%aky_min, unit_internal)
  call self%write_key_val("kyspacing_option", self%kyspacing_option, unit_internal)
  call self%write_key_val("n0_max", self%n0_max, unit_internal)
  call self%write_key_val("n0_min", self%n0_min, unit_internal)
  call self%write_key_val("naky", self%naky, unit_internal)
  call self%write_key_val("nn0", self%nn0, unit_internal)
  call self%write_key_val("ntheta0", self%ntheta0, unit_internal)
  call self%write_key_val("rhostar_range", self%rhostar_range, unit_internal)
  call self%write_key_val("theta0_max", self%theta0_max, unit_internal)
  call self%write_key_val("theta0_min", self%theta0_min, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_kt_grids_range_config

!> Resets the config object to the initial empty state
subroutine reset_kt_grids_range_config(self)
  class(kt_grids_range_config_type), intent(in out) :: self
  type(kt_grids_range_config_type) :: empty
  select type (self)
  type is (kt_grids_range_config_type)
     self = empty
  end select
end subroutine reset_kt_grids_range_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_kt_grids_range_config(self)
  use mp, only: broadcast
  implicit none
  class(kt_grids_range_config_type), intent(in out) :: self
  call broadcast(self%akx_max)
  call broadcast(self%akx_min)
  call broadcast(self%aky_max)
  call broadcast(self%aky_min)
  call broadcast(self%kyspacing_option)
  call broadcast(self%n0_max)
  call broadcast(self%n0_min)
  call broadcast(self%naky)
  call broadcast(self%nn0)
  call broadcast(self%ntheta0)
  call broadcast(self%rhostar_range)
  call broadcast(self%theta0_max)
  call broadcast(self%theta0_min)

  call broadcast(self%exist)
end subroutine broadcast_kt_grids_range_config

!> Gets the default name for this namelist
function get_default_name_kt_grids_range_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_range_config
  get_default_name_kt_grids_range_config = "kt_grids_range_parameters"
end function get_default_name_kt_grids_range_config

!> Gets the default requires index for this namelist
function get_default_requires_index_kt_grids_range_config()
  implicit none
  logical :: get_default_requires_index_kt_grids_range_config
  get_default_requires_index_kt_grids_range_config = .false.
end function get_default_requires_index_kt_grids_range_config

!> Get the module level config instance
pure function get_kt_grids_range_config()
  type(kt_grids_range_config_type) :: get_kt_grids_range_config
  get_kt_grids_range_config = kt_grids_range_config
end function get_kt_grids_range_config
