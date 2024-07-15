
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the kt_grids_single_parameters namelist and populates the member variables
subroutine read_kt_grids_single_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(kt_grids_single_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: n0
  real :: akx, aky, rhostar_single, theta0

  namelist /kt_grids_single_parameters/ akx, aky, n0, rhostar_single, theta0

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  akx = self%akx
  aky = self%aky
  n0 = self%n0
  rhostar_single = self%rhostar_single
  theta0 = self%theta0

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = kt_grids_single_parameters)

  ! Now copy from local variables into type members
  self%akx = akx
  self%aky = aky
  self%n0 = n0
  self%rhostar_single = rhostar_single
  self%theta0 = theta0

  self%exist = exist
end subroutine read_kt_grids_single_config

!> Writes out a namelist representing the current state of the config object
subroutine write_kt_grids_single_config(self, unit)
  implicit none
  class(kt_grids_single_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("akx", self%akx, unit_internal)
  call self%write_key_val("aky", self%aky, unit_internal)
  call self%write_key_val("n0", self%n0, unit_internal)
  call self%write_key_val("rhostar_single", self%rhostar_single, unit_internal)
  call self%write_key_val("theta0", self%theta0, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_kt_grids_single_config

!> Resets the config object to the initial empty state
subroutine reset_kt_grids_single_config(self)
  class(kt_grids_single_config_type), intent(in out) :: self
  type(kt_grids_single_config_type) :: empty
  select type (self)
  type is (kt_grids_single_config_type)
     self = empty
  end select
end subroutine reset_kt_grids_single_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_kt_grids_single_config(self)
  use mp, only: broadcast
  implicit none
  class(kt_grids_single_config_type), intent(in out) :: self
  call broadcast(self%akx)
  call broadcast(self%aky)
  call broadcast(self%n0)
  call broadcast(self%rhostar_single)
  call broadcast(self%theta0)

  call broadcast(self%exist)
end subroutine broadcast_kt_grids_single_config

!> Gets the default name for this namelist
function get_default_name_kt_grids_single_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_kt_grids_single_config
  get_default_name_kt_grids_single_config = "kt_grids_single_parameters"
end function get_default_name_kt_grids_single_config

!> Gets the default requires index for this namelist
function get_default_requires_index_kt_grids_single_config()
  implicit none
  logical :: get_default_requires_index_kt_grids_single_config
  get_default_requires_index_kt_grids_single_config = .false.
end function get_default_requires_index_kt_grids_single_config

!> Get the module level config instance
pure function get_kt_grids_single_config()
  type(kt_grids_single_config_type) :: get_kt_grids_single_config
  get_kt_grids_single_config = kt_grids_single_config
end function get_kt_grids_single_config
