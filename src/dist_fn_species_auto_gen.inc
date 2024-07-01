
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the dist_fn_species_knobs namelist and populates the member variables
subroutine read_dist_fn_species_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(dist_fn_species_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  real :: bakdif, fexpr

  namelist /dist_fn_species_knobs/ bakdif, fexpr

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  bakdif = self%bakdif
  fexpr = self%fexpr

  ! Now read in the main namelist
  call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index, exist)
  if (exist) read(in_file, nml = dist_fn_species_knobs)
  close(unit = in_file)

  ! Now copy from local variables into type members
  self%bakdif = bakdif
  self%fexpr = fexpr

  self%exist = exist
end subroutine read_dist_fn_species_config

!> Writes out a namelist representing the current state of the config object
subroutine write_dist_fn_species_config(self, unit)
  implicit none
  class(dist_fn_species_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("bakdif", self%bakdif, unit_internal)
  call self%write_key_val("fexpr", self%fexpr, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_dist_fn_species_config

!> Resets the config object to the initial empty state
subroutine reset_dist_fn_species_config(self)
  class(dist_fn_species_config_type), intent(in out) :: self
  type(dist_fn_species_config_type) :: empty
  select type (self)
  type is (dist_fn_species_config_type)
     self = empty
  end select
end subroutine reset_dist_fn_species_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_dist_fn_species_config(self)
  use mp, only: broadcast
  implicit none
  class(dist_fn_species_config_type), intent(in out) :: self
  call broadcast(self%bakdif)
  call broadcast(self%fexpr)

  call broadcast(self%exist)
end subroutine broadcast_dist_fn_species_config

!> Gets the default name for this namelist
function get_default_name_dist_fn_species_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_dist_fn_species_config
  get_default_name_dist_fn_species_config = "dist_fn_species_knobs"
end function get_default_name_dist_fn_species_config

!> Gets the default requires index for this namelist
function get_default_requires_index_dist_fn_species_config()
  implicit none
  logical :: get_default_requires_index_dist_fn_species_config
  get_default_requires_index_dist_fn_species_config = .true.
end function get_default_requires_index_dist_fn_species_config

!> Get the array of module level config instances. If it isn't allocated,
!> then return a zero-length array
pure function get_dist_fn_species_config()
  type(dist_fn_species_config_type), allocatable, dimension(:) :: get_dist_fn_species_config
  if (allocated(dist_fn_species_config)) then
    get_dist_fn_species_config = dist_fn_species_config
  else
    allocate(get_dist_fn_species_config(0))
  end if
end function get_dist_fn_species_config
