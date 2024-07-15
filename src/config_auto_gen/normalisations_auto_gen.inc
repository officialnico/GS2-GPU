
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the normalisations_knobs namelist and populates the member variables
subroutine read_normalisations_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(normalisations_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  real :: aref, bref, mref, nref, rhoref, tref, vref, zref

  namelist /normalisations_knobs/ aref, bref, mref, nref, rhoref, tref, vref, zref

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  aref = self%aref
  bref = self%bref
  mref = self%mref
  nref = self%nref
  rhoref = self%rhoref
  tref = self%tref
  vref = self%vref
  zref = self%zref

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = normalisations_knobs)

  ! Now copy from local variables into type members
  self%aref = aref
  self%bref = bref
  self%mref = mref
  self%nref = nref
  self%rhoref = rhoref
  self%tref = tref
  self%vref = vref
  self%zref = zref

  self%exist = exist
end subroutine read_normalisations_config

!> Writes out a namelist representing the current state of the config object
subroutine write_normalisations_config(self, unit)
  implicit none
  class(normalisations_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("aref", self%aref, unit_internal)
  call self%write_key_val("bref", self%bref, unit_internal)
  call self%write_key_val("mref", self%mref, unit_internal)
  call self%write_key_val("nref", self%nref, unit_internal)
  call self%write_key_val("rhoref", self%rhoref, unit_internal)
  call self%write_key_val("tref", self%tref, unit_internal)
  call self%write_key_val("vref", self%vref, unit_internal)
  call self%write_key_val("zref", self%zref, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_normalisations_config

!> Resets the config object to the initial empty state
subroutine reset_normalisations_config(self)
  class(normalisations_config_type), intent(in out) :: self
  type(normalisations_config_type) :: empty
  select type (self)
  type is (normalisations_config_type)
     self = empty
  end select
end subroutine reset_normalisations_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_normalisations_config(self)
  use mp, only: broadcast
  implicit none
  class(normalisations_config_type), intent(in out) :: self
  call broadcast(self%aref)
  call broadcast(self%bref)
  call broadcast(self%mref)
  call broadcast(self%nref)
  call broadcast(self%rhoref)
  call broadcast(self%tref)
  call broadcast(self%vref)
  call broadcast(self%zref)

  call broadcast(self%exist)
end subroutine broadcast_normalisations_config

!> Gets the default name for this namelist
function get_default_name_normalisations_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_normalisations_config
  get_default_name_normalisations_config = "normalisations_knobs"
end function get_default_name_normalisations_config

!> Gets the default requires index for this namelist
function get_default_requires_index_normalisations_config()
  implicit none
  logical :: get_default_requires_index_normalisations_config
  get_default_requires_index_normalisations_config = .false.
end function get_default_requires_index_normalisations_config

!> Get the module level config instance
pure function get_normalisations_config()
  type(normalisations_config_type) :: get_normalisations_config
  get_normalisations_config = normalisations_config
end function get_normalisations_config
