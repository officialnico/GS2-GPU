
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the le_grids_knobs namelist and populates the member variables
subroutine read_le_grids_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(le_grids_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: wfbbc_option
  integer :: negrid, nesub, nesuper, ngauss, nmax, npassing
  logical :: genquad, new_trap_int, new_trap_int_split, radau_gauss_grid, split_passing_region, test, trapped_particles
  real :: bouncefuzz, vcut, wgt_fac

  namelist /le_grids_knobs/ bouncefuzz, genquad, negrid, nesub, nesuper, new_trap_int, new_trap_int_split, ngauss, nmax, npassing, radau_gauss_grid, split_passing_region, test, trapped_particles, vcut, wfbbc_option, wgt_fac

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  bouncefuzz = self%bouncefuzz
  genquad = self%genquad
  negrid = self%negrid
  nesub = self%nesub
  nesuper = self%nesuper
  new_trap_int = self%new_trap_int
  new_trap_int_split = self%new_trap_int_split
  ngauss = self%ngauss
  nmax = self%nmax
  npassing = self%npassing
  radau_gauss_grid = self%radau_gauss_grid
  split_passing_region = self%split_passing_region
  test = self%test
  trapped_particles = self%trapped_particles
  vcut = self%vcut
  wfbbc_option = self%wfbbc_option
  wgt_fac = self%wgt_fac

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = le_grids_knobs)

  ! Now copy from local variables into type members
  self%bouncefuzz = bouncefuzz
  self%genquad = genquad
  self%negrid = negrid
  self%nesub = nesub
  self%nesuper = nesuper
  self%new_trap_int = new_trap_int
  self%new_trap_int_split = new_trap_int_split
  self%ngauss = ngauss
  self%nmax = nmax
  self%npassing = npassing
  self%radau_gauss_grid = radau_gauss_grid
  self%split_passing_region = split_passing_region
  self%test = test
  self%trapped_particles = trapped_particles
  self%vcut = vcut
  self%wfbbc_option = wfbbc_option
  self%wgt_fac = wgt_fac

  self%exist = exist
end subroutine read_le_grids_config

!> Writes out a namelist representing the current state of the config object
subroutine write_le_grids_config(self, unit)
  implicit none
  class(le_grids_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("bouncefuzz", self%bouncefuzz, unit_internal)
  call self%write_key_val("genquad", self%genquad, unit_internal)
  call self%write_key_val("negrid", self%negrid, unit_internal)
  call self%write_key_val("nesub", self%nesub, unit_internal)
  call self%write_key_val("nesuper", self%nesuper, unit_internal)
  call self%write_key_val("new_trap_int", self%new_trap_int, unit_internal)
  call self%write_key_val("new_trap_int_split", self%new_trap_int_split, unit_internal)
  call self%write_key_val("ngauss", self%ngauss, unit_internal)
  call self%write_key_val("nmax", self%nmax, unit_internal)
  call self%write_key_val("npassing", self%npassing, unit_internal)
  call self%write_key_val("radau_gauss_grid", self%radau_gauss_grid, unit_internal)
  call self%write_key_val("split_passing_region", self%split_passing_region, unit_internal)
  call self%write_key_val("test", self%test, unit_internal)
  call self%write_key_val("trapped_particles", self%trapped_particles, unit_internal)
  call self%write_key_val("vcut", self%vcut, unit_internal)
  call self%write_key_val("wfbbc_option", self%wfbbc_option, unit_internal)
  call self%write_key_val("wgt_fac", self%wgt_fac, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_le_grids_config

!> Resets the config object to the initial empty state
subroutine reset_le_grids_config(self)
  class(le_grids_config_type), intent(in out) :: self
  type(le_grids_config_type) :: empty
  select type (self)
  type is (le_grids_config_type)
     self = empty
  end select
end subroutine reset_le_grids_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_le_grids_config(self)
  use mp, only: broadcast
  implicit none
  class(le_grids_config_type), intent(in out) :: self
  call broadcast(self%bouncefuzz)
  call broadcast(self%genquad)
  call broadcast(self%negrid)
  call broadcast(self%nesub)
  call broadcast(self%nesuper)
  call broadcast(self%new_trap_int)
  call broadcast(self%new_trap_int_split)
  call broadcast(self%ngauss)
  call broadcast(self%nmax)
  call broadcast(self%npassing)
  call broadcast(self%radau_gauss_grid)
  call broadcast(self%split_passing_region)
  call broadcast(self%test)
  call broadcast(self%trapped_particles)
  call broadcast(self%vcut)
  call broadcast(self%wfbbc_option)
  call broadcast(self%wgt_fac)

  call broadcast(self%exist)
end subroutine broadcast_le_grids_config

!> Gets the default name for this namelist
function get_default_name_le_grids_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_le_grids_config
  get_default_name_le_grids_config = "le_grids_knobs"
end function get_default_name_le_grids_config

!> Gets the default requires index for this namelist
function get_default_requires_index_le_grids_config()
  implicit none
  logical :: get_default_requires_index_le_grids_config
  get_default_requires_index_le_grids_config = .false.
end function get_default_requires_index_le_grids_config

!> Get the module level config instance
pure function get_le_grids_config()
  type(le_grids_config_type) :: get_le_grids_config
  get_le_grids_config = le_grids_config
end function get_le_grids_config
