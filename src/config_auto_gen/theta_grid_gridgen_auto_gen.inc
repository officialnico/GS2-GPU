
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the theta_grid_gridgen_knobs namelist and populates the member variables
subroutine read_theta_grid_gridgen_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(theta_grid_gridgen_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: npadd
  logical :: skip_gridgen
  real :: alknob, bpknob, deltaw, epsknob, extrknob, regrid_tension, tension, thetamax, widthw

  namelist /theta_grid_gridgen_knobs/ alknob, bpknob, deltaw, epsknob, extrknob, npadd, regrid_tension, skip_gridgen, tension, thetamax, widthw

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  alknob = self%alknob
  bpknob = self%bpknob
  deltaw = self%deltaw
  epsknob = self%epsknob
  extrknob = self%extrknob
  npadd = self%npadd
  regrid_tension = self%regrid_tension
  skip_gridgen = self%skip_gridgen
  tension = self%tension
  thetamax = self%thetamax
  widthw = self%widthw

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = theta_grid_gridgen_knobs)

  ! Now copy from local variables into type members
  self%alknob = alknob
  self%bpknob = bpknob
  self%deltaw = deltaw
  self%epsknob = epsknob
  self%extrknob = extrknob
  self%npadd = npadd
  self%regrid_tension = regrid_tension
  self%skip_gridgen = skip_gridgen
  self%tension = tension
  self%thetamax = thetamax
  self%widthw = widthw

  self%exist = exist
end subroutine read_theta_grid_gridgen_config

!> Writes out a namelist representing the current state of the config object
subroutine write_theta_grid_gridgen_config(self, unit)
  implicit none
  class(theta_grid_gridgen_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("alknob", self%alknob, unit_internal)
  call self%write_key_val("bpknob", self%bpknob, unit_internal)
  call self%write_key_val("deltaw", self%deltaw, unit_internal)
  call self%write_key_val("epsknob", self%epsknob, unit_internal)
  call self%write_key_val("extrknob", self%extrknob, unit_internal)
  call self%write_key_val("npadd", self%npadd, unit_internal)
  call self%write_key_val("regrid_tension", self%regrid_tension, unit_internal)
  call self%write_key_val("skip_gridgen", self%skip_gridgen, unit_internal)
  call self%write_key_val("tension", self%tension, unit_internal)
  call self%write_key_val("thetamax", self%thetamax, unit_internal)
  call self%write_key_val("widthw", self%widthw, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_theta_grid_gridgen_config

!> Resets the config object to the initial empty state
subroutine reset_theta_grid_gridgen_config(self)
  class(theta_grid_gridgen_config_type), intent(in out) :: self
  type(theta_grid_gridgen_config_type) :: empty
  select type (self)
  type is (theta_grid_gridgen_config_type)
     self = empty
  end select
end subroutine reset_theta_grid_gridgen_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_theta_grid_gridgen_config(self)
  use mp, only: broadcast
  implicit none
  class(theta_grid_gridgen_config_type), intent(in out) :: self
  call broadcast(self%alknob)
  call broadcast(self%bpknob)
  call broadcast(self%deltaw)
  call broadcast(self%epsknob)
  call broadcast(self%extrknob)
  call broadcast(self%npadd)
  call broadcast(self%regrid_tension)
  call broadcast(self%skip_gridgen)
  call broadcast(self%tension)
  call broadcast(self%thetamax)
  call broadcast(self%widthw)

  call broadcast(self%exist)
end subroutine broadcast_theta_grid_gridgen_config

!> Gets the default name for this namelist
function get_default_name_theta_grid_gridgen_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_gridgen_config
  get_default_name_theta_grid_gridgen_config = "theta_grid_gridgen_knobs"
end function get_default_name_theta_grid_gridgen_config

!> Gets the default requires index for this namelist
function get_default_requires_index_theta_grid_gridgen_config()
  implicit none
  logical :: get_default_requires_index_theta_grid_gridgen_config
  get_default_requires_index_theta_grid_gridgen_config = .false.
end function get_default_requires_index_theta_grid_gridgen_config

!> Get the module level config instance
pure function get_theta_grid_gridgen_config()
  type(theta_grid_gridgen_config_type) :: get_theta_grid_gridgen_config
  get_theta_grid_gridgen_config = theta_grid_gridgen_config
end function get_theta_grid_gridgen_config
