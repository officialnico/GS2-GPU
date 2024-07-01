
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the dist_fn_knobs namelist and populates the member variables
subroutine read_dist_fn_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(dist_fn_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: boundary_option, wfbbc_option
  character(len = 30) :: adiabatic_option
  integer :: g_exb_start_timestep
  logical :: boundary_off_grid, def_parity, esv, even, exponential_boundary, gf_lo_integrate, hyper_in_initialisation, lf_decompose, lf_default, mult_imp, nonad_zero, opt_source, start_from_previous_solution
  logical :: zero_forbid
  real :: afilter, apfac, btor_slab, driftknob, exponential_boundary_factor, g_exb, g_exb_error_limit, g_exb_start_time, g_exbfac, gridfac, mach, omprimfac, poisfac, tpdriftknob, vparknob, wfb

  namelist /dist_fn_knobs/ adiabatic_option, afilter, apfac, boundary_off_grid, boundary_option, btor_slab, def_parity, driftknob, esv, even, exponential_boundary, exponential_boundary_factor, g_exb, g_exb_error_limit, &
    g_exb_start_time, g_exb_start_timestep, g_exbfac, gf_lo_integrate, gridfac, hyper_in_initialisation, lf_decompose, lf_default, mach, mult_imp, nonad_zero, omprimfac, opt_source, poisfac, &
    start_from_previous_solution, tpdriftknob, vparknob, wfb, wfbbc_option, zero_forbid

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  adiabatic_option = self%adiabatic_option
  afilter = self%afilter
  apfac = self%apfac
  boundary_off_grid = self%boundary_off_grid
  boundary_option = self%boundary_option
  btor_slab = self%btor_slab
  def_parity = self%def_parity
  driftknob = self%driftknob
  esv = self%esv
  even = self%even
  exponential_boundary = self%exponential_boundary
  exponential_boundary_factor = self%exponential_boundary_factor
  g_exb = self%g_exb
  g_exb_error_limit = self%g_exb_error_limit
  g_exb_start_time = self%g_exb_start_time
  g_exb_start_timestep = self%g_exb_start_timestep
  g_exbfac = self%g_exbfac
  gf_lo_integrate = self%gf_lo_integrate
  gridfac = self%gridfac
  hyper_in_initialisation = self%hyper_in_initialisation
  lf_decompose = self%lf_decompose
  lf_default = self%lf_default
  mach = self%mach
  mult_imp = self%mult_imp
  nonad_zero = self%nonad_zero
  omprimfac = self%omprimfac
  opt_source = self%opt_source
  poisfac = self%poisfac
  start_from_previous_solution = self%start_from_previous_solution
  tpdriftknob = self%tpdriftknob
  vparknob = self%vparknob
  wfb = self%wfb
  wfbbc_option = self%wfbbc_option
  zero_forbid = self%zero_forbid

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = dist_fn_knobs)

  ! Now copy from local variables into type members
  self%adiabatic_option = adiabatic_option
  self%afilter = afilter
  self%apfac = apfac
  self%boundary_off_grid = boundary_off_grid
  self%boundary_option = boundary_option
  self%btor_slab = btor_slab
  self%def_parity = def_parity
  self%driftknob = driftknob
  self%esv = esv
  self%even = even
  self%exponential_boundary = exponential_boundary
  self%exponential_boundary_factor = exponential_boundary_factor
  self%g_exb = g_exb
  self%g_exb_error_limit = g_exb_error_limit
  self%g_exb_start_time = g_exb_start_time
  self%g_exb_start_timestep = g_exb_start_timestep
  self%g_exbfac = g_exbfac
  self%gf_lo_integrate = gf_lo_integrate
  self%gridfac = gridfac
  self%hyper_in_initialisation = hyper_in_initialisation
  self%lf_decompose = lf_decompose
  self%lf_default = lf_default
  self%mach = mach
  self%mult_imp = mult_imp
  self%nonad_zero = nonad_zero
  self%omprimfac = omprimfac
  self%opt_source = opt_source
  self%poisfac = poisfac
  self%start_from_previous_solution = start_from_previous_solution
  self%tpdriftknob = tpdriftknob
  self%vparknob = vparknob
  self%wfb = wfb
  self%wfbbc_option = wfbbc_option
  self%zero_forbid = zero_forbid

  self%exist = exist
end subroutine read_dist_fn_config

!> Writes out a namelist representing the current state of the config object
subroutine write_dist_fn_config(self, unit)
  implicit none
  class(dist_fn_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("adiabatic_option", self%adiabatic_option, unit_internal)
  call self%write_key_val("afilter", self%afilter, unit_internal)
  call self%write_key_val("apfac", self%apfac, unit_internal)
  call self%write_key_val("boundary_off_grid", self%boundary_off_grid, unit_internal)
  call self%write_key_val("boundary_option", self%boundary_option, unit_internal)
  call self%write_key_val("btor_slab", self%btor_slab, unit_internal)
  call self%write_key_val("def_parity", self%def_parity, unit_internal)
  call self%write_key_val("driftknob", self%driftknob, unit_internal)
  call self%write_key_val("esv", self%esv, unit_internal)
  call self%write_key_val("even", self%even, unit_internal)
  call self%write_key_val("exponential_boundary", self%exponential_boundary, unit_internal)
  call self%write_key_val("exponential_boundary_factor", self%exponential_boundary_factor, unit_internal)
  call self%write_key_val("g_exb", self%g_exb, unit_internal)
  call self%write_key_val("g_exb_error_limit", self%g_exb_error_limit, unit_internal)
  call self%write_key_val("g_exb_start_time", self%g_exb_start_time, unit_internal)
  call self%write_key_val("g_exb_start_timestep", self%g_exb_start_timestep, unit_internal)
  call self%write_key_val("g_exbfac", self%g_exbfac, unit_internal)
  call self%write_key_val("gf_lo_integrate", self%gf_lo_integrate, unit_internal)
  call self%write_key_val("gridfac", self%gridfac, unit_internal)
  call self%write_key_val("hyper_in_initialisation", self%hyper_in_initialisation, unit_internal)
  call self%write_key_val("lf_decompose", self%lf_decompose, unit_internal)
  call self%write_key_val("lf_default", self%lf_default, unit_internal)
  call self%write_key_val("mach", self%mach, unit_internal)
  call self%write_key_val("mult_imp", self%mult_imp, unit_internal)
  call self%write_key_val("nonad_zero", self%nonad_zero, unit_internal)
  call self%write_key_val("omprimfac", self%omprimfac, unit_internal)
  call self%write_key_val("opt_source", self%opt_source, unit_internal)
  call self%write_key_val("poisfac", self%poisfac, unit_internal)
  call self%write_key_val("start_from_previous_solution", self%start_from_previous_solution, unit_internal)
  call self%write_key_val("tpdriftknob", self%tpdriftknob, unit_internal)
  call self%write_key_val("vparknob", self%vparknob, unit_internal)
  call self%write_key_val("wfb", self%wfb, unit_internal)
  call self%write_key_val("wfbbc_option", self%wfbbc_option, unit_internal)
  call self%write_key_val("zero_forbid", self%zero_forbid, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_dist_fn_config

!> Resets the config object to the initial empty state
subroutine reset_dist_fn_config(self)
  class(dist_fn_config_type), intent(in out) :: self
  type(dist_fn_config_type) :: empty
  select type (self)
  type is (dist_fn_config_type)
     self = empty
  end select
end subroutine reset_dist_fn_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_dist_fn_config(self)
  use mp, only: broadcast
  implicit none
  class(dist_fn_config_type), intent(in out) :: self
  call broadcast(self%adiabatic_option)
  call broadcast(self%afilter)
  call broadcast(self%apfac)
  call broadcast(self%boundary_off_grid)
  call broadcast(self%boundary_option)
  call broadcast(self%btor_slab)
  call broadcast(self%def_parity)
  call broadcast(self%driftknob)
  call broadcast(self%esv)
  call broadcast(self%even)
  call broadcast(self%exponential_boundary)
  call broadcast(self%exponential_boundary_factor)
  call broadcast(self%g_exb)
  call broadcast(self%g_exb_error_limit)
  call broadcast(self%g_exb_start_time)
  call broadcast(self%g_exb_start_timestep)
  call broadcast(self%g_exbfac)
  call broadcast(self%gf_lo_integrate)
  call broadcast(self%gridfac)
  call broadcast(self%hyper_in_initialisation)
  call broadcast(self%lf_decompose)
  call broadcast(self%lf_default)
  call broadcast(self%mach)
  call broadcast(self%mult_imp)
  call broadcast(self%nonad_zero)
  call broadcast(self%omprimfac)
  call broadcast(self%opt_source)
  call broadcast(self%poisfac)
  call broadcast(self%start_from_previous_solution)
  call broadcast(self%tpdriftknob)
  call broadcast(self%vparknob)
  call broadcast(self%wfb)
  call broadcast(self%wfbbc_option)
  call broadcast(self%zero_forbid)

  call broadcast(self%exist)
end subroutine broadcast_dist_fn_config

!> Gets the default name for this namelist
function get_default_name_dist_fn_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_dist_fn_config
  get_default_name_dist_fn_config = "dist_fn_knobs"
end function get_default_name_dist_fn_config

!> Gets the default requires index for this namelist
function get_default_requires_index_dist_fn_config()
  implicit none
  logical :: get_default_requires_index_dist_fn_config
  get_default_requires_index_dist_fn_config = .false.
end function get_default_requires_index_dist_fn_config

!> Get the module level config instance
pure function get_dist_fn_config()
  type(dist_fn_config_type) :: get_dist_fn_config
  get_dist_fn_config = dist_fn_config
end function get_dist_fn_config
