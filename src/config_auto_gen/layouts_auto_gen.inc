
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the layouts_knobs namelist and populates the member variables
subroutine read_layouts_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(layouts_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 5) :: layout
  character(len = run_name_size) :: fft_wisdom_file
  integer :: nproc_e_lo, nproc_g_lo, nproc_le_lo, nproc_lz_lo, nproc_xxf_lo, nproc_yxf_lo
  logical :: fft_measure_plan, fft_use_wisdom, gf_local_fields, intmom_sub, intspec_sub, local_field_solve, opt_local_copy, opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap, simple_gf_decomposition
  logical :: unbalanced_xxf, unbalanced_yxf
  real :: max_unbalanced_xxf, max_unbalanced_yxf

  namelist /layouts_knobs/ fft_measure_plan, fft_use_wisdom, fft_wisdom_file, gf_local_fields, intmom_sub, intspec_sub, layout, local_field_solve, max_unbalanced_xxf, max_unbalanced_yxf, nproc_e_lo, nproc_g_lo, nproc_le_lo, &
    nproc_lz_lo, nproc_xxf_lo, nproc_yxf_lo, opt_local_copy, opt_redist_nbk, opt_redist_persist, opt_redist_persist_overlap, simple_gf_decomposition, unbalanced_xxf, unbalanced_yxf

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  fft_measure_plan = self%fft_measure_plan
  fft_use_wisdom = self%fft_use_wisdom
  fft_wisdom_file = self%fft_wisdom_file
  gf_local_fields = self%gf_local_fields
  intmom_sub = self%intmom_sub
  intspec_sub = self%intspec_sub
  layout = self%layout
  local_field_solve = self%local_field_solve
  max_unbalanced_xxf = self%max_unbalanced_xxf
  max_unbalanced_yxf = self%max_unbalanced_yxf
  nproc_e_lo = self%nproc_e_lo
  nproc_g_lo = self%nproc_g_lo
  nproc_le_lo = self%nproc_le_lo
  nproc_lz_lo = self%nproc_lz_lo
  nproc_xxf_lo = self%nproc_xxf_lo
  nproc_yxf_lo = self%nproc_yxf_lo
  opt_local_copy = self%opt_local_copy
  opt_redist_nbk = self%opt_redist_nbk
  opt_redist_persist = self%opt_redist_persist
  opt_redist_persist_overlap = self%opt_redist_persist_overlap
  simple_gf_decomposition = self%simple_gf_decomposition
  unbalanced_xxf = self%unbalanced_xxf
  unbalanced_yxf = self%unbalanced_yxf

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = layouts_knobs)

  ! Now copy from local variables into type members
  self%fft_measure_plan = fft_measure_plan
  self%fft_use_wisdom = fft_use_wisdom
  self%fft_wisdom_file = fft_wisdom_file
  self%gf_local_fields = gf_local_fields
  self%intmom_sub = intmom_sub
  self%intspec_sub = intspec_sub
  self%layout = layout
  self%local_field_solve = local_field_solve
  self%max_unbalanced_xxf = max_unbalanced_xxf
  self%max_unbalanced_yxf = max_unbalanced_yxf
  self%nproc_e_lo = nproc_e_lo
  self%nproc_g_lo = nproc_g_lo
  self%nproc_le_lo = nproc_le_lo
  self%nproc_lz_lo = nproc_lz_lo
  self%nproc_xxf_lo = nproc_xxf_lo
  self%nproc_yxf_lo = nproc_yxf_lo
  self%opt_local_copy = opt_local_copy
  self%opt_redist_nbk = opt_redist_nbk
  self%opt_redist_persist = opt_redist_persist
  self%opt_redist_persist_overlap = opt_redist_persist_overlap
  self%simple_gf_decomposition = simple_gf_decomposition
  self%unbalanced_xxf = unbalanced_xxf
  self%unbalanced_yxf = unbalanced_yxf

  self%exist = exist
end subroutine read_layouts_config

!> Writes out a namelist representing the current state of the config object
subroutine write_layouts_config(self, unit)
  implicit none
  class(layouts_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("fft_measure_plan", self%fft_measure_plan, unit_internal)
  call self%write_key_val("fft_use_wisdom", self%fft_use_wisdom, unit_internal)
  call self%write_key_val("fft_wisdom_file", self%fft_wisdom_file, unit_internal)
  call self%write_key_val("gf_local_fields", self%gf_local_fields, unit_internal)
  call self%write_key_val("intmom_sub", self%intmom_sub, unit_internal)
  call self%write_key_val("intspec_sub", self%intspec_sub, unit_internal)
  call self%write_key_val("layout", self%layout, unit_internal)
  call self%write_key_val("local_field_solve", self%local_field_solve, unit_internal)
  call self%write_key_val("max_unbalanced_xxf", self%max_unbalanced_xxf, unit_internal)
  call self%write_key_val("max_unbalanced_yxf", self%max_unbalanced_yxf, unit_internal)
  call self%write_key_val("nproc_e_lo", self%nproc_e_lo, unit_internal)
  call self%write_key_val("nproc_g_lo", self%nproc_g_lo, unit_internal)
  call self%write_key_val("nproc_le_lo", self%nproc_le_lo, unit_internal)
  call self%write_key_val("nproc_lz_lo", self%nproc_lz_lo, unit_internal)
  call self%write_key_val("nproc_xxf_lo", self%nproc_xxf_lo, unit_internal)
  call self%write_key_val("nproc_yxf_lo", self%nproc_yxf_lo, unit_internal)
  call self%write_key_val("opt_local_copy", self%opt_local_copy, unit_internal)
  call self%write_key_val("opt_redist_nbk", self%opt_redist_nbk, unit_internal)
  call self%write_key_val("opt_redist_persist", self%opt_redist_persist, unit_internal)
  call self%write_key_val("opt_redist_persist_overlap", self%opt_redist_persist_overlap, unit_internal)
  call self%write_key_val("simple_gf_decomposition", self%simple_gf_decomposition, unit_internal)
  call self%write_key_val("unbalanced_xxf", self%unbalanced_xxf, unit_internal)
  call self%write_key_val("unbalanced_yxf", self%unbalanced_yxf, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_layouts_config

!> Resets the config object to the initial empty state
subroutine reset_layouts_config(self)
  class(layouts_config_type), intent(in out) :: self
  type(layouts_config_type) :: empty
  select type (self)
  type is (layouts_config_type)
     self = empty
  end select
end subroutine reset_layouts_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_layouts_config(self)
  use mp, only: broadcast
  implicit none
  class(layouts_config_type), intent(in out) :: self
  call broadcast(self%fft_measure_plan)
  call broadcast(self%fft_use_wisdom)
  call broadcast(self%fft_wisdom_file)
  call broadcast(self%gf_local_fields)
  call broadcast(self%intmom_sub)
  call broadcast(self%intspec_sub)
  call broadcast(self%layout)
  call broadcast(self%local_field_solve)
  call broadcast(self%max_unbalanced_xxf)
  call broadcast(self%max_unbalanced_yxf)
  call broadcast(self%nproc_e_lo)
  call broadcast(self%nproc_g_lo)
  call broadcast(self%nproc_le_lo)
  call broadcast(self%nproc_lz_lo)
  call broadcast(self%nproc_xxf_lo)
  call broadcast(self%nproc_yxf_lo)
  call broadcast(self%opt_local_copy)
  call broadcast(self%opt_redist_nbk)
  call broadcast(self%opt_redist_persist)
  call broadcast(self%opt_redist_persist_overlap)
  call broadcast(self%simple_gf_decomposition)
  call broadcast(self%unbalanced_xxf)
  call broadcast(self%unbalanced_yxf)

  call broadcast(self%exist)
end subroutine broadcast_layouts_config

!> Gets the default name for this namelist
function get_default_name_layouts_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_layouts_config
  get_default_name_layouts_config = "layouts_knobs"
end function get_default_name_layouts_config

!> Gets the default requires index for this namelist
function get_default_requires_index_layouts_config()
  implicit none
  logical :: get_default_requires_index_layouts_config
  get_default_requires_index_layouts_config = .false.
end function get_default_requires_index_layouts_config

!> Get the module level config instance
pure function get_layouts_config()
  type(layouts_config_type) :: get_layouts_config
  get_layouts_config = layouts_config
end function get_layouts_config
