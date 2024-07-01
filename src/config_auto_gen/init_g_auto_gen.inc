
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the init_g_knobs namelist and populates the member variables
subroutine read_init_g_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(init_g_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 150) :: restart_dir
  character(len = 20) :: eq_type, ginit_option
  character(len = run_name_size) :: restart_file
  complex, dimension(3) :: aparamp, phiamp, ukxy_pt
  integer :: adj_spec, eq_mode_n, eq_mode_u, ikpar_init, ikx_init, max_mode, proc_to_save_fields, restart_eig_id
  integer, dimension(2) :: ikk, itt
  integer, dimension(3) :: ikkk, ittt
  logical :: chop_side, clean_init, constant_random_flag, even, include_explicit_source_in_restart, initial_condition_is_nonadiabatic_dfn, input_check_recon, left, new_field_init, null_apar, null_bpar, null_phi
  logical :: read_many
  real :: a0, apar0, b0, den0, den1, den2, dphiinit, imfac, init_spec_pow, kpar_init, phifrac, phiinit, phiinit0, phiinit_rand, prof_width, refac, scale, tpar0, tpar1, tpar2, tperp0, tperp1, tperp2, tstart
  real :: upar0, upar1, upar2, width0, zf_init

  namelist /init_g_knobs/ a0, adj_spec, apar0, aparamp, b0, chop_side, clean_init, constant_random_flag, den0, den1, den2, dphiinit, eq_mode_n, eq_mode_u, eq_type, even, ginit_option, ikk, ikkk, ikpar_init, ikx_init, imfac, &
    include_explicit_source_in_restart, init_spec_pow, initial_condition_is_nonadiabatic_dfn, input_check_recon, itt, ittt, kpar_init, left, max_mode, new_field_init, null_apar, null_bpar, null_phi, &
    phiamp, phifrac, phiinit, phiinit0, phiinit_rand, proc_to_save_fields, prof_width, read_many, refac, restart_dir, restart_eig_id, restart_file, scale, tpar0, tpar1, tpar2, tperp0, tperp1, tperp2, &
    tstart, ukxy_pt, upar0, upar1, upar2, width0, zf_init

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  a0 = self%a0
  adj_spec = self%adj_spec
  apar0 = self%apar0
  aparamp = self%aparamp
  b0 = self%b0
  chop_side = self%chop_side
  clean_init = self%clean_init
  constant_random_flag = self%constant_random_flag
  den0 = self%den0
  den1 = self%den1
  den2 = self%den2
  dphiinit = self%dphiinit
  eq_mode_n = self%eq_mode_n
  eq_mode_u = self%eq_mode_u
  eq_type = self%eq_type
  even = self%even
  ginit_option = self%ginit_option
  ikk = self%ikk
  ikkk = self%ikkk
  ikpar_init = self%ikpar_init
  ikx_init = self%ikx_init
  imfac = self%imfac
  include_explicit_source_in_restart = self%include_explicit_source_in_restart
  init_spec_pow = self%init_spec_pow
  initial_condition_is_nonadiabatic_dfn = self%initial_condition_is_nonadiabatic_dfn
  input_check_recon = self%input_check_recon
  itt = self%itt
  ittt = self%ittt
  kpar_init = self%kpar_init
  left = self%left
  max_mode = self%max_mode
  new_field_init = self%new_field_init
  null_apar = self%null_apar
  null_bpar = self%null_bpar
  null_phi = self%null_phi
  phiamp = self%phiamp
  phifrac = self%phifrac
  phiinit = self%phiinit
  phiinit0 = self%phiinit0
  phiinit_rand = self%phiinit_rand
  proc_to_save_fields = self%proc_to_save_fields
  prof_width = self%prof_width
  read_many = self%read_many
  refac = self%refac
  restart_dir = self%restart_dir
  restart_eig_id = self%restart_eig_id
  restart_file = self%restart_file
  scale = self%scale
  tpar0 = self%tpar0
  tpar1 = self%tpar1
  tpar2 = self%tpar2
  tperp0 = self%tperp0
  tperp1 = self%tperp1
  tperp2 = self%tperp2
  tstart = self%tstart
  ukxy_pt = self%ukxy_pt
  upar0 = self%upar0
  upar1 = self%upar1
  upar2 = self%upar2
  width0 = self%width0
  zf_init = self%zf_init

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = init_g_knobs)

  ! Now copy from local variables into type members
  self%a0 = a0
  self%adj_spec = adj_spec
  self%apar0 = apar0
  self%aparamp = aparamp
  self%b0 = b0
  self%chop_side = chop_side
  self%clean_init = clean_init
  self%constant_random_flag = constant_random_flag
  self%den0 = den0
  self%den1 = den1
  self%den2 = den2
  self%dphiinit = dphiinit
  self%eq_mode_n = eq_mode_n
  self%eq_mode_u = eq_mode_u
  self%eq_type = eq_type
  self%even = even
  self%ginit_option = ginit_option
  self%ikk = ikk
  self%ikkk = ikkk
  self%ikpar_init = ikpar_init
  self%ikx_init = ikx_init
  self%imfac = imfac
  self%include_explicit_source_in_restart = include_explicit_source_in_restart
  self%init_spec_pow = init_spec_pow
  self%initial_condition_is_nonadiabatic_dfn = initial_condition_is_nonadiabatic_dfn
  self%input_check_recon = input_check_recon
  self%itt = itt
  self%ittt = ittt
  self%kpar_init = kpar_init
  self%left = left
  self%max_mode = max_mode
  self%new_field_init = new_field_init
  self%null_apar = null_apar
  self%null_bpar = null_bpar
  self%null_phi = null_phi
  self%phiamp = phiamp
  self%phifrac = phifrac
  self%phiinit = phiinit
  self%phiinit0 = phiinit0
  self%phiinit_rand = phiinit_rand
  self%proc_to_save_fields = proc_to_save_fields
  self%prof_width = prof_width
  self%read_many = read_many
  self%refac = refac
  self%restart_dir = restart_dir
  self%restart_eig_id = restart_eig_id
  self%restart_file = restart_file
  self%scale = scale
  self%tpar0 = tpar0
  self%tpar1 = tpar1
  self%tpar2 = tpar2
  self%tperp0 = tperp0
  self%tperp1 = tperp1
  self%tperp2 = tperp2
  self%tstart = tstart
  self%ukxy_pt = ukxy_pt
  self%upar0 = upar0
  self%upar1 = upar1
  self%upar2 = upar2
  self%width0 = width0
  self%zf_init = zf_init

  self%exist = exist
end subroutine read_init_g_config

!> Writes out a namelist representing the current state of the config object
subroutine write_init_g_config(self, unit)
  implicit none
  class(init_g_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("a0", self%a0, unit_internal)
  call self%write_key_val("adj_spec", self%adj_spec, unit_internal)
  call self%write_key_val("apar0", self%apar0, unit_internal)
  call self%write_key_val("aparamp", self%aparamp, unit_internal)
  call self%write_key_val("b0", self%b0, unit_internal)
  call self%write_key_val("chop_side", self%chop_side, unit_internal)
  call self%write_key_val("clean_init", self%clean_init, unit_internal)
  call self%write_key_val("constant_random_flag", self%constant_random_flag, unit_internal)
  call self%write_key_val("den0", self%den0, unit_internal)
  call self%write_key_val("den1", self%den1, unit_internal)
  call self%write_key_val("den2", self%den2, unit_internal)
  call self%write_key_val("dphiinit", self%dphiinit, unit_internal)
  call self%write_key_val("eq_mode_n", self%eq_mode_n, unit_internal)
  call self%write_key_val("eq_mode_u", self%eq_mode_u, unit_internal)
  call self%write_key_val("eq_type", self%eq_type, unit_internal)
  call self%write_key_val("even", self%even, unit_internal)
  call self%write_key_val("ginit_option", self%ginit_option, unit_internal)
  call self%write_key_val("ikk", self%ikk, unit_internal)
  call self%write_key_val("ikkk", self%ikkk, unit_internal)
  call self%write_key_val("ikpar_init", self%ikpar_init, unit_internal)
  call self%write_key_val("ikx_init", self%ikx_init, unit_internal)
  call self%write_key_val("imfac", self%imfac, unit_internal)
  call self%write_key_val("include_explicit_source_in_restart", self%include_explicit_source_in_restart, unit_internal)
  call self%write_key_val("init_spec_pow", self%init_spec_pow, unit_internal)
  call self%write_key_val("initial_condition_is_nonadiabatic_dfn", self%initial_condition_is_nonadiabatic_dfn, unit_internal)
  call self%write_key_val("input_check_recon", self%input_check_recon, unit_internal)
  call self%write_key_val("itt", self%itt, unit_internal)
  call self%write_key_val("ittt", self%ittt, unit_internal)
  call self%write_key_val("kpar_init", self%kpar_init, unit_internal)
  call self%write_key_val("left", self%left, unit_internal)
  call self%write_key_val("max_mode", self%max_mode, unit_internal)
  call self%write_key_val("new_field_init", self%new_field_init, unit_internal)
  call self%write_key_val("null_apar", self%null_apar, unit_internal)
  call self%write_key_val("null_bpar", self%null_bpar, unit_internal)
  call self%write_key_val("null_phi", self%null_phi, unit_internal)
  call self%write_key_val("phiamp", self%phiamp, unit_internal)
  call self%write_key_val("phifrac", self%phifrac, unit_internal)
  call self%write_key_val("phiinit", self%phiinit, unit_internal)
  call self%write_key_val("phiinit0", self%phiinit0, unit_internal)
  call self%write_key_val("phiinit_rand", self%phiinit_rand, unit_internal)
  call self%write_key_val("proc_to_save_fields", self%proc_to_save_fields, unit_internal)
  call self%write_key_val("prof_width", self%prof_width, unit_internal)
  call self%write_key_val("read_many", self%read_many, unit_internal)
  call self%write_key_val("refac", self%refac, unit_internal)
  call self%write_key_val("restart_dir", self%restart_dir, unit_internal)
  call self%write_key_val("restart_eig_id", self%restart_eig_id, unit_internal)
  call self%write_key_val("restart_file", self%restart_file, unit_internal)
  call self%write_key_val("scale", self%scale, unit_internal)
  call self%write_key_val("tpar0", self%tpar0, unit_internal)
  call self%write_key_val("tpar1", self%tpar1, unit_internal)
  call self%write_key_val("tpar2", self%tpar2, unit_internal)
  call self%write_key_val("tperp0", self%tperp0, unit_internal)
  call self%write_key_val("tperp1", self%tperp1, unit_internal)
  call self%write_key_val("tperp2", self%tperp2, unit_internal)
  call self%write_key_val("tstart", self%tstart, unit_internal)
  call self%write_key_val("ukxy_pt", self%ukxy_pt, unit_internal)
  call self%write_key_val("upar0", self%upar0, unit_internal)
  call self%write_key_val("upar1", self%upar1, unit_internal)
  call self%write_key_val("upar2", self%upar2, unit_internal)
  call self%write_key_val("width0", self%width0, unit_internal)
  call self%write_key_val("zf_init", self%zf_init, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_init_g_config

!> Resets the config object to the initial empty state
subroutine reset_init_g_config(self)
  class(init_g_config_type), intent(in out) :: self
  type(init_g_config_type) :: empty
  select type (self)
  type is (init_g_config_type)
     self = empty
  end select
end subroutine reset_init_g_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_init_g_config(self)
  use mp, only: broadcast
  implicit none
  class(init_g_config_type), intent(in out) :: self
  call broadcast(self%a0)
  call broadcast(self%adj_spec)
  call broadcast(self%apar0)
  call broadcast(self%aparamp)
  call broadcast(self%b0)
  call broadcast(self%chop_side)
  call broadcast(self%clean_init)
  call broadcast(self%constant_random_flag)
  call broadcast(self%den0)
  call broadcast(self%den1)
  call broadcast(self%den2)
  call broadcast(self%dphiinit)
  call broadcast(self%eq_mode_n)
  call broadcast(self%eq_mode_u)
  call broadcast(self%eq_type)
  call broadcast(self%even)
  call broadcast(self%ginit_option)
  call broadcast(self%ikk)
  call broadcast(self%ikkk)
  call broadcast(self%ikpar_init)
  call broadcast(self%ikx_init)
  call broadcast(self%imfac)
  call broadcast(self%include_explicit_source_in_restart)
  call broadcast(self%init_spec_pow)
  call broadcast(self%initial_condition_is_nonadiabatic_dfn)
  call broadcast(self%input_check_recon)
  call broadcast(self%itt)
  call broadcast(self%ittt)
  call broadcast(self%kpar_init)
  call broadcast(self%left)
  call broadcast(self%max_mode)
  call broadcast(self%new_field_init)
  call broadcast(self%null_apar)
  call broadcast(self%null_bpar)
  call broadcast(self%null_phi)
  call broadcast(self%phiamp)
  call broadcast(self%phifrac)
  call broadcast(self%phiinit)
  call broadcast(self%phiinit0)
  call broadcast(self%phiinit_rand)
  call broadcast(self%proc_to_save_fields)
  call broadcast(self%prof_width)
  call broadcast(self%read_many)
  call broadcast(self%refac)
  call broadcast(self%restart_dir)
  call broadcast(self%restart_eig_id)
  call broadcast(self%restart_file)
  call broadcast(self%scale)
  call broadcast(self%tpar0)
  call broadcast(self%tpar1)
  call broadcast(self%tpar2)
  call broadcast(self%tperp0)
  call broadcast(self%tperp1)
  call broadcast(self%tperp2)
  call broadcast(self%tstart)
  call broadcast(self%ukxy_pt)
  call broadcast(self%upar0)
  call broadcast(self%upar1)
  call broadcast(self%upar2)
  call broadcast(self%width0)
  call broadcast(self%zf_init)

  call broadcast(self%exist)
end subroutine broadcast_init_g_config

!> Gets the default name for this namelist
function get_default_name_init_g_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_init_g_config
  get_default_name_init_g_config = "init_g_knobs"
end function get_default_name_init_g_config

!> Gets the default requires index for this namelist
function get_default_requires_index_init_g_config()
  implicit none
  logical :: get_default_requires_index_init_g_config
  get_default_requires_index_init_g_config = .false.
end function get_default_requires_index_init_g_config

!> Get the module level config instance
pure function get_init_g_config()
  type(init_g_config_type) :: get_init_g_config
  get_init_g_config = init_g_config
end function get_init_g_config
