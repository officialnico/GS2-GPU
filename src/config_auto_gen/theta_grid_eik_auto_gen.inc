
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the theta_grid_eik_knobs namelist and populates the member variables
subroutine read_theta_grid_eik_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(theta_grid_eik_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = eqfile_length) :: eqfile, eqnormfile
  integer :: bishop, iflux, irho, isym, itor, ntheta_geometry
  logical :: chs_eq, dfit_eq, efit_eq, equal_arc, force_sym, gen_eq, gs2d_eq, idfit_eq, local_eq, ppl_eq, transp_eq, use_large_aspect, writelots
  real :: alpha_input, beta_prime_input, delrho, dp_mult, invlp_input, s_hat_input

  namelist /theta_grid_eik_knobs/ alpha_input, beta_prime_input, bishop, chs_eq, delrho, dfit_eq, dp_mult, efit_eq, eqfile, eqnormfile, equal_arc, force_sym, gen_eq, gs2d_eq, idfit_eq, iflux, invlp_input, irho, isym, itor, local_eq, &
    ntheta_geometry, ppl_eq, s_hat_input, transp_eq, use_large_aspect, writelots

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  alpha_input = self%alpha_input
  beta_prime_input = self%beta_prime_input
  bishop = self%bishop
  chs_eq = self%chs_eq
  delrho = self%delrho
  dfit_eq = self%dfit_eq
  dp_mult = self%dp_mult
  efit_eq = self%efit_eq
  eqfile = self%eqfile
  eqnormfile = self%eqnormfile
  equal_arc = self%equal_arc
  force_sym = self%force_sym
  gen_eq = self%gen_eq
  gs2d_eq = self%gs2d_eq
  idfit_eq = self%idfit_eq
  iflux = self%iflux
  invlp_input = self%invlp_input
  irho = self%irho
  isym = self%isym
  itor = self%itor
  local_eq = self%local_eq
  ntheta_geometry = self%ntheta_geometry
  ppl_eq = self%ppl_eq
  s_hat_input = self%s_hat_input
  transp_eq = self%transp_eq
  use_large_aspect = self%use_large_aspect
  writelots = self%writelots

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = theta_grid_eik_knobs)

  ! Now copy from local variables into type members
  self%alpha_input = alpha_input
  self%beta_prime_input = beta_prime_input
  self%bishop = bishop
  self%chs_eq = chs_eq
  self%delrho = delrho
  self%dfit_eq = dfit_eq
  self%dp_mult = dp_mult
  self%efit_eq = efit_eq
  self%eqfile = eqfile
  self%eqnormfile = eqnormfile
  self%equal_arc = equal_arc
  self%force_sym = force_sym
  self%gen_eq = gen_eq
  self%gs2d_eq = gs2d_eq
  self%idfit_eq = idfit_eq
  self%iflux = iflux
  self%invlp_input = invlp_input
  self%irho = irho
  self%isym = isym
  self%itor = itor
  self%local_eq = local_eq
  self%ntheta_geometry = ntheta_geometry
  self%ppl_eq = ppl_eq
  self%s_hat_input = s_hat_input
  self%transp_eq = transp_eq
  self%use_large_aspect = use_large_aspect
  self%writelots = writelots

  self%exist = exist
end subroutine read_theta_grid_eik_config

!> Writes out a namelist representing the current state of the config object
subroutine write_theta_grid_eik_config(self, unit)
  implicit none
  class(theta_grid_eik_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("alpha_input", self%alpha_input, unit_internal)
  call self%write_key_val("beta_prime_input", self%beta_prime_input, unit_internal)
  call self%write_key_val("bishop", self%bishop, unit_internal)
  call self%write_key_val("chs_eq", self%chs_eq, unit_internal)
  call self%write_key_val("delrho", self%delrho, unit_internal)
  call self%write_key_val("dfit_eq", self%dfit_eq, unit_internal)
  call self%write_key_val("dp_mult", self%dp_mult, unit_internal)
  call self%write_key_val("efit_eq", self%efit_eq, unit_internal)
  call self%write_key_val("eqfile", self%eqfile, unit_internal)
  call self%write_key_val("eqnormfile", self%eqnormfile, unit_internal)
  call self%write_key_val("equal_arc", self%equal_arc, unit_internal)
  call self%write_key_val("force_sym", self%force_sym, unit_internal)
  call self%write_key_val("gen_eq", self%gen_eq, unit_internal)
  call self%write_key_val("gs2d_eq", self%gs2d_eq, unit_internal)
  call self%write_key_val("idfit_eq", self%idfit_eq, unit_internal)
  call self%write_key_val("iflux", self%iflux, unit_internal)
  call self%write_key_val("invlp_input", self%invlp_input, unit_internal)
  call self%write_key_val("irho", self%irho, unit_internal)
  call self%write_key_val("isym", self%isym, unit_internal)
  call self%write_key_val("itor", self%itor, unit_internal)
  call self%write_key_val("local_eq", self%local_eq, unit_internal)
  call self%write_key_val("ntheta_geometry", self%ntheta_geometry, unit_internal)
  call self%write_key_val("ppl_eq", self%ppl_eq, unit_internal)
  call self%write_key_val("s_hat_input", self%s_hat_input, unit_internal)
  call self%write_key_val("transp_eq", self%transp_eq, unit_internal)
  call self%write_key_val("use_large_aspect", self%use_large_aspect, unit_internal)
  call self%write_key_val("writelots", self%writelots, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_theta_grid_eik_config

!> Resets the config object to the initial empty state
subroutine reset_theta_grid_eik_config(self)
  class(theta_grid_eik_config_type), intent(in out) :: self
  type(theta_grid_eik_config_type) :: empty
  select type (self)
  type is (theta_grid_eik_config_type)
     self = empty
  end select
end subroutine reset_theta_grid_eik_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_theta_grid_eik_config(self)
  use mp, only: broadcast
  implicit none
  class(theta_grid_eik_config_type), intent(in out) :: self
  call broadcast(self%alpha_input)
  call broadcast(self%beta_prime_input)
  call broadcast(self%bishop)
  call broadcast(self%chs_eq)
  call broadcast(self%delrho)
  call broadcast(self%dfit_eq)
  call broadcast(self%dp_mult)
  call broadcast(self%efit_eq)
  call broadcast(self%eqfile)
  call broadcast(self%eqnormfile)
  call broadcast(self%equal_arc)
  call broadcast(self%force_sym)
  call broadcast(self%gen_eq)
  call broadcast(self%gs2d_eq)
  call broadcast(self%idfit_eq)
  call broadcast(self%iflux)
  call broadcast(self%invlp_input)
  call broadcast(self%irho)
  call broadcast(self%isym)
  call broadcast(self%itor)
  call broadcast(self%local_eq)
  call broadcast(self%ntheta_geometry)
  call broadcast(self%ppl_eq)
  call broadcast(self%s_hat_input)
  call broadcast(self%transp_eq)
  call broadcast(self%use_large_aspect)
  call broadcast(self%writelots)

  call broadcast(self%exist)
end subroutine broadcast_theta_grid_eik_config

!> Gets the default name for this namelist
function get_default_name_theta_grid_eik_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_eik_config
  get_default_name_theta_grid_eik_config = "theta_grid_eik_knobs"
end function get_default_name_theta_grid_eik_config

!> Gets the default requires index for this namelist
function get_default_requires_index_theta_grid_eik_config()
  implicit none
  logical :: get_default_requires_index_theta_grid_eik_config
  get_default_requires_index_theta_grid_eik_config = .false.
end function get_default_requires_index_theta_grid_eik_config

!> Get the module level config instance
pure function get_theta_grid_eik_config()
  type(theta_grid_eik_config_type) :: get_theta_grid_eik_config
  get_theta_grid_eik_config = theta_grid_eik_config
end function get_theta_grid_eik_config
