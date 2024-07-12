
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the hyper_knobs namelist and populates the member variables
subroutine read_hyper_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(hyper_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 9) :: hyper_option
  integer :: nexp
  logical :: const_amp, damp_zonal_only, gridnorm, include_kpar, isotropic_model, isotropic_shear
  real :: d_hyper, d_hyper3d, d_hyperres, d_hypervisc, kx_cut, ky_cut, omega_osc, p_hyper3d

  namelist /hyper_knobs/ const_amp, d_hyper, d_hyper3d, d_hyperres, d_hypervisc, damp_zonal_only, gridnorm, hyper_option, include_kpar, isotropic_model, isotropic_shear, kx_cut, ky_cut, nexp, omega_osc, p_hyper3d

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  const_amp = self%const_amp
  d_hyper = self%d_hyper
  d_hyper3d = self%d_hyper3d
  d_hyperres = self%d_hyperres
  d_hypervisc = self%d_hypervisc
  damp_zonal_only = self%damp_zonal_only
  gridnorm = self%gridnorm
  hyper_option = self%hyper_option
  include_kpar = self%include_kpar
  isotropic_model = self%isotropic_model
  isotropic_shear = self%isotropic_shear
  kx_cut = self%kx_cut
  ky_cut = self%ky_cut
  nexp = self%nexp
  omega_osc = self%omega_osc
  p_hyper3d = self%p_hyper3d

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = hyper_knobs)

  ! Now copy from local variables into type members
  self%const_amp = const_amp
  self%d_hyper = d_hyper
  self%d_hyper3d = d_hyper3d
  self%d_hyperres = d_hyperres
  self%d_hypervisc = d_hypervisc
  self%damp_zonal_only = damp_zonal_only
  self%gridnorm = gridnorm
  self%hyper_option = hyper_option
  self%include_kpar = include_kpar
  self%isotropic_model = isotropic_model
  self%isotropic_shear = isotropic_shear
  self%kx_cut = kx_cut
  self%ky_cut = ky_cut
  self%nexp = nexp
  self%omega_osc = omega_osc
  self%p_hyper3d = p_hyper3d

  self%exist = exist
end subroutine read_hyper_config

!> Writes out a namelist representing the current state of the config object
subroutine write_hyper_config(self, unit)
  implicit none
  class(hyper_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("const_amp", self%const_amp, unit_internal)
  call self%write_key_val("d_hyper", self%d_hyper, unit_internal)
  call self%write_key_val("d_hyper3d", self%d_hyper3d, unit_internal)
  call self%write_key_val("d_hyperres", self%d_hyperres, unit_internal)
  call self%write_key_val("d_hypervisc", self%d_hypervisc, unit_internal)
  call self%write_key_val("damp_zonal_only", self%damp_zonal_only, unit_internal)
  call self%write_key_val("gridnorm", self%gridnorm, unit_internal)
  call self%write_key_val("hyper_option", self%hyper_option, unit_internal)
  call self%write_key_val("include_kpar", self%include_kpar, unit_internal)
  call self%write_key_val("isotropic_model", self%isotropic_model, unit_internal)
  call self%write_key_val("isotropic_shear", self%isotropic_shear, unit_internal)
  call self%write_key_val("kx_cut", self%kx_cut, unit_internal)
  call self%write_key_val("ky_cut", self%ky_cut, unit_internal)
  call self%write_key_val("nexp", self%nexp, unit_internal)
  call self%write_key_val("omega_osc", self%omega_osc, unit_internal)
  call self%write_key_val("p_hyper3d", self%p_hyper3d, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_hyper_config

!> Resets the config object to the initial empty state
subroutine reset_hyper_config(self)
  class(hyper_config_type), intent(in out) :: self
  type(hyper_config_type) :: empty
  select type (self)
  type is (hyper_config_type)
     self = empty
  end select
end subroutine reset_hyper_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_hyper_config(self)
  use mp, only: broadcast
  implicit none
  class(hyper_config_type), intent(in out) :: self
  call broadcast(self%const_amp)
  call broadcast(self%d_hyper)
  call broadcast(self%d_hyper3d)
  call broadcast(self%d_hyperres)
  call broadcast(self%d_hypervisc)
  call broadcast(self%damp_zonal_only)
  call broadcast(self%gridnorm)
  call broadcast(self%hyper_option)
  call broadcast(self%include_kpar)
  call broadcast(self%isotropic_model)
  call broadcast(self%isotropic_shear)
  call broadcast(self%kx_cut)
  call broadcast(self%ky_cut)
  call broadcast(self%nexp)
  call broadcast(self%omega_osc)
  call broadcast(self%p_hyper3d)

  call broadcast(self%exist)
end subroutine broadcast_hyper_config

!> Gets the default name for this namelist
function get_default_name_hyper_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_hyper_config
  get_default_name_hyper_config = "hyper_knobs"
end function get_default_name_hyper_config

!> Gets the default requires index for this namelist
function get_default_requires_index_hyper_config()
  implicit none
  logical :: get_default_requires_index_hyper_config
  get_default_requires_index_hyper_config = .false.
end function get_default_requires_index_hyper_config

!> Get the module level config instance
pure function get_hyper_config()
  type(hyper_config_type) :: get_hyper_config
  get_hyper_config = hyper_config
end function get_hyper_config
