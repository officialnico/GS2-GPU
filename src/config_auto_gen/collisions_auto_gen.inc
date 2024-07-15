
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the collisions_knobs namelist and populates the member variables
subroutine read_collisions_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(collisions_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  character(len = 20) :: collision_model, ediff_scheme, lorentz_scheme
  integer :: ncheck, timesteps_between_collisions
  logical :: adjust, conservative, conserve_forbid_zero, conserve_moments, const_v, ei_coll_only, force_collisions, heating, hypermult, resistivity, special_wfb_lorentz, split_collisions, test, use_le_layout
  logical :: vary_vnew, vpar_zero_mean
  real :: cfac, etol, etola, ewindow, ewindowa, vnfac, vnslow

  namelist /collisions_knobs/ adjust, cfac, collision_model, conservative, conserve_forbid_zero, conserve_moments, const_v, ediff_scheme, ei_coll_only, etol, etola, ewindow, ewindowa, force_collisions, heating, hypermult, &
    lorentz_scheme, ncheck, resistivity, special_wfb_lorentz, split_collisions, test, timesteps_between_collisions, use_le_layout, vary_vnew, vnfac, vnslow, vpar_zero_mean

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  adjust = self%adjust
  cfac = self%cfac
  collision_model = self%collision_model
  conservative = self%conservative
  conserve_forbid_zero = self%conserve_forbid_zero
  conserve_moments = self%conserve_moments
  const_v = self%const_v
  ediff_scheme = self%ediff_scheme
  ei_coll_only = self%ei_coll_only
  etol = self%etol
  etola = self%etola
  ewindow = self%ewindow
  ewindowa = self%ewindowa
  force_collisions = self%force_collisions
  heating = self%heating
  hypermult = self%hypermult
  lorentz_scheme = self%lorentz_scheme
  ncheck = self%ncheck
  resistivity = self%resistivity
  special_wfb_lorentz = self%special_wfb_lorentz
  split_collisions = self%split_collisions
  test = self%test
  timesteps_between_collisions = self%timesteps_between_collisions
  use_le_layout = self%use_le_layout
  vary_vnew = self%vary_vnew
  vnfac = self%vnfac
  vnslow = self%vnslow
  vpar_zero_mean = self%vpar_zero_mean

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = collisions_knobs)

  ! Now copy from local variables into type members
  self%adjust = adjust
  self%cfac = cfac
  self%collision_model = collision_model
  self%conservative = conservative
  self%conserve_forbid_zero = conserve_forbid_zero
  self%conserve_moments = conserve_moments
  self%const_v = const_v
  self%ediff_scheme = ediff_scheme
  self%ei_coll_only = ei_coll_only
  self%etol = etol
  self%etola = etola
  self%ewindow = ewindow
  self%ewindowa = ewindowa
  self%force_collisions = force_collisions
  self%heating = heating
  self%hypermult = hypermult
  self%lorentz_scheme = lorentz_scheme
  self%ncheck = ncheck
  self%resistivity = resistivity
  self%special_wfb_lorentz = special_wfb_lorentz
  self%split_collisions = split_collisions
  self%test = test
  self%timesteps_between_collisions = timesteps_between_collisions
  self%use_le_layout = use_le_layout
  self%vary_vnew = vary_vnew
  self%vnfac = vnfac
  self%vnslow = vnslow
  self%vpar_zero_mean = vpar_zero_mean

  self%exist = exist
end subroutine read_collisions_config

!> Writes out a namelist representing the current state of the config object
subroutine write_collisions_config(self, unit)
  implicit none
  class(collisions_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("adjust", self%adjust, unit_internal)
  call self%write_key_val("cfac", self%cfac, unit_internal)
  call self%write_key_val("collision_model", self%collision_model, unit_internal)
  call self%write_key_val("conservative", self%conservative, unit_internal)
  call self%write_key_val("conserve_forbid_zero", self%conserve_forbid_zero, unit_internal)
  call self%write_key_val("conserve_moments", self%conserve_moments, unit_internal)
  call self%write_key_val("const_v", self%const_v, unit_internal)
  call self%write_key_val("ediff_scheme", self%ediff_scheme, unit_internal)
  call self%write_key_val("ei_coll_only", self%ei_coll_only, unit_internal)
  call self%write_key_val("etol", self%etol, unit_internal)
  call self%write_key_val("etola", self%etola, unit_internal)
  call self%write_key_val("ewindow", self%ewindow, unit_internal)
  call self%write_key_val("ewindowa", self%ewindowa, unit_internal)
  call self%write_key_val("force_collisions", self%force_collisions, unit_internal)
  call self%write_key_val("heating", self%heating, unit_internal)
  call self%write_key_val("hypermult", self%hypermult, unit_internal)
  call self%write_key_val("lorentz_scheme", self%lorentz_scheme, unit_internal)
  call self%write_key_val("ncheck", self%ncheck, unit_internal)
  call self%write_key_val("resistivity", self%resistivity, unit_internal)
  call self%write_key_val("special_wfb_lorentz", self%special_wfb_lorentz, unit_internal)
  call self%write_key_val("split_collisions", self%split_collisions, unit_internal)
  call self%write_key_val("test", self%test, unit_internal)
  call self%write_key_val("timesteps_between_collisions", self%timesteps_between_collisions, unit_internal)
  call self%write_key_val("use_le_layout", self%use_le_layout, unit_internal)
  call self%write_key_val("vary_vnew", self%vary_vnew, unit_internal)
  call self%write_key_val("vnfac", self%vnfac, unit_internal)
  call self%write_key_val("vnslow", self%vnslow, unit_internal)
  call self%write_key_val("vpar_zero_mean", self%vpar_zero_mean, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_collisions_config

!> Resets the config object to the initial empty state
subroutine reset_collisions_config(self)
  class(collisions_config_type), intent(in out) :: self
  type(collisions_config_type) :: empty
  select type (self)
  type is (collisions_config_type)
     self = empty
  end select
end subroutine reset_collisions_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_collisions_config(self)
  use mp, only: broadcast
  implicit none
  class(collisions_config_type), intent(in out) :: self
  call broadcast(self%adjust)
  call broadcast(self%cfac)
  call broadcast(self%collision_model)
  call broadcast(self%conservative)
  call broadcast(self%conserve_forbid_zero)
  call broadcast(self%conserve_moments)
  call broadcast(self%const_v)
  call broadcast(self%ediff_scheme)
  call broadcast(self%ei_coll_only)
  call broadcast(self%etol)
  call broadcast(self%etola)
  call broadcast(self%ewindow)
  call broadcast(self%ewindowa)
  call broadcast(self%force_collisions)
  call broadcast(self%heating)
  call broadcast(self%hypermult)
  call broadcast(self%lorentz_scheme)
  call broadcast(self%ncheck)
  call broadcast(self%resistivity)
  call broadcast(self%special_wfb_lorentz)
  call broadcast(self%split_collisions)
  call broadcast(self%test)
  call broadcast(self%timesteps_between_collisions)
  call broadcast(self%use_le_layout)
  call broadcast(self%vary_vnew)
  call broadcast(self%vnfac)
  call broadcast(self%vnslow)
  call broadcast(self%vpar_zero_mean)

  call broadcast(self%exist)
end subroutine broadcast_collisions_config

!> Gets the default name for this namelist
function get_default_name_collisions_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_collisions_config
  get_default_name_collisions_config = "collisions_knobs"
end function get_default_name_collisions_config

!> Gets the default requires index for this namelist
function get_default_requires_index_collisions_config()
  implicit none
  logical :: get_default_requires_index_collisions_config
  get_default_requires_index_collisions_config = .false.
end function get_default_requires_index_collisions_config

!> Get the module level config instance
pure function get_collisions_config()
  type(collisions_config_type) :: get_collisions_config
  get_collisions_config = collisions_config
end function get_collisions_config
