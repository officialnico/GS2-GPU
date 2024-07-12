
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the theta_grid_parameters namelist and populates the member variables
subroutine read_theta_grid_parameters_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(theta_grid_parameters_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: geotype, mmode, n_mxh, nmode, nperiod, ntheta
  real :: akappa, akappri, alpmhd, asurf, btor_slab, c0_mxh, delta2, delta3, deltam, deltampri, deltan, deltanpri, eps, epsl, kp, pk, qinp, r_geo, raxis, rhoc, rmaj, shat, shift, shiftvert, theta2, theta3
  real :: thetad, thetak, thetam, thetan, tri, tripri, zaxis
  real, dimension(mxh_max_moments) :: c_mxh, dc_mxh_dr, ds_mxh_dr, s_mxh

  namelist /theta_grid_parameters/ akappa, akappri, alpmhd, asurf, btor_slab, c0_mxh, c_mxh, dc_mxh_dr, delta2, delta3, deltam, deltampri, deltan, deltanpri, ds_mxh_dr, eps, epsl, geotype, kp, mmode, n_mxh, nmode, nperiod, ntheta, pk, &
    qinp, r_geo, raxis, rhoc, rmaj, s_mxh, shat, shift, shiftvert, theta2, theta3, thetad, thetak, thetam, thetan, tri, tripri, zaxis

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  akappa = self%akappa
  akappri = self%akappri
  alpmhd = self%alpmhd
  asurf = self%asurf
  btor_slab = self%btor_slab
  c0_mxh = self%c0_mxh
  c_mxh = self%c_mxh
  dc_mxh_dr = self%dc_mxh_dr
  delta2 = self%delta2
  delta3 = self%delta3
  deltam = self%deltam
  deltampri = self%deltampri
  deltan = self%deltan
  deltanpri = self%deltanpri
  ds_mxh_dr = self%ds_mxh_dr
  eps = self%eps
  epsl = self%epsl
  geotype = self%geotype
  kp = self%kp
  mmode = self%mmode
  n_mxh = self%n_mxh
  nmode = self%nmode
  nperiod = self%nperiod
  ntheta = self%ntheta
  pk = self%pk
  qinp = self%qinp
  r_geo = self%r_geo
  raxis = self%raxis
  rhoc = self%rhoc
  rmaj = self%rmaj
  s_mxh = self%s_mxh
  shat = self%shat
  shift = self%shift
  shiftvert = self%shiftvert
  theta2 = self%theta2
  theta3 = self%theta3
  thetad = self%thetad
  thetak = self%thetak
  thetam = self%thetam
  thetan = self%thetan
  tri = self%tri
  tripri = self%tripri
  zaxis = self%zaxis

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = theta_grid_parameters)

  ! Now copy from local variables into type members
  self%akappa = akappa
  self%akappri = akappri
  self%alpmhd = alpmhd
  self%asurf = asurf
  self%btor_slab = btor_slab
  self%c0_mxh = c0_mxh
  self%c_mxh = c_mxh
  self%dc_mxh_dr = dc_mxh_dr
  self%delta2 = delta2
  self%delta3 = delta3
  self%deltam = deltam
  self%deltampri = deltampri
  self%deltan = deltan
  self%deltanpri = deltanpri
  self%ds_mxh_dr = ds_mxh_dr
  self%eps = eps
  self%epsl = epsl
  self%geotype = geotype
  self%kp = kp
  self%mmode = mmode
  self%n_mxh = n_mxh
  self%nmode = nmode
  self%nperiod = nperiod
  self%ntheta = ntheta
  self%pk = pk
  self%qinp = qinp
  self%r_geo = r_geo
  self%raxis = raxis
  self%rhoc = rhoc
  self%rmaj = rmaj
  self%s_mxh = s_mxh
  self%shat = shat
  self%shift = shift
  self%shiftvert = shiftvert
  self%theta2 = theta2
  self%theta3 = theta3
  self%thetad = thetad
  self%thetak = thetak
  self%thetam = thetam
  self%thetan = thetan
  self%tri = tri
  self%tripri = tripri
  self%zaxis = zaxis

  self%exist = exist
end subroutine read_theta_grid_parameters_config

!> Writes out a namelist representing the current state of the config object
subroutine write_theta_grid_parameters_config(self, unit)
  implicit none
  class(theta_grid_parameters_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("akappa", self%akappa, unit_internal)
  call self%write_key_val("akappri", self%akappri, unit_internal)
  call self%write_key_val("alpmhd", self%alpmhd, unit_internal)
  call self%write_key_val("asurf", self%asurf, unit_internal)
  call self%write_key_val("btor_slab", self%btor_slab, unit_internal)
  call self%write_key_val("c0_mxh", self%c0_mxh, unit_internal)
  call self%write_key_val("c_mxh", self%c_mxh, unit_internal)
  call self%write_key_val("dc_mxh_dr", self%dc_mxh_dr, unit_internal)
  call self%write_key_val("delta2", self%delta2, unit_internal)
  call self%write_key_val("delta3", self%delta3, unit_internal)
  call self%write_key_val("deltam", self%deltam, unit_internal)
  call self%write_key_val("deltampri", self%deltampri, unit_internal)
  call self%write_key_val("deltan", self%deltan, unit_internal)
  call self%write_key_val("deltanpri", self%deltanpri, unit_internal)
  call self%write_key_val("ds_mxh_dr", self%ds_mxh_dr, unit_internal)
  call self%write_key_val("eps", self%eps, unit_internal)
  call self%write_key_val("epsl", self%epsl, unit_internal)
  call self%write_key_val("geotype", self%geotype, unit_internal)
  call self%write_key_val("kp", self%kp, unit_internal)
  call self%write_key_val("mmode", self%mmode, unit_internal)
  call self%write_key_val("n_mxh", self%n_mxh, unit_internal)
  call self%write_key_val("nmode", self%nmode, unit_internal)
  call self%write_key_val("nperiod", self%nperiod, unit_internal)
  call self%write_key_val("ntheta", self%ntheta, unit_internal)
  call self%write_key_val("pk", self%pk, unit_internal)
  call self%write_key_val("qinp", self%qinp, unit_internal)
  call self%write_key_val("r_geo", self%r_geo, unit_internal)
  call self%write_key_val("raxis", self%raxis, unit_internal)
  call self%write_key_val("rhoc", self%rhoc, unit_internal)
  call self%write_key_val("rmaj", self%rmaj, unit_internal)
  call self%write_key_val("s_mxh", self%s_mxh, unit_internal)
  call self%write_key_val("shat", self%shat, unit_internal)
  call self%write_key_val("shift", self%shift, unit_internal)
  call self%write_key_val("shiftvert", self%shiftvert, unit_internal)
  call self%write_key_val("theta2", self%theta2, unit_internal)
  call self%write_key_val("theta3", self%theta3, unit_internal)
  call self%write_key_val("thetad", self%thetad, unit_internal)
  call self%write_key_val("thetak", self%thetak, unit_internal)
  call self%write_key_val("thetam", self%thetam, unit_internal)
  call self%write_key_val("thetan", self%thetan, unit_internal)
  call self%write_key_val("tri", self%tri, unit_internal)
  call self%write_key_val("tripri", self%tripri, unit_internal)
  call self%write_key_val("zaxis", self%zaxis, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_theta_grid_parameters_config

!> Resets the config object to the initial empty state
subroutine reset_theta_grid_parameters_config(self)
  class(theta_grid_parameters_config_type), intent(in out) :: self
  type(theta_grid_parameters_config_type) :: empty
  select type (self)
  type is (theta_grid_parameters_config_type)
     self = empty
  end select
end subroutine reset_theta_grid_parameters_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_theta_grid_parameters_config(self)
  use mp, only: broadcast
  implicit none
  class(theta_grid_parameters_config_type), intent(in out) :: self
  call broadcast(self%akappa)
  call broadcast(self%akappri)
  call broadcast(self%alpmhd)
  call broadcast(self%asurf)
  call broadcast(self%btor_slab)
  call broadcast(self%c0_mxh)
  call broadcast(self%c_mxh)
  call broadcast(self%dc_mxh_dr)
  call broadcast(self%delta2)
  call broadcast(self%delta3)
  call broadcast(self%deltam)
  call broadcast(self%deltampri)
  call broadcast(self%deltan)
  call broadcast(self%deltanpri)
  call broadcast(self%ds_mxh_dr)
  call broadcast(self%eps)
  call broadcast(self%epsl)
  call broadcast(self%geotype)
  call broadcast(self%kp)
  call broadcast(self%mmode)
  call broadcast(self%n_mxh)
  call broadcast(self%nmode)
  call broadcast(self%nperiod)
  call broadcast(self%ntheta)
  call broadcast(self%pk)
  call broadcast(self%qinp)
  call broadcast(self%r_geo)
  call broadcast(self%raxis)
  call broadcast(self%rhoc)
  call broadcast(self%rmaj)
  call broadcast(self%s_mxh)
  call broadcast(self%shat)
  call broadcast(self%shift)
  call broadcast(self%shiftvert)
  call broadcast(self%theta2)
  call broadcast(self%theta3)
  call broadcast(self%thetad)
  call broadcast(self%thetak)
  call broadcast(self%thetam)
  call broadcast(self%thetan)
  call broadcast(self%tri)
  call broadcast(self%tripri)
  call broadcast(self%zaxis)

  call broadcast(self%exist)
end subroutine broadcast_theta_grid_parameters_config

!> Gets the default name for this namelist
function get_default_name_theta_grid_parameters_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_theta_grid_parameters_config
  get_default_name_theta_grid_parameters_config = "theta_grid_parameters"
end function get_default_name_theta_grid_parameters_config

!> Gets the default requires index for this namelist
function get_default_requires_index_theta_grid_parameters_config()
  implicit none
  logical :: get_default_requires_index_theta_grid_parameters_config
  get_default_requires_index_theta_grid_parameters_config = .false.
end function get_default_requires_index_theta_grid_parameters_config

!> Get the module level config instance
pure function get_theta_grid_parameters_config()
  type(theta_grid_parameters_config_type) :: get_theta_grid_parameters_config
  get_theta_grid_parameters_config = theta_grid_parameters_config
end function get_theta_grid_parameters_config
