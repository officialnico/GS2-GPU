
!---------------------------------------
! Following is for the config_type
!---------------------------------------
        
!> Reads in the gs2_diagnostics_knobs namelist and populates the member variables
subroutine read_diagnostics_base_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class(diagnostics_base_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
  integer :: conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, igomega, navg, nc_sync_freq, ncheck, nmovie, nsave, nwrite, nwrite_mult
  logical :: append_old, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, file_safety_check, make_movie, ob_midplane, print_flux_line, print_line, save_distfn
  logical :: save_for_restart, save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, write_bpar_over_time
  logical :: write_cerr, write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db
  logical :: write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gs, write_gyx, write_heating, write_jext
  logical :: write_kpar, write_line, write_lorentzian, write_max_verr, write_moments, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega, write_parity, write_pflux_sym, write_pflux_tormom
  logical :: write_phi_over_time, write_ql_metric, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer
  real :: conv_test_multiplier, omegatinst, omegatol

  namelist /gs2_diagnostics_knobs/ append_old, conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, conv_test_multiplier, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, &
    file_safety_check, igomega, make_movie, navg, nc_sync_freq, ncheck, nmovie, nsave, nwrite, nwrite_mult, ob_midplane, omegatinst, omegatol, print_flux_line, print_line, save_distfn, &
    save_for_restart, save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, &
    write_bpar_over_time, write_cerr, write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, &
    write_final_db, write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gs, write_gyx, &
    write_heating, write_jext, write_kpar, write_line, write_lorentzian, write_max_verr, write_moments, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega, write_parity, &
    write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_ql_metric, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
  append_old = self%append_old
  conv_max_step = self%conv_max_step
  conv_min_step = self%conv_min_step
  conv_nstep_av = self%conv_nstep_av
  conv_nsteps_converged = self%conv_nsteps_converged
  conv_test_multiplier = self%conv_test_multiplier
  dump_check1 = self%dump_check1
  dump_check2 = self%dump_check2
  dump_fields_periodically = self%dump_fields_periodically
  enable_parallel = self%enable_parallel
  exit_when_converged = self%exit_when_converged
  file_safety_check = self%file_safety_check
  igomega = self%igomega
  make_movie = self%make_movie
  navg = self%navg
  nc_sync_freq = self%nc_sync_freq
  ncheck = self%ncheck
  nmovie = self%nmovie
  nsave = self%nsave
  nwrite = self%nwrite
  nwrite_mult = self%nwrite_mult
  ob_midplane = self%ob_midplane
  omegatinst = self%omegatinst
  omegatol = self%omegatol
  print_flux_line = self%print_flux_line
  print_line = self%print_line
  save_distfn = self%save_distfn
  save_for_restart = self%save_for_restart
  save_glo_info_and_grids = self%save_glo_info_and_grids
  save_many = self%save_many
  save_velocities = self%save_velocities
  serial_netcdf4 = self%serial_netcdf4
  use_nonlin_convergence = self%use_nonlin_convergence
  write_any = self%write_any
  write_apar_over_time = self%write_apar_over_time
  write_ascii = self%write_ascii
  write_avg_moments = self%write_avg_moments
  write_bpar_over_time = self%write_bpar_over_time
  write_cerr = self%write_cerr
  write_collisional = self%write_collisional
  write_correlation = self%write_correlation
  write_correlation_extend = self%write_correlation_extend
  write_cross_phase = self%write_cross_phase
  write_density_over_time = self%write_density_over_time
  write_eigenfunc = self%write_eigenfunc
  write_fields = self%write_fields
  write_final_antot = self%write_final_antot
  write_final_db = self%write_final_db
  write_final_epar = self%write_final_epar
  write_final_fields = self%write_final_fields
  write_final_moments = self%write_final_moments
  write_flux_line = self%write_flux_line
  write_fluxes = self%write_fluxes
  write_fluxes_by_mode = self%write_fluxes_by_mode
  write_full_moments_notgc = self%write_full_moments_notgc
  write_g = self%write_g
  write_gs = self%write_gs
  write_gyx = self%write_gyx
  write_heating = self%write_heating
  write_jext = self%write_jext
  write_kpar = self%write_kpar
  write_line = self%write_line
  write_lorentzian = self%write_lorentzian
  write_max_verr = self%write_max_verr
  write_moments = self%write_moments
  write_nl_flux_dist = self%write_nl_flux_dist
  write_ntot_over_time = self%write_ntot_over_time
  write_omavg = self%write_omavg
  write_omega = self%write_omega
  write_parity = self%write_parity
  write_pflux_sym = self%write_pflux_sym
  write_pflux_tormom = self%write_pflux_tormom
  write_phi_over_time = self%write_phi_over_time
  write_ql_metric = self%write_ql_metric
  write_symmetry = self%write_symmetry
  write_tperp_over_time = self%write_tperp_over_time
  write_upar_over_time = self%write_upar_over_time
  write_verr = self%write_verr
  write_zonal_transfer = self%write_zonal_transfer

  ! Now read in the main namelist
  in_file = input_unit_exist(self%get_name(), exist)
  if (exist) read(in_file, nml = gs2_diagnostics_knobs)

  ! Now copy from local variables into type members
  self%append_old = append_old
  self%conv_max_step = conv_max_step
  self%conv_min_step = conv_min_step
  self%conv_nstep_av = conv_nstep_av
  self%conv_nsteps_converged = conv_nsteps_converged
  self%conv_test_multiplier = conv_test_multiplier
  self%dump_check1 = dump_check1
  self%dump_check2 = dump_check2
  self%dump_fields_periodically = dump_fields_periodically
  self%enable_parallel = enable_parallel
  self%exit_when_converged = exit_when_converged
  self%file_safety_check = file_safety_check
  self%igomega = igomega
  self%make_movie = make_movie
  self%navg = navg
  self%nc_sync_freq = nc_sync_freq
  self%ncheck = ncheck
  self%nmovie = nmovie
  self%nsave = nsave
  self%nwrite = nwrite
  self%nwrite_mult = nwrite_mult
  self%ob_midplane = ob_midplane
  self%omegatinst = omegatinst
  self%omegatol = omegatol
  self%print_flux_line = print_flux_line
  self%print_line = print_line
  self%save_distfn = save_distfn
  self%save_for_restart = save_for_restart
  self%save_glo_info_and_grids = save_glo_info_and_grids
  self%save_many = save_many
  self%save_velocities = save_velocities
  self%serial_netcdf4 = serial_netcdf4
  self%use_nonlin_convergence = use_nonlin_convergence
  self%write_any = write_any
  self%write_apar_over_time = write_apar_over_time
  self%write_ascii = write_ascii
  self%write_avg_moments = write_avg_moments
  self%write_bpar_over_time = write_bpar_over_time
  self%write_cerr = write_cerr
  self%write_collisional = write_collisional
  self%write_correlation = write_correlation
  self%write_correlation_extend = write_correlation_extend
  self%write_cross_phase = write_cross_phase
  self%write_density_over_time = write_density_over_time
  self%write_eigenfunc = write_eigenfunc
  self%write_fields = write_fields
  self%write_final_antot = write_final_antot
  self%write_final_db = write_final_db
  self%write_final_epar = write_final_epar
  self%write_final_fields = write_final_fields
  self%write_final_moments = write_final_moments
  self%write_flux_line = write_flux_line
  self%write_fluxes = write_fluxes
  self%write_fluxes_by_mode = write_fluxes_by_mode
  self%write_full_moments_notgc = write_full_moments_notgc
  self%write_g = write_g
  self%write_gs = write_gs
  self%write_gyx = write_gyx
  self%write_heating = write_heating
  self%write_jext = write_jext
  self%write_kpar = write_kpar
  self%write_line = write_line
  self%write_lorentzian = write_lorentzian
  self%write_max_verr = write_max_verr
  self%write_moments = write_moments
  self%write_nl_flux_dist = write_nl_flux_dist
  self%write_ntot_over_time = write_ntot_over_time
  self%write_omavg = write_omavg
  self%write_omega = write_omega
  self%write_parity = write_parity
  self%write_pflux_sym = write_pflux_sym
  self%write_pflux_tormom = write_pflux_tormom
  self%write_phi_over_time = write_phi_over_time
  self%write_ql_metric = write_ql_metric
  self%write_symmetry = write_symmetry
  self%write_tperp_over_time = write_tperp_over_time
  self%write_upar_over_time = write_upar_over_time
  self%write_verr = write_verr
  self%write_zonal_transfer = write_zonal_transfer

  self%exist = exist
end subroutine read_diagnostics_base_config

!> Writes out a namelist representing the current state of the config object
subroutine write_diagnostics_base_config(self, unit)
  implicit none
  class(diagnostics_base_config_type), intent(in) :: self
  integer, intent(in), optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
  call self%write_key_val("append_old", self%append_old, unit_internal)
  call self%write_key_val("conv_max_step", self%conv_max_step, unit_internal)
  call self%write_key_val("conv_min_step", self%conv_min_step, unit_internal)
  call self%write_key_val("conv_nstep_av", self%conv_nstep_av, unit_internal)
  call self%write_key_val("conv_nsteps_converged", self%conv_nsteps_converged, unit_internal)
  call self%write_key_val("conv_test_multiplier", self%conv_test_multiplier, unit_internal)
  call self%write_key_val("dump_check1", self%dump_check1, unit_internal)
  call self%write_key_val("dump_check2", self%dump_check2, unit_internal)
  call self%write_key_val("dump_fields_periodically", self%dump_fields_periodically, unit_internal)
  call self%write_key_val("enable_parallel", self%enable_parallel, unit_internal)
  call self%write_key_val("exit_when_converged", self%exit_when_converged, unit_internal)
  call self%write_key_val("file_safety_check", self%file_safety_check, unit_internal)
  call self%write_key_val("igomega", self%igomega, unit_internal)
  call self%write_key_val("make_movie", self%make_movie, unit_internal)
  call self%write_key_val("navg", self%navg, unit_internal)
  call self%write_key_val("nc_sync_freq", self%nc_sync_freq, unit_internal)
  call self%write_key_val("ncheck", self%ncheck, unit_internal)
  call self%write_key_val("nmovie", self%nmovie, unit_internal)
  call self%write_key_val("nsave", self%nsave, unit_internal)
  call self%write_key_val("nwrite", self%nwrite, unit_internal)
  call self%write_key_val("nwrite_mult", self%nwrite_mult, unit_internal)
  call self%write_key_val("ob_midplane", self%ob_midplane, unit_internal)
  call self%write_key_val("omegatinst", self%omegatinst, unit_internal)
  call self%write_key_val("omegatol", self%omegatol, unit_internal)
  call self%write_key_val("print_flux_line", self%print_flux_line, unit_internal)
  call self%write_key_val("print_line", self%print_line, unit_internal)
  call self%write_key_val("save_distfn", self%save_distfn, unit_internal)
  call self%write_key_val("save_for_restart", self%save_for_restart, unit_internal)
  call self%write_key_val("save_glo_info_and_grids", self%save_glo_info_and_grids, unit_internal)
  call self%write_key_val("save_many", self%save_many, unit_internal)
  call self%write_key_val("save_velocities", self%save_velocities, unit_internal)
  call self%write_key_val("serial_netcdf4", self%serial_netcdf4, unit_internal)
  call self%write_key_val("use_nonlin_convergence", self%use_nonlin_convergence, unit_internal)
  call self%write_key_val("write_any", self%write_any, unit_internal)
  call self%write_key_val("write_apar_over_time", self%write_apar_over_time, unit_internal)
  call self%write_key_val("write_ascii", self%write_ascii, unit_internal)
  call self%write_key_val("write_avg_moments", self%write_avg_moments, unit_internal)
  call self%write_key_val("write_bpar_over_time", self%write_bpar_over_time, unit_internal)
  call self%write_key_val("write_cerr", self%write_cerr, unit_internal)
  call self%write_key_val("write_collisional", self%write_collisional, unit_internal)
  call self%write_key_val("write_correlation", self%write_correlation, unit_internal)
  call self%write_key_val("write_correlation_extend", self%write_correlation_extend, unit_internal)
  call self%write_key_val("write_cross_phase", self%write_cross_phase, unit_internal)
  call self%write_key_val("write_density_over_time", self%write_density_over_time, unit_internal)
  call self%write_key_val("write_eigenfunc", self%write_eigenfunc, unit_internal)
  call self%write_key_val("write_fields", self%write_fields, unit_internal)
  call self%write_key_val("write_final_antot", self%write_final_antot, unit_internal)
  call self%write_key_val("write_final_db", self%write_final_db, unit_internal)
  call self%write_key_val("write_final_epar", self%write_final_epar, unit_internal)
  call self%write_key_val("write_final_fields", self%write_final_fields, unit_internal)
  call self%write_key_val("write_final_moments", self%write_final_moments, unit_internal)
  call self%write_key_val("write_flux_line", self%write_flux_line, unit_internal)
  call self%write_key_val("write_fluxes", self%write_fluxes, unit_internal)
  call self%write_key_val("write_fluxes_by_mode", self%write_fluxes_by_mode, unit_internal)
  call self%write_key_val("write_full_moments_notgc", self%write_full_moments_notgc, unit_internal)
  call self%write_key_val("write_g", self%write_g, unit_internal)
  call self%write_key_val("write_gs", self%write_gs, unit_internal)
  call self%write_key_val("write_gyx", self%write_gyx, unit_internal)
  call self%write_key_val("write_heating", self%write_heating, unit_internal)
  call self%write_key_val("write_jext", self%write_jext, unit_internal)
  call self%write_key_val("write_kpar", self%write_kpar, unit_internal)
  call self%write_key_val("write_line", self%write_line, unit_internal)
  call self%write_key_val("write_lorentzian", self%write_lorentzian, unit_internal)
  call self%write_key_val("write_max_verr", self%write_max_verr, unit_internal)
  call self%write_key_val("write_moments", self%write_moments, unit_internal)
  call self%write_key_val("write_nl_flux_dist", self%write_nl_flux_dist, unit_internal)
  call self%write_key_val("write_ntot_over_time", self%write_ntot_over_time, unit_internal)
  call self%write_key_val("write_omavg", self%write_omavg, unit_internal)
  call self%write_key_val("write_omega", self%write_omega, unit_internal)
  call self%write_key_val("write_parity", self%write_parity, unit_internal)
  call self%write_key_val("write_pflux_sym", self%write_pflux_sym, unit_internal)
  call self%write_key_val("write_pflux_tormom", self%write_pflux_tormom, unit_internal)
  call self%write_key_val("write_phi_over_time", self%write_phi_over_time, unit_internal)
  call self%write_key_val("write_ql_metric", self%write_ql_metric, unit_internal)
  call self%write_key_val("write_symmetry", self%write_symmetry, unit_internal)
  call self%write_key_val("write_tperp_over_time", self%write_tperp_over_time, unit_internal)
  call self%write_key_val("write_upar_over_time", self%write_upar_over_time, unit_internal)
  call self%write_key_val("write_verr", self%write_verr, unit_internal)
  call self%write_key_val("write_zonal_transfer", self%write_zonal_transfer, unit_internal)
  call self%write_namelist_footer(unit_internal)
end subroutine write_diagnostics_base_config

!> Resets the config object to the initial empty state
subroutine reset_diagnostics_base_config(self)
  class(diagnostics_base_config_type), intent(in out) :: self
  type(diagnostics_base_config_type) :: empty
  select type (self)
  type is (diagnostics_base_config_type)
     self = empty
  end select
end subroutine reset_diagnostics_base_config

!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_diagnostics_base_config(self)
  use mp, only: broadcast
  implicit none
  class(diagnostics_base_config_type), intent(in out) :: self
  call broadcast(self%append_old)
  call broadcast(self%conv_max_step)
  call broadcast(self%conv_min_step)
  call broadcast(self%conv_nstep_av)
  call broadcast(self%conv_nsteps_converged)
  call broadcast(self%conv_test_multiplier)
  call broadcast(self%dump_check1)
  call broadcast(self%dump_check2)
  call broadcast(self%dump_fields_periodically)
  call broadcast(self%enable_parallel)
  call broadcast(self%exit_when_converged)
  call broadcast(self%file_safety_check)
  call broadcast(self%igomega)
  call broadcast(self%make_movie)
  call broadcast(self%navg)
  call broadcast(self%nc_sync_freq)
  call broadcast(self%ncheck)
  call broadcast(self%nmovie)
  call broadcast(self%nsave)
  call broadcast(self%nwrite)
  call broadcast(self%nwrite_mult)
  call broadcast(self%ob_midplane)
  call broadcast(self%omegatinst)
  call broadcast(self%omegatol)
  call broadcast(self%print_flux_line)
  call broadcast(self%print_line)
  call broadcast(self%save_distfn)
  call broadcast(self%save_for_restart)
  call broadcast(self%save_glo_info_and_grids)
  call broadcast(self%save_many)
  call broadcast(self%save_velocities)
  call broadcast(self%serial_netcdf4)
  call broadcast(self%use_nonlin_convergence)
  call broadcast(self%write_any)
  call broadcast(self%write_apar_over_time)
  call broadcast(self%write_ascii)
  call broadcast(self%write_avg_moments)
  call broadcast(self%write_bpar_over_time)
  call broadcast(self%write_cerr)
  call broadcast(self%write_collisional)
  call broadcast(self%write_correlation)
  call broadcast(self%write_correlation_extend)
  call broadcast(self%write_cross_phase)
  call broadcast(self%write_density_over_time)
  call broadcast(self%write_eigenfunc)
  call broadcast(self%write_fields)
  call broadcast(self%write_final_antot)
  call broadcast(self%write_final_db)
  call broadcast(self%write_final_epar)
  call broadcast(self%write_final_fields)
  call broadcast(self%write_final_moments)
  call broadcast(self%write_flux_line)
  call broadcast(self%write_fluxes)
  call broadcast(self%write_fluxes_by_mode)
  call broadcast(self%write_full_moments_notgc)
  call broadcast(self%write_g)
  call broadcast(self%write_gs)
  call broadcast(self%write_gyx)
  call broadcast(self%write_heating)
  call broadcast(self%write_jext)
  call broadcast(self%write_kpar)
  call broadcast(self%write_line)
  call broadcast(self%write_lorentzian)
  call broadcast(self%write_max_verr)
  call broadcast(self%write_moments)
  call broadcast(self%write_nl_flux_dist)
  call broadcast(self%write_ntot_over_time)
  call broadcast(self%write_omavg)
  call broadcast(self%write_omega)
  call broadcast(self%write_parity)
  call broadcast(self%write_pflux_sym)
  call broadcast(self%write_pflux_tormom)
  call broadcast(self%write_phi_over_time)
  call broadcast(self%write_ql_metric)
  call broadcast(self%write_symmetry)
  call broadcast(self%write_tperp_over_time)
  call broadcast(self%write_upar_over_time)
  call broadcast(self%write_verr)
  call broadcast(self%write_zonal_transfer)

  call broadcast(self%exist)
end subroutine broadcast_diagnostics_base_config

!> Gets the default name for this namelist
function get_default_name_diagnostics_base_config()
  implicit none
  character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_diagnostics_base_config
  get_default_name_diagnostics_base_config = "gs2_diagnostics_knobs"
end function get_default_name_diagnostics_base_config

!> Gets the default requires index for this namelist
function get_default_requires_index_diagnostics_base_config()
  implicit none
  logical :: get_default_requires_index_diagnostics_base_config
  get_default_requires_index_diagnostics_base_config = .false.
end function get_default_requires_index_diagnostics_base_config

!> Get the module level config instance
pure function get_diagnostics_base_config()
  type(diagnostics_base_config_type) :: get_diagnostics_base_config
  get_diagnostics_base_config = diagnostics_base_config
end function get_diagnostics_base_config
