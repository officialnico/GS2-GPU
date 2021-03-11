!> A module for handling the configuration of the diagnostics
!! module via the namelist diagnostics_config.
module diagnostics_config
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use simpledataio, only: sdf=>sdatio_file
  use diagnostics_ascii, only: diagnostics_ascii_type
  use diagnostics_dimensions, only: diagnostics_dimension_list_type
  use diagnostics_base_config, only: new_gs2_diagnostics_config_type => diagnostics_base_config_type
  implicit none

  private

  public :: init_diagnostics_config
  public :: finish_diagnostics_config
  public :: diagnostics_type
  public :: results_summary_type
  public :: override_screen_printout_options

  public :: new_gs2_diagnostics_config_type
  public :: set_diagnostics_config_config
  public :: get_diagnostics_config_config
  public :: apply_new_defaults

  real, parameter :: initial_value = -1.0
  
  !> A type for storing the current results of the simulation
  type results_summary_type
     real :: phi2 = initial_value
     real :: apar2 = initial_value
     real :: bpar2 = initial_value
     real :: total_heat_flux = initial_value
     real :: total_momentum_flux = initial_value
     real :: total_particle_flux = initial_value
     real :: max_growth_rate = initial_value
     real :: diffusivity = initial_value

     ! Individual heat fluxes
     real, dimension(:), pointer :: species_es_heat_flux
     real, dimension(:), pointer :: species_apar_heat_flux
     real, dimension(:), pointer :: species_bpar_heat_flux

     ! Total fluxes
     real, dimension(:), pointer :: species_heat_flux
     real, dimension(:), pointer :: species_momentum_flux
     real, dimension(:), pointer :: species_particle_flux
     real, dimension(:), pointer :: species_energy_exchange

     ! Average total fluxes
     real, dimension(:), pointer :: species_heat_flux_avg
     real, dimension(:), pointer :: species_momentum_flux_avg
     real, dimension(:), pointer :: species_particle_flux_avg

     ! Heating
     real, dimension(:), pointer :: species_heating
     real, dimension(:), pointer :: species_heating_avg

     ! Growth rates
     complex, dimension(:,:), pointer :: omega_average

     ! free energy transfer due to nonlinearity
     complex, dimension(:,:), pointer :: zonal_transfer

  end type results_summary_type

  !> A type for storing the diagnostics configuration,
  !! a reference to the output file, and current 
  !! results of the simulation
  type diagnostics_type
     type(sdf) :: sfile
     type(diagnostics_ascii_type) :: ascii_files
     type(results_summary_type) :: current_results
     type(diagnostics_dimension_list_type) :: dims
     !> Integer below gives the sdatio type 
     !! which corresponds to a gs2 real
     integer :: rtype
     integer :: itype
     integer :: ctype
     integer :: istep
     integer :: verbosity = 3
     logical :: create
     logical :: wryte
     logical :: reed
     logical :: replay
     logical :: distributed
     logical :: parallel
     logical :: exit
     logical :: vary_vnew_only
     logical :: calculate_fluxes
     logical :: is_trinity_run
     logical :: appending
     real :: user_time
     real :: user_time_old
     real :: start_time
     real, dimension(:), pointer :: fluxfac
     integer :: nwrite
     integer :: nwrite_mult
     logical :: write_any
     logical :: append_old
     logical :: enable_parallel
     logical :: serial_netcdf4
     integer :: igomega
     logical :: print_line
     logical :: print_flux_line
     logical :: write_line
     logical :: write_flux_line
     logical :: write_fields
     logical :: write_phi_over_time
     logical :: write_apar_over_time
     logical :: write_bpar_over_time
     logical :: make_movie
     logical :: dump_fields_periodically
     logical :: write_moments
     logical :: write_full_moments_notgc
     logical :: write_ntot_over_time
     logical :: write_density_over_time
     logical :: write_upar_over_time
     logical :: write_tperp_over_time
     logical :: write_fluxes
     logical :: write_fluxes_by_mode
     logical :: write_symmetry
     logical :: write_nl_flux_dist
     logical :: write_parity
     logical :: write_pflux_sym !Not implemented?
     logical :: write_pflux_tormom !Not implemented?
     logical :: write_omega
     integer :: navg
     real :: omegatinst
     real :: omegatol
     logical :: exit_when_converged
     logical :: write_verr
     logical :: write_cerr
     logical :: write_max_verr
     integer :: ncheck
     logical :: write_heating
     logical :: write_ascii
     logical :: write_gyx
     logical :: write_g
     logical :: write_lpoly
     integer :: conv_nstep_av
     real :: conv_test_multiplier
     integer :: conv_min_step
     integer :: conv_max_step
     integer :: conv_nsteps_converged
     logical :: use_nonlin_convergence
     logical :: write_cross_phase
     logical :: write_correlation
     logical :: write_correlation_extend
     logical :: write_jext
     logical :: write_lorentzian
     logical :: write_eigenfunc
     logical :: write_final_fields
     logical :: write_kpar
     logical :: write_final_epar
     logical :: write_final_db
     logical :: write_final_moments
     logical :: write_final_antot
     logical :: write_gs
     integer :: nsave
     logical :: save_for_restart
     logical :: save_many
     logical :: file_safety_check
     logical :: save_distfn
     logical :: save_glo_info_and_grids
     logical :: save_velocities
     logical :: write_zonal_transfer
     logical :: write_collisional
     logical :: write_omavg
     logical :: write_gg
     logical :: ob_midplane
     logical :: write_nl_flux
     logical :: write_hrate
     logical :: write_avg_moments
     logical :: dump_check1
     logical :: dump_check2
     integer :: nmovie
  end type diagnostics_type

  !> Used for testing... causes screen printout to be 
  !! generated regardless of the values of print_line 
  !! and print_flux_line if set to true
  logical :: override_screen_printout_options = .false.

  logical :: initialized = .false.

  type(new_gs2_diagnostics_config_type) :: new_gs2_diagnostics_config  
  
contains
  subroutine init_diagnostics_config(gnostics, new_gs2_diagnostics_config_in)
    use unit_tests, only: debug_message
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    type(new_gs2_diagnostics_config_type), intent(in), optional :: new_gs2_diagnostics_config_in
    if(initialized) return
    initialized = .true.
    call debug_message(3, 'diagnostics_config::init_diagnostics_config &
      & starting')
    call read_parameters(gnostics, new_gs2_diagnostics_config_in)
    call debug_message(3, 'diagnostics_config::init_diagnostics_config &
      & read_parameters')
    call allocate_current_results(gnostics)
  end subroutine init_diagnostics_config

  subroutine finish_diagnostics_config(gnostics)
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    initialized = .false.
    call deallocate_current_results(gnostics)
    call new_gs2_diagnostics_config%reset()
  end subroutine finish_diagnostics_config

  subroutine allocate_current_results(gnostics)
    use species, only: nspec
    use kt_grids, only: naky, ntheta0
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    allocate(gnostics%current_results%species_es_heat_flux(nspec))
    allocate(gnostics%current_results%species_apar_heat_flux(nspec))
    allocate(gnostics%current_results%species_bpar_heat_flux(nspec))
    allocate(gnostics%current_results%species_heat_flux(nspec))
    allocate(gnostics%current_results%species_momentum_flux(nspec))
    allocate(gnostics%current_results%species_particle_flux(nspec))
    allocate(gnostics%current_results%species_energy_exchange(nspec))
    allocate(gnostics%current_results%species_heat_flux_avg(nspec))
    allocate(gnostics%current_results%species_momentum_flux_avg(nspec))
    allocate(gnostics%current_results%species_particle_flux_avg(nspec))
    allocate(gnostics%current_results%species_heating(nspec))
    allocate(gnostics%current_results%species_heating_avg(nspec))
    allocate(gnostics%current_results%omega_average(ntheta0, naky))
    allocate(gnostics%current_results%zonal_transfer(ntheta0,naky))

  end subroutine allocate_current_results

  subroutine deallocate_current_results(gnostics)
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
   
    deallocate(gnostics%current_results%species_es_heat_flux)
    deallocate(gnostics%current_results%species_apar_heat_flux)
    deallocate(gnostics%current_results%species_bpar_heat_flux)
    deallocate(gnostics%current_results%species_heat_flux)
    deallocate(gnostics%current_results%species_momentum_flux)
    deallocate(gnostics%current_results%species_particle_flux)
    deallocate(gnostics%current_results%species_heat_flux_avg)
    deallocate(gnostics%current_results%species_momentum_flux_avg)
    deallocate(gnostics%current_results%species_particle_flux_avg)
    deallocate(gnostics%current_results%species_heating)
    deallocate(gnostics%current_results%species_heating_avg)
    deallocate(gnostics%current_results%omega_average)
    deallocate(gnostics%current_results%zonal_transfer)
  end subroutine deallocate_current_results

  !> Set values according to the new diagnostics defaults.
  !>
  !> Will overwrite existing values, unless the config object has
  !> already been initialised, in which case returns without changing
  !> any values
  subroutine apply_new_defaults(new_gs2_diagnostics_config_in)
    type(new_gs2_diagnostics_config_type), intent(inout) :: new_gs2_diagnostics_config_in

    if (new_gs2_diagnostics_config_in%is_initialised()) return

    new_gs2_diagnostics_config_in%navg = 10
    new_gs2_diagnostics_config_in%nmovie = -1
    new_gs2_diagnostics_config_in%nwrite = 10
    new_gs2_diagnostics_config_in%ob_midplane = .false.
    new_gs2_diagnostics_config_in%omegatinst = 1.0e6
    new_gs2_diagnostics_config_in%print_line = .false.
    new_gs2_diagnostics_config_in%write_correlation = .true.
    new_gs2_diagnostics_config_in%write_fields = .true.
    new_gs2_diagnostics_config_in%write_moments = .true.
  end subroutine apply_new_defaults

  subroutine read_parameters(gnostics, new_gs2_diagnostics_config_in)
    use diagnostics_base_config, only: warn_about_nonfunctional_selection
    use file_utils, only: input_unit, error_unit, input_unit_exist
    use text_options, only: text_option, get_option_value
    use mp, only: proc0, broadcast, nproc, iproc
    use nonlinear_terms, only: nonlin
    use unit_tests, only: debug_message
    implicit none
    type(diagnostics_type), intent(out) :: gnostics
    type(new_gs2_diagnostics_config_type), intent(in), optional :: new_gs2_diagnostics_config_in
    logical :: exist

    if (present(new_gs2_diagnostics_config_in)) new_gs2_diagnostics_config = new_gs2_diagnostics_config_in

    ! Smart defaults
    if (.not.new_gs2_diagnostics_config%is_initialised()) then
      if(nonlin) then
        new_gs2_diagnostics_config%write_fluxes = .true.
        new_gs2_diagnostics_config%write_fluxes_by_mode = .true.
      end if
      call apply_new_defaults(new_gs2_diagnostics_config)
    end if

    call new_gs2_diagnostics_config%init(name = 'gs2_diagnostics_knobs', requires_index = .false.)

    ! Print some health warnings if switches are not their default
    ! values and are not available in this diagnostics module
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%dump_check1, "dump_check1")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%dump_check2, "dump_check2")
    call warn_about_nonfunctional_selection(.not. new_gs2_diagnostics_config%file_safety_check, "file_safety_check")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%nc_sync_freq /= 10, "nc_sync_freq")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%ob_midplane, "ob_midplane")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%write_hrate, "write_hrate")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%write_nl_flux, "write_nl_flux")
    call warn_about_nonfunctional_selection(.not. new_gs2_diagnostics_config%write_omavg, "write_omavg")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%write_pflux_sym, "write_pflux_sym")
    call warn_about_nonfunctional_selection(new_gs2_diagnostics_config%write_pflux_tormom, "write_pflux_tormom")

    ! Copy out internal values into module level parameters
    gnostics%append_old = new_gs2_diagnostics_config%append_old
    gnostics%conv_max_step = new_gs2_diagnostics_config%conv_max_step
    gnostics%conv_min_step = new_gs2_diagnostics_config%conv_min_step
    gnostics%conv_nstep_av = new_gs2_diagnostics_config%conv_nstep_av
    gnostics%conv_nsteps_converged = new_gs2_diagnostics_config%conv_nsteps_converged
    gnostics%conv_test_multiplier = new_gs2_diagnostics_config%conv_test_multiplier
    gnostics%dump_check1 = new_gs2_diagnostics_config%dump_check1
    gnostics%dump_check2 = new_gs2_diagnostics_config%dump_check2
    gnostics%dump_fields_periodically = new_gs2_diagnostics_config%dump_fields_periodically
    gnostics%enable_parallel = new_gs2_diagnostics_config%enable_parallel
    gnostics%exit_when_converged = new_gs2_diagnostics_config%exit_when_converged
    gnostics%file_safety_check = new_gs2_diagnostics_config%file_safety_check
    gnostics%igomega = new_gs2_diagnostics_config%igomega
    gnostics%make_movie = new_gs2_diagnostics_config%make_movie
    gnostics%navg = new_gs2_diagnostics_config%navg
    gnostics%ncheck = new_gs2_diagnostics_config%ncheck
    gnostics%nmovie = new_gs2_diagnostics_config%nmovie
    gnostics%nsave = new_gs2_diagnostics_config%nsave
    gnostics%nwrite = new_gs2_diagnostics_config%nwrite
    gnostics%nwrite_mult = new_gs2_diagnostics_config%nwrite_mult
    gnostics%ob_midplane = new_gs2_diagnostics_config%ob_midplane
    gnostics%omegatinst = new_gs2_diagnostics_config%omegatinst
    gnostics%omegatol = new_gs2_diagnostics_config%omegatol
    gnostics%print_flux_line = new_gs2_diagnostics_config%print_flux_line
    gnostics%print_line = new_gs2_diagnostics_config%print_line
    gnostics%save_distfn = new_gs2_diagnostics_config%save_distfn
    gnostics%save_for_restart = new_gs2_diagnostics_config%save_for_restart
    gnostics%save_glo_info_and_grids = new_gs2_diagnostics_config%save_glo_info_and_grids
    gnostics%save_many = new_gs2_diagnostics_config%save_many
    gnostics%save_velocities = new_gs2_diagnostics_config%save_velocities
    gnostics%serial_netcdf4 = new_gs2_diagnostics_config%serial_netcdf4
    gnostics%use_nonlin_convergence = new_gs2_diagnostics_config%use_nonlin_convergence
    gnostics%write_any = new_gs2_diagnostics_config%write_any
    gnostics%write_apar_over_time = new_gs2_diagnostics_config%write_apar_over_time
    gnostics%write_ascii = new_gs2_diagnostics_config%write_ascii
    gnostics%write_avg_moments = new_gs2_diagnostics_config%write_avg_moments
    gnostics%write_bpar_over_time = new_gs2_diagnostics_config%write_bpar_over_time
    gnostics%write_cerr = new_gs2_diagnostics_config%write_cerr
    gnostics%write_collisional = new_gs2_diagnostics_config%write_collisional
    gnostics%write_correlation = new_gs2_diagnostics_config%write_correlation
    gnostics%write_correlation_extend = new_gs2_diagnostics_config%write_correlation_extend
    gnostics%write_cross_phase = new_gs2_diagnostics_config%write_cross_phase
    gnostics%write_density_over_time = new_gs2_diagnostics_config%write_density_over_time
    gnostics%write_eigenfunc = new_gs2_diagnostics_config%write_eigenfunc
    gnostics%write_fields = new_gs2_diagnostics_config%write_fields
    gnostics%write_final_antot = new_gs2_diagnostics_config%write_final_antot
    gnostics%write_final_db = new_gs2_diagnostics_config%write_final_db
    gnostics%write_final_epar = new_gs2_diagnostics_config%write_final_epar
    gnostics%write_final_fields = new_gs2_diagnostics_config%write_final_fields
    gnostics%write_final_moments = new_gs2_diagnostics_config%write_final_moments
    gnostics%write_flux_line = new_gs2_diagnostics_config%write_flux_line
    gnostics%write_fluxes = new_gs2_diagnostics_config%write_fluxes
    gnostics%write_fluxes_by_mode = new_gs2_diagnostics_config%write_fluxes_by_mode
    gnostics%write_full_moments_notgc = new_gs2_diagnostics_config%write_full_moments_notgc
    gnostics%write_g = new_gs2_diagnostics_config%write_g
    gnostics%write_gg = new_gs2_diagnostics_config%write_gg
    gnostics%write_gs = new_gs2_diagnostics_config%write_gs
    gnostics%write_gyx = new_gs2_diagnostics_config%write_gyx
    gnostics%write_heating = new_gs2_diagnostics_config%write_heating
    gnostics%write_hrate = new_gs2_diagnostics_config%write_hrate
    gnostics%write_jext = new_gs2_diagnostics_config%write_jext
    gnostics%write_kpar = new_gs2_diagnostics_config%write_kpar
    gnostics%write_line = new_gs2_diagnostics_config%write_line
    gnostics%write_lorentzian = new_gs2_diagnostics_config%write_lorentzian
    gnostics%write_lpoly = new_gs2_diagnostics_config%write_lpoly
    gnostics%write_max_verr = new_gs2_diagnostics_config%write_max_verr
    gnostics%write_moments = new_gs2_diagnostics_config%write_moments
    gnostics%write_nl_flux = new_gs2_diagnostics_config%write_nl_flux
    gnostics%write_nl_flux_dist = new_gs2_diagnostics_config%write_nl_flux_dist
    gnostics%write_ntot_over_time = new_gs2_diagnostics_config%write_ntot_over_time
    gnostics%write_omavg = new_gs2_diagnostics_config%write_omavg
    gnostics%write_omega = new_gs2_diagnostics_config%write_omega
    gnostics%write_parity = new_gs2_diagnostics_config%write_parity
    gnostics%write_pflux_sym = new_gs2_diagnostics_config%write_pflux_sym
    gnostics%write_pflux_tormom = new_gs2_diagnostics_config%write_pflux_tormom
    gnostics%write_phi_over_time = new_gs2_diagnostics_config%write_phi_over_time
    gnostics%write_symmetry = new_gs2_diagnostics_config%write_symmetry
    gnostics%write_tperp_over_time = new_gs2_diagnostics_config%write_tperp_over_time
    gnostics%write_upar_over_time = new_gs2_diagnostics_config%write_upar_over_time
    gnostics%write_verr = new_gs2_diagnostics_config%write_verr
    gnostics%write_zonal_transfer = new_gs2_diagnostics_config%write_zonal_transfer

    exist = new_gs2_diagnostics_config%exist

    if (.not.(gnostics%save_for_restart.or.gnostics%save_distfn)) then
       gnostics%nsave = -1
    endif

    if (override_screen_printout_options) then 
       gnostics%print_line = .true.
       gnostics%print_flux_line = .true.
    end if
  end subroutine read_parameters

  !> Set the module level config type
  !> Will abort if the module has already been initialised to avoid
  !> inconsistencies.
  subroutine set_diagnostics_config_config(diagnostics_config_config_in)
    use mp, only: mp_abort
    type(new_gs2_diagnostics_config_type), intent(in), optional :: diagnostics_config_config_in
    if (initialized) then
       call mp_abort("Trying to set diagnostics_config config when already initialized.", to_screen = .true.)
    end if
    if (present(diagnostics_config_config_in)) then
       new_gs2_diagnostics_config = diagnostics_config_config_in
    end if
  end subroutine set_diagnostics_config_config

  !> Get the module level config instance
  function get_diagnostics_config_config()
    type(new_gs2_diagnostics_config_type) :: get_diagnostics_config_config
    get_diagnostics_config_config = new_gs2_diagnostics_config
  end function get_diagnostics_config_config

end module diagnostics_config
