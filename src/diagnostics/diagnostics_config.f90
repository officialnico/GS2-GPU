!> A module for handling the configuration of the diagnostics
!! module via the namelist diagnostics_config.
module diagnostics_config
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  use simpledataio, only: sdf=>sdatio_file
  use diagnostics_ascii, only: diagnostics_ascii_type
  use diagnostics_dimensions, only: diagnostics_dimension_list_type
  implicit none

  private

  public :: init_diagnostics_config
  public :: finish_diagnostics_config
  public :: diagnostics_type
  public :: results_summary_type
  public :: override_screen_printout_options

  public :: new_gs2_diagnostics_config_type
  public :: set_diagnostics_config_config
  
  !> A type for storing the current results of the simulation
  type results_summary_type
     real :: phi2
     real :: apar2
     real :: bpar2
     real :: total_heat_flux
     real :: total_momentum_flux
     real :: total_particle_flux
     real :: max_growth_rate
     real :: diffusivity

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
  
  !> Used to represent the input configuration of new_gs2_diagnostics
  !> Very similar to original diagnostics config but with some new inputs
  !> and several different default values. Really would be good to just
  !> have a single diagnostics config type.
  type, extends(abstract_config_type) :: new_gs2_diagnostics_config_type
     ! namelist : gs2_diagnostics_knobs
     ! indexed : false     
     !> If true, and a copy of the .out.nc already exists,
     !> open it and append to it. Obviously grid sizes must be
     !> unchanged.
     logical :: append_old = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     integer :: conv_max_step = 80000
     !> FIXME: Add documentation
     integer :: conv_min_step = 4000
     !> FIXME: Add documentation
     integer :: conv_nstep_av = 4000
     !> FIXME: Add documentation
     integer :: conv_nsteps_converged = 10000
     !> FIXME: Add documentation
     real :: conv_test_multiplier = 4e-1
     !> FIXME: Add documentation
     logical :: dump_check1 = .false.
     !> FIXME: Add documentation
     logical :: dump_check2 = .false.
     !> FIXME: Add documentation
     logical :: dump_fields_periodically = .false.
     !> If built with parallel IO capability
     !> enable parallel IO. Currently disabled by
     !> default because of problems on Helios:
     !> parallel IO seems to work on 1 node but not
     !> on >1. This problem does not affect parallel
     !> IO for the single restart file because it only
     !> affects output with infinite dimensions
     logical :: enable_parallel = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: exit_when_converged = .true.
     !> FIXME: Add documentation
     logical :: file_safety_check = .true.
     !> FIXME: Add documentation
     integer :: igomega = 0
     !> FIXME: Add documentation
     logical :: make_movie = .false.
     !> FIXME: Add documentation
     integer :: navg = 10
     !> FIXME: Add documentation
     integer :: ncheck = 10 ! Only in new diagnostics
     !> FIXME: Add documentation
     integer :: nmovie = -1
     !> FIXME: Add documentation
     integer :: nsave = -1
     !> FIXME: Add documentation
     integer :: nwrite = 10
     !> FIXME: Add documentation
     integer :: nwrite_mult = 10
     !> FIXME: Add documentation
     logical :: ob_midplane = .false.
     !> FIXME: Add documentation
     real :: omegatinst = 1.0e6
     !> FIXME: Add documentation
     real :: omegatol = 1e-3
     !> FIXME: Add documentation
     logical :: print_flux_line = .false.
     !> FIXME: Add documentation
     logical :: print_line = .false.
     !> FIXME: Add documentation
     logical :: save_distfn = .false.
     !> FIXME: Add documentation
     logical :: save_for_restart = .false.
     !> FIXME: Add documentation
     logical :: save_glo_info_and_grids = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: save_many = .false.
     !> FIXME: Add documentation
     logical :: save_velocities = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: serial_netcdf4 = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: use_nonlin_convergence = .false.
     !> FIXME: Add documentation
     logical :: write_any = .true. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_apar_over_time = .false.
     !> FIXME: Add documentation
     logical :: write_ascii = .true.
     !> FIXME: Add documentation
     logical :: write_avg_moments = .false.
     !> FIXME: Add documentation
     logical :: write_bpar_over_time = .false.
     !> FIXME: Add documentation
     logical :: write_cerr = .false.
     !> FIXME: Add documentation
     logical :: write_collisional = .false. ! Only in new diagnostics     
     !> FIXME: Add documentation
     logical :: write_correlation = .true.
     !> FIXME: Add documentation
     logical :: write_correlation_extend = .false.
     !> FIXME: Add documentation
     logical :: write_cross_phase = .false.
     !> FIXME: Add documentation
     logical :: write_density_over_time = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_eigenfunc = .false.
     !> FIXME: Add documentation
     logical :: write_fields = .true.
     !> FIXME: Add documentation
     logical :: write_final_antot = .false.
     !> FIXME: Add documentation
     logical :: write_final_db = .false.
     !> FIXME: Add documentation
     logical :: write_final_epar = .false.
     !> FIXME: Add documentation
     logical :: write_final_fields = .false.
     !> FIXME: Add documentation
     logical :: write_final_moments = .false.
     !> FIXME: Add documentation
     logical :: write_flux_line = .true.
     !> Write fluxes over time
     !> This gets a smart default
     !> If run is nonlinear, this defaults to true
     logical :: write_fluxes = .false. ! Only in new diagnostics
     !> Write fluxes over time as a function of Fourier mode
     !> This gets a smart default
     !> If run is nonlinear, this defaults to true
     logical :: write_fluxes_by_mode = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_full_moments_notgc = .false.
     !> FIXME: Add documentation
     logical :: write_g = .false.
     !> FIXME: Add documentation
     logical :: write_gg = .false.
     !> FIXME: Add documentation
     logical :: write_gs = .false.
     !> FIXME: Add documentation
     logical :: write_gyx = .false.
     !> FIXME: Add documentation
     logical :: write_heating = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_hrate = .false.
     !> FIXME: Add documentation
     logical :: write_jext = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_kpar = .false.
     !> FIXME: Add documentation
     logical :: write_line = .true.
     !> FIXME: Add documentation
     logical :: write_lorentzian = .false.
     !> FIXME: Add documentation
     logical :: write_lpoly = .false.
     !> FIXME: Add documentation
     logical :: write_max_verr = .false.
     !> FIXME: Add documentation
     logical :: write_moments = .true.
     !> FIXME: Add documentation
     logical :: write_nl_flux = .true.
     !> FIXME: Add documentation
     logical :: write_nl_flux_dist = .false.
     !> FIXME: Add documentation
     logical :: write_ntot_over_time = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_omavg = .true.
     !> FIXME: Add documentation
     logical :: write_omega = .true.
     !> FIXME: Add documentation
     logical :: write_parity = .false.
     !> FIXME: Add documentation
     logical :: write_pflux_sym = .false. 
     !> FIXME: Add documentation
     logical :: write_pflux_tormom = .false. 
     !> FIXME: Add documentation
     logical :: write_phi_over_time = .false.
     !> FIXME: Add documentation
     logical :: write_symmetry = .false.
     !> FIXME: Add documentation
     logical :: write_tperp_over_time = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_upar_over_time = .false. ! Only in new diagnostics
     !> FIXME: Add documentation
     logical :: write_verr = .false.
     !> FIXME: Add documentation
     logical :: write_zonal_transfer = .false. ! Only in new diagnostics
   contains
     procedure, public :: read => read_new_gs2_diagnostics_config
     procedure, public :: write => write_new_gs2_diagnostics_config
     procedure, public :: reset => reset_new_gs2_diagnostics_config
     procedure, public :: broadcast => broadcast_new_gs2_diagnostics_config
     procedure, public, nopass :: get_default_name => get_default_name_new_gs2_diagnostics_config
     procedure, public, nopass :: get_default_requires_index => get_default_requires_index_new_gs2_diagnostics_config
  end type new_gs2_diagnostics_config_type
  
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

  subroutine read_parameters(gnostics, new_gs2_diagnostics_config_in)
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
    end if

    call new_gs2_diagnostics_config%init(name = 'gs2_diagnostics_knobs', requires_index = .false.)

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
  
  !---------------------------------------
  ! Following is for the config_type
  !---------------------------------------
  
  !> Reads in the gs2_diagnostics_knobs namelist and populates the member variables
  subroutine read_new_gs2_diagnostics_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(new_gs2_diagnostics_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.    
    logical :: append_old, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, file_safety_check, make_movie, ob_midplane, print_flux_line, print_line, save_distfn, save_for_restart, save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, write_bpar_over_time, write_cerr, write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db, write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating, write_hrate, write_jext, write_kpar, write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega, write_parity, write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer
    real :: conv_test_multiplier, omegatinst, omegatol
    integer :: conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, igomega, navg, ncheck, nmovie, nsave, nwrite, nwrite_mult

    namelist /gs2_diagnostics_knobs/ append_old, conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, conv_test_multiplier, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, file_safety_check, igomega, make_movie, navg, ncheck, nmovie, nsave, nwrite, nwrite_mult, ob_midplane, omegatinst, omegatol, print_flux_line, print_line, save_distfn, save_for_restart, save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, write_bpar_over_time, write_cerr, write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db, write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating, write_hrate, write_jext, write_kpar, write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega, write_parity, write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer

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
    write_gg = self%write_gg
    write_gs = self%write_gs
    write_gyx = self%write_gyx
    write_heating = self%write_heating
    write_hrate = self%write_hrate
    write_jext = self%write_jext
    write_kpar = self%write_kpar
    write_line = self%write_line
    write_lorentzian = self%write_lorentzian
    write_lpoly = self%write_lpoly
    write_max_verr = self%write_max_verr
    write_moments = self%write_moments
    write_nl_flux = self%write_nl_flux
    write_nl_flux_dist = self%write_nl_flux_dist
    write_ntot_over_time = self%write_ntot_over_time
    write_omavg = self%write_omavg
    write_omega = self%write_omega
    write_parity = self%write_parity
    write_pflux_sym = self%write_pflux_sym
    write_pflux_tormom = self%write_pflux_tormom
    write_phi_over_time = self%write_phi_over_time
    write_symmetry = self%write_symmetry
    write_tperp_over_time = self%write_tperp_over_time
    write_upar_over_time = self%write_upar_over_time
    write_verr = self%write_verr
    write_zonal_transfer = self%write_zonal_transfer

    ! Now read in the main namelist
    if (self%get_requires_index()) then
       call get_indexed_namelist_unit(in_file, trim(self%get_name()), self%index)
       read(in_file, nml = gs2_diagnostics_knobs, iostat = iostat)
       exist = .true.
       if (iostat .ne.0) exist = .false.
       close(unit = in_file, iostat = iostat)
    else
       in_file = input_unit_exist(self%get_name(), exist)
       if (exist) read(in_file, nml = gs2_diagnostics_knobs)
    endif

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
    self%write_gg = write_gg
    self%write_gs = write_gs
    self%write_gyx = write_gyx
    self%write_heating = write_heating
    self%write_hrate = write_hrate
    self%write_jext = write_jext
    self%write_kpar = write_kpar
    self%write_line = write_line
    self%write_lorentzian = write_lorentzian
    self%write_lpoly = write_lpoly
    self%write_max_verr = write_max_verr
    self%write_moments = write_moments
    self%write_nl_flux = write_nl_flux
    self%write_nl_flux_dist = write_nl_flux_dist
    self%write_ntot_over_time = write_ntot_over_time
    self%write_omavg = write_omavg
    self%write_omega = write_omega
    self%write_parity = write_parity
    self%write_pflux_sym = write_pflux_sym
    self%write_pflux_tormom = write_pflux_tormom
    self%write_phi_over_time = write_phi_over_time
    self%write_symmetry = write_symmetry
    self%write_tperp_over_time = write_tperp_over_time
    self%write_upar_over_time = write_upar_over_time
    self%write_verr = write_verr
    self%write_zonal_transfer = write_zonal_transfer

    self%exist = exist
  end subroutine read_new_gs2_diagnostics_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_new_gs2_diagnostics_config(self, unit)
    implicit none
    class(new_gs2_diagnostics_config_type), intent(in) :: self
    integer, intent(in) , optional:: unit
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
    call self%write_key_val("write_gg", self%write_gg, unit_internal)
    call self%write_key_val("write_gs", self%write_gs, unit_internal)
    call self%write_key_val("write_gyx", self%write_gyx, unit_internal)
    call self%write_key_val("write_heating", self%write_heating, unit_internal)
    call self%write_key_val("write_hrate", self%write_hrate, unit_internal)
    call self%write_key_val("write_jext", self%write_jext, unit_internal)
    call self%write_key_val("write_kpar", self%write_kpar, unit_internal)
    call self%write_key_val("write_line", self%write_line, unit_internal)
    call self%write_key_val("write_lorentzian", self%write_lorentzian, unit_internal)
    call self%write_key_val("write_lpoly", self%write_lpoly, unit_internal)
    call self%write_key_val("write_max_verr", self%write_max_verr, unit_internal)
    call self%write_key_val("write_moments", self%write_moments, unit_internal)
    call self%write_key_val("write_nl_flux", self%write_nl_flux, unit_internal)
    call self%write_key_val("write_nl_flux_dist", self%write_nl_flux_dist, unit_internal)
    call self%write_key_val("write_ntot_over_time", self%write_ntot_over_time, unit_internal)
    call self%write_key_val("write_omavg", self%write_omavg, unit_internal)
    call self%write_key_val("write_omega", self%write_omega, unit_internal)
    call self%write_key_val("write_parity", self%write_parity, unit_internal)
    call self%write_key_val("write_pflux_sym", self%write_pflux_sym, unit_internal)
    call self%write_key_val("write_pflux_tormom", self%write_pflux_tormom, unit_internal)
    call self%write_key_val("write_phi_over_time", self%write_phi_over_time, unit_internal)
    call self%write_key_val("write_symmetry", self%write_symmetry, unit_internal)
    call self%write_key_val("write_tperp_over_time", self%write_tperp_over_time, unit_internal)
    call self%write_key_val("write_upar_over_time", self%write_upar_over_time, unit_internal)
    call self%write_key_val("write_verr", self%write_verr, unit_internal)
    call self%write_key_val("write_zonal_transfer", self%write_zonal_transfer, unit_internal)
    call self%write_namelist_footer(unit_internal)
  end subroutine write_new_gs2_diagnostics_config

  !> Resets the config object to the initial empty state
  subroutine reset_new_gs2_diagnostics_config(self)
    class(new_gs2_diagnostics_config_type), intent(in out) :: self
    type(new_gs2_diagnostics_config_type) :: empty
    select type (self)
    type is (new_gs2_diagnostics_config_type)
       self = empty
    end select
  end subroutine reset_new_gs2_diagnostics_config

  !> Broadcasts all config parameters so object is populated identically on
  !! all processors
  subroutine broadcast_new_gs2_diagnostics_config(self)
    use mp, only: broadcast
    implicit none
    class(new_gs2_diagnostics_config_type), intent(in out) :: self
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
    call broadcast(self%write_gg)
    call broadcast(self%write_gs)
    call broadcast(self%write_gyx)
    call broadcast(self%write_heating)
    call broadcast(self%write_hrate)
    call broadcast(self%write_jext)
    call broadcast(self%write_kpar)
    call broadcast(self%write_line)
    call broadcast(self%write_lorentzian)
    call broadcast(self%write_lpoly)
    call broadcast(self%write_max_verr)
    call broadcast(self%write_moments)
    call broadcast(self%write_nl_flux)
    call broadcast(self%write_nl_flux_dist)
    call broadcast(self%write_ntot_over_time)
    call broadcast(self%write_omavg)
    call broadcast(self%write_omega)
    call broadcast(self%write_parity)
    call broadcast(self%write_pflux_sym)
    call broadcast(self%write_pflux_tormom)
    call broadcast(self%write_phi_over_time)
    call broadcast(self%write_symmetry)
    call broadcast(self%write_tperp_over_time)
    call broadcast(self%write_upar_over_time)
    call broadcast(self%write_verr)
    call broadcast(self%write_zonal_transfer)

    call broadcast(self%exist)
  end subroutine broadcast_new_gs2_diagnostics_config

  !> Gets the default name for this namelist
  function get_default_name_new_gs2_diagnostics_config()
    implicit none
    character(len = CONFIG_MAX_NAME_LEN) :: get_default_name_new_gs2_diagnostics_config
    get_default_name_new_gs2_diagnostics_config = "gs2_diagnostics_knobs"
  end function get_default_name_new_gs2_diagnostics_config

  !> Gets the default requires index for this namelist
  function get_default_requires_index_new_gs2_diagnostics_config()
    implicit none
    logical :: get_default_requires_index_new_gs2_diagnostics_config
    get_default_requires_index_new_gs2_diagnostics_config = .false.
  end function get_default_requires_index_new_gs2_diagnostics_config
end module diagnostics_config

