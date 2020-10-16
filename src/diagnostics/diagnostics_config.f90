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
     !>
     !> @warning Only present in new diagnostics
     logical :: append_old = .false.
     !> The maximum number of steps allowed before declaring convergence with a warning
     integer :: conv_max_step = 80000
     !> The minimum number of steps before the convergence condition will be met
     integer :: conv_min_step = 4000
     !> The number of timesteps the convergence condition averages over
     integer :: conv_nstep_av = 4000
     !> The number of steps where convergence is true before convergence is accepted
     integer :: conv_nsteps_converged = 10000
     !> Multiplier for the cumulative average of the heat flux when testing
     !> convergence with [[gs2_diagnostics_knobs:use_nonlin_convergence]]
     real :: conv_test_multiplier = 4e-1
     !> Write out the field-line average of \(\phi\) to `dump.check1`. This option
     !> is usually used for Rosenbluth-Hinton calculations.
     !>
     !> @warning You probably don't want this?
     !>
     !> @warning Non-functional in new diagnostics
     logical :: dump_check1 = .false.
     !> Write out \(A_\parallel(k_x, k_y)\) at `igomega` to `<run name>.dc2`.
     !>
     !> @warning You probably don't want this?
     !>
     !> @warning Non-functional in new diagnostics
     logical :: dump_check2 = .false.
     !> Write out \(\phi, A_\parallel, B_\parallel\) to `dump.fields.t=(time)`
     !> every `10 * nwrite` timesteps.
     !>
     !> @warning You probably don't want this?  Expensive.
     logical :: dump_fields_periodically = .false.
     !> If built with parallel IO capability enable parallel IO. Currently
     !> disabled by default because of problems on Helios: parallel IO seems to
     !> work on 1 node but not on >1. This problem does not affect parallel IO
     !> for the single restart file because it only affects output with infinite
     !> dimensions
     !>
     !> @warning Only in new diagnostics
     logical :: enable_parallel = .false.
     !> Exit when the run has reached convergence
     logical :: exit_when_converged = .true.
     !> Verify that the restart file(s) can be written before starting the
     !> run.
     !>
     !> @warning Currently non-functional in new diagnostics
     logical :: file_safety_check = .true.
     !> Index in \(\theta\) to measure various diagnostics at
     integer :: igomega = 0
     !> Write some fields out in real space over time, suitable for creating
     !> animations. Timestep period is controlled with
     !> [[gs2_diagnostics_knobs:nmovie]]
     !>
     !>
     !> @warning This option can write a lot of data! Consider doing the Fourier
     !> transforms of the fields in post-processing instead.
     !>
     !> @warning In new diagnostics, the timestep period is controlled with
     !> [[gs2_diagnostics_knobs:nwrite]]
     logical :: make_movie = .false.
     !> Number of timesteps to average over for time-averaged diagnostics
     integer :: navg = 10
     !> FIXME: Timestep period of which to check velocity space resolution and
     !> correction by varying collisionality. But doesn't happen on timesteps
     !> when diagnostics are written.
     !>
     !> @warning This value "shadows" the `ncheck` input of the
     !> collisions namelist.
     integer :: ncheck = 10 ! Only in new diagnostics
     !> Sets the number of output steps between syncing the netcdf file.
     !> Higher values may reduce disk overheads but increase the risk of
     !> data loss on unclean exit.
     !>
     !> @warning Not present in new diagnostics
     ! integer :: nc_sync_freq = 10

     !> Timestep period to write real space fields \(\phi, A_\parallel,
     !> B_\parallel\)
     !>
     !> @warning Non-functional in new diagnostics. Instead use
     !> [[gs2_diagnostics_knobs:nwrite]]
     integer :: nmovie = -1
     !> Timestep period for writing restart data if
     !> [[gs2_diagnostics_knobs:save_for_restart]] is `.true.`
     integer :: nsave = -1
     !> Timestep period for writing outputs
     !>
     !> @warning Default value differs between old (`== 100`) and new (`== 10`)
     !> diagnostics!
     integer :: nwrite = 10
     !> Timestep period multiplier for certain "large" diagnostics, which are
     !> written every `nwrite_mult * nwrite` timesteps.
     !>
     !> FIXME: What datasets? `phicorr_sum`, `phiextend_sum` in old diagnostics,
     !> `f` and `fyx` in new diagnostics?
     !>
     !> @warning Controls different diagnostics in new diagnostics
     integer :: nwrite_mult = 10
     !> If `.true.`, write moments at the outboard midplane only, as opposed to
     !> along the entire flux surface
     !>
     !> @warning Non-functional in new diagnostics
     !>
     !> @warning Default value differs between old (.true.) and new (.false.) diagnostics
     logical :: ob_midplane = .false.
     !> Threshold frequency (\(\omega\)) for detecting a numerical
     !> instability. If \(\omega\) averaged over [[gs2_diagnostics_knobs:navg]]
     !> timesteps is greater than `omegatinst`, abort the run
     !>
     !> @warning Default value differs between old (1.0) and new (1.0e6)
     !> diagnostics
     real :: omegatinst = 1.0e6
     !> Frequency (\(\omega\)) convergence tolerance. Consider the simulation
     !> converged and finish the run if \(\omega\) has changed by less than
     !> `omegatol` over the last [[gs2_diagnostics_knobs:navg]] timesteps.
     !>
     !> More explicitly, the convergence criterion is:
     !>
     !> $$\sqrt{\sum_t^{n_{avg}} |\bar{\omega} - \omega|^2} < \min{\bar{\omega}, 1.0}\cdot\omega_{tol}$$
     !>
     !> where \(\bar{\omega}\) is \(\omega\) averaged over the last `navg`
     !> timesteps.
     real :: omegatol = 1e-3
     !> Print the instantaneous fluxes to screen/stdout every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: print_flux_line = .false.
     !> Print estimated frequencies and growth rates to screen/stdout every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: print_line = .false.
     !> Write the distribution function to `<rootname>.nc.dfn.<processor>`
     !> files. Only written at end of simulation
     logical :: save_distfn = .false.
     !> Write restart files every [[gs2_diagnostics_knobs:nwrite]] timesteps and
     !> at end of simulation
     logical :: save_for_restart = .false.
     !> Save some layout and distribution information in restart files
     !>
     !> @warning Only in new diagnostics
     logical :: save_glo_info_and_grids = .false.
     !> If `.true.`, write multiple restart files. Requires `gs2` to have been
     !> built with parallel netCDF
     logical :: save_many = .false.
     !> Save parallel and perpendicular velocities in final restart and/or
     !> distribution function files
     !>
     !> @warning Only in new diagnostics
     logical :: save_velocities = .false.
     !> @warning Not used
     !>
     !> @warning Only in new diagnostics
     logical :: serial_netcdf4 = .false.
     !> For nonlinear runs, stop the run when the averaged differential of the
     !> summed averaged heat flux drops below a threshold, controlled by
     !> [[gs2_diagnostics_knobs:conv_test_multiplier]]
     logical :: use_nonlin_convergence = .false.
     !> If `.false.`, skip writing any diagnostics
     !>
     !> @warning Only in new diagnostics
     logical :: write_any = .true.
     !> Write the entire \(A_\parallel\) field every
     !> [gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning New diagnostics requires
     !> [[gs2_diagnostics_knobs:write_fields]] to also be `.true.`
     !> (the default) to enable this
     logical :: write_apar_over_time = .false.
     !> Write diagnostics to text files. Generally this creates a different text
     !> file for each diagnostic. Note that this independent of whether or not
     !> netCDF files are used.
     !>
     !> FIXME: Verify old and new diagnostics write these files in same format
     logical :: write_ascii = .true.
     !> Write flux surface averaged low-order moments of \(g\) every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Prints a warning that this is ignored unless
     !> `grid_option==box`, but this doesn't appear to be the case?
     logical :: write_avg_moments = .false.
     !> Write the entire \(B_\parallel\) field every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning New diagnostics requires
     !> [[gs2_diagnostics_knobs:write_fields]] to also be `.true.`
     !> (the default) to enable this
     logical :: write_bpar_over_time = .false.
     !> Write the collision error every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps
     !>
     !> FIXME: What does this mean?
     !>
     !> @warning In new diagnostics, only writes to ascii
     logical :: write_cerr = .false.
     !> Write collisional heating (collisional and hyper viscous rate of loss of
     !> free energy for each mode) every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps
     !>
     !> @warning Only in new diagnostics
     logical :: write_collisional = .false.
     !> Write two point parallel correlation function calculated from the
     !> electric potential as a function of \(k_y\) and parallel separation
     !> every [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: write_correlation = .true.
     !> Write two point parallel correlation function calculated from the
     !> electric potential as a function of \(k_y\) and parallel separation,
     !> time averaged and calculated along the extended domain every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: write_correlation_extend = .false.
     !> Write the cross phase between electron density and perpendicular
     !> electron temperature every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps. Calculated at both the outboard midplane, averaged across
     !> \(x\) and \(y\); and averaged across all space
     !>
     !> @warning In old diagnostics, only written to text files
     logical :: write_cross_phase = .false.
     !> Write the whole non-adiabatic part of the density moment every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Requires [[gs2_diagnostics_knobs:write_moments]] to also be
     !> `.true.`
     !>
     !> @warning Only in new diagnostics
     logical :: write_density_over_time = .false.
     !> Write \(\phi, A_\parallel, B_\parallel\) normalised to the value of
     !> \(\phi\) at the outboard midplane every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps.
     !>
     !> If [[gs2_diagnostics_knobs:write_ascii]] is enabled, the text file is
     !> `<runname>.eigenfunc`.
     !>
     !> @note The normalising factor for a given \((\theta, k_y)\) point may not
     !> be exactly \(\phi\) at the outboard midplane if this value is zero
     !> there. If the adjacent point in \(\theta\) is also zero, then the fields
     !> will be unnormalised at that \((\theta, k_y)\).
     !>
     !> @warning If this option is turned on, the same normalising factor will
     !> also be used for the text output in
     !> [[gs2_diagnostics_knobs:write_final_moments]]
     !>
     !> @warning In old diagnostics, this is only done on exit, not every
     !> `nwrite` timesteps!
     !>
     !> @warning The text output is normalised (see note above), but the netCDF
     !> output is not
     logical :: write_eigenfunc = .false.
     !> Enable writing \(\phi, A_\parallel, B_\parallel\) every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning In old diagnostics, this only writes them on exit and
     !> [[gs2_diagnostics_knobs:write_phi_over_time]],
     !> [[gs2_diagnostics_knobs:write_apar_over_time]] and
     !> [[gs2_diagnostics_knobs:write_bpar_over_time]] instead control writing
     !> the fields over time, respectively. In new diagnostics, this will disable
     !> writing all three fields if set to `.false.`
     !>
     !> @warning Default value differs between old (`.false.`) and new
     !> (`.true.`) diagnostics
     logical :: write_fields = .true.
     !> Write the right-hand sides of the field equations at the final
     !> timestep. If [[gs2_diagnostics_knobs:write_ascii]] is enabled, the file
     !> suffix is `.antot`
     !>
     !> @warning In new diagnostics, this is not written to netCDF!
     logical :: write_final_antot = .false.
     !> Write \(\delta B\) at the final timestep. If
     !> [[gs2_diagnostics_knobs:write_ascii]] is enabled, the file suffix is
     !> `.db`
     !>
     !> @warning In new diagnostics, this is not written to netCDF!
     logical :: write_final_db = .false.
     !> Write \(E_\parallel\) at the final timestep. If
     !> [[gs2_diagnostics_knobs:write_ascii]] is enabled, the file suffix is
     !> `.epar`
     !>
     !> @warning In new diagnostics, this is not written to netCDF!
     logical :: write_final_epar = .false.
     !> Write \(\phi, A_\parallel, B_\parallel\) at the final timestep. If
     !> [[gs2_diagnostics_knobs:write_ascii]] is enabled, the file suffix is
     !> `.fields`
     !>
     !> @warning In new diagnostics, this is not written to netCDF!
     logical :: write_final_fields = .false.
     !> Write various moments (densities, parallel and perpendicular
     !> velocities and temperatures, and heat and momentum fluxes) at
     !> the final timestep. If [[gs2_diagnostics_knobs:write_ascii]]
     !> is enabled, the file suffix is `.fields` and contains the
     !> moments, their magnitudes, and field-line averages. The netCDF
     !> output has only the values.
     !>
     !> @warning If [[gs2_diagnostics_knobs:write_eigenfunc]] is
     !> enabled, then the same normalising factor is used for the
     !> final moments. Otherwise, the moments are unnormalised.
     !>
     !> @warning In new diagnostics, this is not written to netCDF!
     !>
     !> @warning The text output is normalised (see note above), but
     !> the netCDF output is not
     logical :: write_final_moments = .false.
     !> Write instantaneous fluxes to output file every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Output formats and quantities differ between old and new
     !> diagnostics
     logical :: write_flux_line = .true.
     !> Write fluxes (heat, momentum and particle; total and per-species) every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps.
     !>
     !> If run is nonlinear, this defaults to true
     !>
     !> @warning Only in new diagnostics. In old diagnostics, this is combined
     !> with write_fluxes_by_mode and called write_nl_flux.
     logical :: write_fluxes = .false.
     !> Write fluxes (heat, momentum and particle; total and per-species) as a
     !> function of Fourier mode every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps.
     !>
     !> If run is nonlinear, this defaults to true
     !>
     !> @warning Only in new diagnostics. In old diagnostics, this is combined
     !> with write_fluxes and called write_nl_flux.
     logical :: write_fluxes_by_mode = .false.
     !> Write moments (density, parallel flow, and parallel and perpendicular
     !> temperatures) in non-guiding centre coordinates every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: write_full_moments_notgc = .false.
     !> FIXME: Write either \(g\) or \(f\)? every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Only to text files
     logical :: write_g = .false.
     !> @warning Unused
     logical :: write_gg = .false.
     !> FIXME: Add documentation
     !>
     !> @warning Only outputs to text files
     logical :: write_gs = .false.
     !> FIXME: Add documentation
     !>
     !> @warning Only outputs to text files
     logical :: write_gyx = .false.
     !> Write out various heating, free energy and energy injection
     !> diagnostics. Text file extension is `.heat`
     !>
     !> @warning Only in new diagnostics. Equivalent to
     !> [[gs2_diagnostics_knobs:write_hrate]] in old diagnostics
     logical :: write_heating = .false.
     !> Write out various heating, free energy and energy injection
     !> diagnostics. Text file extension is `.heat`
     !>
     !> @warning Non-functional in new diagnostics. Equivalent to
     !> [[gs2_diagnostics_knobs:write_heating]] in new diagnostics
     !>
     !> @warning Old diagnostics only writes to text files
     logical :: write_hrate = .false.
     !> Write time averaged external current in the antenna,
     !> \(\operatorname{Re}(k_\perp^2 A_\mathrm{antenna})\), as a function of
     !> \(k_x, k_y\). File suffix is `.jext`
     !>
     !> @warning Only in new diagnostics
     logical :: write_jext = .false.
     !> Write the parallel spectrum of \(\phi, A_\parallel,
     !> B_\parallel\) at final timestep
     !>
     !> @warning Output only to text files
     logical :: write_kpar = .false.
     !> Print estimated frequencies and growth rates to text file every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     logical :: write_line = .true.
     !> Write `antenna_w` every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps
     !>
     !> FIXME: Define `antenna_w`
     !>
     !> @warning Old diagnostics only writes to text files
     !>
     !> @warning New diagnostics only writes to netCDF
     logical :: write_lorentzian = .false.
     !> Write Lagrange interpolating function value at interpolating
     !> points specified by `xloc` to text file with suffix `.interp`
     !>
     !> FIXME: Define `xloc`
     logical :: write_lpoly = .false.
     !> Write the estimated maximum error from velocity space
     !> integrals for various quantities
     !>
     !> @warning Old diagnostics only outputs to text file
     logical :: write_max_verr = .false.
     !> Write various moments (total and non-adiabatic part of perturbed species
     !> density, perturbed parallel flow, perturbed parallel and perpendicular
     !> temperatures, parallel heat flux, particule flux and heat flux) every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Only outputs to netCDF files
     !>
     !> @warning New diagnostics writes various averages of the moments. To get
     !> the same output as old diagnostics, **also** set
     !> [[gs2_diagnostics_knobs:write_ntot_over_time]],
     !> [[gs2_diagnostics_knobs:write_density_over_time]],
     !> [[gs2_diagnostics_knobs:write_upar_over_time]],
     !> [[gs2_diagnostics_knobs:write_tperp_over_time]] to `.true.`
     logical :: write_moments = .true.
     !> Write various fluxes every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps
     !>
     !> @warning Non-functional in new diagnostics
     logical :: write_nl_flux = .false.
     !> Write the poloidal distributions of the fluxes of particles,
     !> parallel momentum, perpendicular momentum, and energy every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps.
     !>
     !> See section 3.1 and appendix A of Ball et al. PPCF 58 (2016) 045023 as well
     !> as section 5 of "GS2 analytic geometry specification"
     logical :: write_nl_flux_dist = .false.
     !> Write the whole total density moment every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Requires [[gs2_diagnostics_knobs:write_moments]] to also be
     !> `.true.`
     !>
     !> @warning Only in new diagnostics
     logical :: write_ntot_over_time = .false.
     !> Write time-averaged growth rate and frequency to the output text file
     !> every [[gs2_diagnostics_knobs:nwrite]] timesteps. Time average is
     !> rolling window over the previous [[gs2_diagnostics_knobs:navg]]
     !> timesteps
     !>
     !> @warning Non-functional in new diagnostics, time-average is also written
     !> with [[gs2_diagnostics_knobs:write_omega]] instead
     logical :: write_omavg = .true.
     !> Write instantaneous growth rate and frequency every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Old diagnostics only writes to output text file
     !>
     !> @warning New diagnostics only writes to netCDF file
     !>
     !> @warning New diagnostics also writes time-average with this
     !> option. Output is the complex frequency, rather than separate growth
     !> rate and real frequency
     logical :: write_omega = .true.
     !> Write parities in distribution and particle fluxes to text file with the
     !> suffix `.parity`
     !>
     !> FIXME: Clarify what this means
     logical :: write_parity = .false.
     !> Write particle flux density as a function of \(\theta\) and velocity
     !> space every [[gs2_diagnostics_knobs:nwrite]] timesteps. Used for looking
     !> at the effect of asymmetry. See Parra et al POP 18 062501 2011 and ask
     !> Jung-Pyo Lee
     !>
     !> @warning Only outputs to netCDF
     !>
     !> @warning Non-functional in new diagnostics, use
     !> [[gs2_diagnostics_knobs:write_symmetry]] instead
     logical :: write_pflux_sym = .false. 
     !> Write toroidal angular momentum flux carried by particle flux every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps. Only calculated if `gs2` is
     !> built with `LOWFLOW` on.
     !>
     !> @warning Only outputs to netCDF
     !>
     !> @warning Non-functional in new diagnostics
     logical :: write_pflux_tormom = .false. 
     !> Write the entire \(\phi\) field every [gs2_diagnostics_knobs:nwrite]]
     !> timesteps
     !>
     !> @warning New diagnostics requires
     !> [[gs2_diagnostics_knobs:write_fields]] to also be `.true.`
     !> (the default) to enable this
     logical :: write_phi_over_time = .false.
     !> Write the particle and momentum flux as a function of \(\theta\) and
     !> velocity space. See [[gs2_diagnostics_knobs:write_pflux_sym]] and
     !> [[gs2_diagnostics_knobs:write_pflux_tormom]]
     !>
     !> @warning Only outputs to netCDF
     !>
     !> @warning Old diagnostics does not write the particle flux, use
     !> [[gs2_diagnostics_knobs:write_pflux_sym]] instead
     logical :: write_symmetry = .false.
     !> Write the whole perturbed perpendicular temperature moment every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Requires [[gs2_diagnostics_knobs:write_moments]] to also be
     !> `.true.`
     !>
     !> @warning Only in new diagnostics
     logical :: write_tperp_over_time = .false.
     !> Write the whole perturbed parallel velocity moment every
     !> [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Requires [[gs2_diagnostics_knobs:write_moments]] to also be
     !> `.true.`
     !>
     !> @warning Only in new diagnostics
     logical :: write_upar_over_time = .false.
     !> Write estimates of error resulting from velocity space integrals in the
     !> calculation of various quantities every [[gs2_diagnostics_knobs:nwrite]]
     !> timesteps.
     !>
     !> @warning Output only to text file
     !>
     !> @warning New diagnostics also writes to netCDF every
     !> [[gs2_diagnostics_knobs:ncheck]] timesteps, along with some other
     !> quantities such as collisionality
     logical :: write_verr = .false.
     !> Write the transfer of free energy, \(\tau\), as a function of \((k_x,
     !> k_y)\) every [[gs2_diagnostics_knobs:nwrite]] timesteps
     !>
     !> @warning Only in new diagnostics
     logical :: write_zonal_transfer = .false.
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
    integer :: conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, igomega, navg, ncheck, nmovie, nsave, nwrite, nwrite_mult
    logical :: append_old, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, file_safety_check, make_movie, ob_midplane, print_flux_line, print_line, save_distfn
    logical :: save_for_restart, save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, write_bpar_over_time
    logical :: write_cerr, write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db
    logical :: write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating
    logical :: write_hrate, write_jext, write_kpar, write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega
    logical :: write_parity, write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer
    real :: conv_test_multiplier, omegatinst, omegatol

    namelist /gs2_diagnostics_knobs/ append_old, conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, conv_test_multiplier, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, &
         file_safety_check, igomega, make_movie, navg, ncheck, nmovie, nsave, nwrite, nwrite_mult, ob_midplane, omegatinst, omegatol, print_flux_line, print_line, save_distfn, save_for_restart, &
         save_glo_info_and_grids, save_many, save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, write_avg_moments, write_bpar_over_time, write_cerr, &
         write_collisional, write_correlation, write_correlation_extend, write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db, write_final_epar, &
         write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating, write_hrate, &
         write_jext, write_kpar, write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist, write_ntot_over_time, write_omavg, write_omega, write_parity, &
         write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry, write_tperp_over_time, write_upar_over_time, write_verr, write_zonal_transfer

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

