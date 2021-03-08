module diagnostics_base_config
  use abstract_config, only: abstract_config_type, CONFIG_MAX_NAME_LEN
  implicit none
  !> Used to represent the input configuration of
  !> gs2_diagnostics. This is the version used by the original
  !> diagnostics module.
  type, extends(abstract_config_type) :: diagnostics_base_config_type
    ! namelist : gs2_diagnostics_knobs
    ! indexed : false

    !> If true, and a copy of the .out.nc already exists,
    !> open it and append to it. Obviously grid sizes must be
    !> unchanged.
    !>
    !> @warning Only present in new diagnostics
    logical :: append_old = .false.

    !> Used for Trinity nonlinear convergence
    !> ([[gs2_diagnostics_knobs:use_nonlin_convergence]]). The maximum number
    !> of Trinity steps, after which consider the run converged regardless
    integer :: conv_max_step = 80000
    !> Used for the Trinity nonlinear convergence check
    !> ([[gs2_diagnostics_knobs:use_nonlin_convergence]]). The minimum number
    !> of Trinity steps before checking convergence condition
    integer :: conv_min_step = 4000
    !> Used for the Trinity nonlinear convergence check
    !> ([[gs2_diagnostics_knobs:use_nonlin_convergence]]). The number of
    !> timesteps the convergence condition averages over
    integer :: conv_nstep_av = 4000
    !> Used for the Trinity nonlinear convergence check
    !> ([[gs2_diagnostics_knobs:use_nonlin_convergence]]). The number of steps
    !> where convergence is true before convergence is accepted
    integer :: conv_nsteps_converged = 10000
    !> Used for the Trinity nonlinear convergence check
    !> ([[gs2_diagnostics_knobs:use_nonlin_convergence]]). Multiplier for the
    !> cumulative average of the heat flux
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
    !> Index in \(\theta\) to measure various diagnostics at. By default, this
    !> is the outboard midplane \(\theta = 0\), but it may be useful to set
    !> this to a non-zero value for tearing parity modes.
    integer :: igomega = 0
    !> Write out \(\phi, A_\parallel, B_\parallel\) in real space over time,
    !> suitable for creating animations. Timestep period is controlled with
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
    !>
    !> @warning Default value differs between old (100) and new (10)
    !> diagnostics
    integer :: navg = 100

    !> FIXME: Timestep period of which to check velocity space resolution and
    !> correction by varying collisionality. But doesn't happen on timesteps
    !> when diagnostics are written.
    !>
    !> @warning Only in new diagnostics
    !>
    !> @warning This value "shadows" the `ncheck` input of the
    !> collisions namelist.
    integer :: ncheck = 10

    !> Sets the number of output steps between syncing the netcdf file.
    !> Higher values may reduce disk overheads but increase the risk of
    !> data loss on unclean exit.
    !>
    !> @warning Not present in new diagnostics
    integer :: nc_sync_freq = 10
    !> Timestep period to write real space fields \(\phi, A_\parallel,
    !> B_\parallel\)
    !>
    !> @warning Non-functional in new diagnostics. Instead use
    !> [[gs2_diagnostics_knobs:nwrite]]
    !>
    !> @warning Default value differs between old (1000) and new (-1)
    !> diagnostics
    integer :: nmovie = 1000
    !> Timestep period for writing restart data if
    !> [[gs2_diagnostics_knobs:save_for_restart]] is `.true.`. Negative values
    !> disable the periodic checkpoints.
    integer :: nsave = -1
    !> Timestep period for writing outputs
    !>
    !> @warning Default value differs between old (`== 100`) and new (`== 10`)
    !> diagnostics!
    integer :: nwrite = 100
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
    logical :: ob_midplane = .true.
    !> Threshold complex frequency (\(\Omega\)) for detecting a numerical
    !> instability. If \(\abs(\Omega)\) averaged over
    !> [[gs2_diagnostics_knobs:navg]] timesteps is greater than `omegatinst`,
    !> abort the run.
    !>
    !> @warning Default value differs between old (1.0) and new (1.0e6)
    !> diagnostics
    real :: omegatinst = 1.0
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
    !>
    !> @warning Default value differs between old (`.true.`) and new
    !> (`.false.`) diagnostics
    logical :: print_line = .true.
    !> Write the distribution function to `<rootname>.nc.dfn.<processor>`
    !> files. Only written at end of simulation
    logical :: save_distfn = .false.
    !> Write restart files at the end of the simulation. If
    !> [[gs2_diagnostics_knobs:nsave]] is positive, then also enable
    !> checkpoints by writing restart files every `nsave` timesteps.
    logical :: save_for_restart = .false.

    !> Save some layout and distribution information in restart files
    !>
    !> @warning Only in new diagnostics
    logical :: save_glo_info_and_grids = .false.

    !> If `.true.`, write one restart file per processor, otherwise write a
    !> single restart file.
    !>
    !> If `gs2` has not been built with parallel netCDF, `save_many` is ignored
    !> and there is always one file per processor (equivalent to `save_many =
    !> .true.`).
    !>
    !> If `write_many` is enabled, you probably want to also set
    !> [[init_g_knobs:read_many]] in order to restart from multiple files.
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
    !> summed averaged heat flux drops below a threshold for a period of time,
    !> controlled by [[gs2_diagnostics_knobs:conv_test_multiplier]],
    !> [[gs2_diagnostics_knobs:conv_nsteps_converged]],
    !> [[gs2_diagnostics_knobs:conv_nsteps_av]],
    !> [[gs2_diagnostics_knobs:conv_max_step]], and
    !> [[gs2_diagnostics_knobs:conv_min_step]].
    logical :: use_nonlin_convergence = .false.

    !> If `.false.`, skip writing any diagnostics
    !>
    !> @warning Only in new diagnostics
    !>
    !> @warning This also turns off checking the linear convergence
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
    !> timesteps to text file with suffix `.cres`
    !>
    !> FIXME: What does this mean?
    !>
    !> @warning In new diagnostics, only writes to ascii
    !>
    !> @warning This is expensive
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
    !>
    !> @warning Default value differs between old (`.false.`) and new
    !> (`.true.`) diagnostics
    logical :: write_correlation = .false.
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
    logical :: write_fields = .false.
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
    !> Write \(g(v_\parallel,v_\perp)\) at `ik=it=is=1, ig=0` to text file
    !> `<runname>.dist`
    !>
    !> @warning Only to text files
    logical :: write_g = .false.
    !> @warning Unused
    !>
    !> FIXME: Delete
    logical :: write_gg = .false.
    !> FIXME: Add documentation
    !>
    !> @warning Only outputs to text files
    logical :: write_gs = .false.
    !> Write \(g\) as a function of real space every
    !> [[gs2_diagnostics_knobs:nmovie]] timesteps to text file
    !> "<runname>.yxdist"
    !>
    !> @warning Only to text files
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
    !> B_\parallel\) at final timestep. File suffix is `.kpar`
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
    !>
    !> @warning Default value differs between old (`.false.`) and new
    !> (`.true.`) diagnostics
    logical :: write_moments = .false.
    !> Write various fluxes every [[gs2_diagnostics_knobs:nwrite]]
    !> timesteps
    !>
    !> @warning Non-functional in new diagnostics
    logical :: write_nl_flux = .false.
    !> Write the poloidal distributions of the fluxes of particles,
    !> parallel momentum, perpendicular momentum, and energy every
    !> [[gs2_diagnostics_knobs:nwrite]] timesteps.
    !>
    !> See section 3.1 and appendix A of [Ball et al. PPCF 58 (2016)
    !> 045023](https://doi.org/10.1088/0741-3335/58/4/045023) as well as
    !> section 5 of "GS2 analytic geometry specification"
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
    !> at the effect of asymmetry. See [Parra et al POP 18 062501
    !> 2011](https://doi.org/10.1063/1.3586332) and ask Jung-Pyo Lee
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
    !>
    !> @warning This is expensive
    logical :: write_verr = .false.

    !> Write the transfer of free energy, \(\tau\), as a function of \((k_x,
    !> k_y)\), averaged over \(\theta\), every [[gs2_diagnostics_knobs:nwrite]]
    !> timesteps
    !>
    !> @warning Only in new diagnostics
    logical :: write_zonal_transfer = .false.
  contains
    procedure, public :: read => read_diagnostics_base_config
    procedure, public :: write => write_diagnostics_base_config
    procedure, public :: reset => reset_diagnostics_base_config
    procedure, public :: broadcast => broadcast_diagnostics_base_config
    procedure, public, nopass :: get_default_name => get_default_name_diagnostics_base_config
    procedure, public, nopass :: get_default_requires_index => get_default_requires_index_diagnostics_base_config
  end type diagnostics_base_config_type

contains

  !> Prints a warning about a non-functional diagnostic
  subroutine warn_about_nonfunctional_selection(diagnostic_switch, name)
    use run_parameters, only: use_old_diagnostics
    !> Diagnostic switch
    logical, intent(in) :: diagnostic_switch
    !> Name of the switch
    character(len=*), intent(in) :: name

    character(len=3) :: diagnostics_module

    if (use_old_diagnostics) then
      diagnostics_module = "old"
    else
      diagnostics_module = "new"
    end if

    if (diagnostic_switch) then
      write(*, '(a, a, a, a, a)') "WARNING: diagnostic '", name, &
           "' is (currently) non-functional in ", diagnostics_module, &
           " diagnostics module. This will have no effect."
    end if
  end subroutine warn_about_nonfunctional_selection

  !> Reads in the gs2_diagnostics_knobs namelist and populates the member variables
  subroutine read_diagnostics_base_config(self)
    use file_utils, only: input_unit_exist, get_indexed_namelist_unit
    use mp, only: proc0
    implicit none
    class(diagnostics_base_config_type), intent(in out) :: self
    logical :: exist
    integer :: in_file, iostat

    ! Note: When this routine is in the module where these variables live
    ! we shadow the module level variables here. This is intentional to provide
    ! isolation and ensure we can move this routine to another module easily.
    integer :: conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, igomega, navg, nc_sync_freq, ncheck, nmovie, nsave
    integer :: nwrite, nwrite_mult
    logical :: append_old, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, file_safety_check
    logical :: make_movie, ob_midplane, print_flux_line, print_line, save_distfn, save_for_restart, save_glo_info_and_grids, save_many
    logical :: save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii
    logical :: write_avg_moments, write_bpar_over_time, write_cerr, write_collisional, write_correlation, write_correlation_extend
    logical :: write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db
    logical :: write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode
    logical :: write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating, write_hrate, write_jext, write_kpar
    logical :: write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist
    logical :: write_omavg, write_omega, write_parity, write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry
    logical :: write_upar_over_time, write_verr, write_zonal_transfer
    real :: conv_test_multiplier, omegatinst, omegatol

    namelist /gs2_diagnostics_knobs/ append_old, conv_max_step, conv_min_step, conv_nstep_av, conv_nsteps_converged, &
         conv_test_multiplier, dump_check1, dump_check2, dump_fields_periodically, enable_parallel, exit_when_converged, &
         file_safety_check, igomega, make_movie, navg, nc_sync_freq, ncheck, nmovie, nsave, nwrite, nwrite_mult, ob_midplane, &
         omegatinst, omegatol, print_flux_line, print_line, save_distfn, save_for_restart, save_glo_info_and_grids, save_many, &
         save_velocities, serial_netcdf4, use_nonlin_convergence, write_any, write_apar_over_time, write_ascii, &
         write_avg_moments, write_bpar_over_time, write_cerr, write_collisional, write_correlation, write_correlation_extend, &
         write_cross_phase, write_density_over_time, write_eigenfunc, write_fields, write_final_antot, write_final_db, &
         write_final_epar, write_final_fields, write_final_moments, write_flux_line, write_fluxes, write_fluxes_by_mode, &
         write_full_moments_notgc, write_g, write_gg, write_gs, write_gyx, write_heating, write_hrate, write_jext, write_kpar, &
         write_line, write_lorentzian, write_lpoly, write_max_verr, write_moments, write_nl_flux, write_nl_flux_dist, &
         write_omavg, write_omega, write_parity, write_pflux_sym, write_pflux_tormom, write_phi_over_time, write_symmetry, &
         write_upar_over_time, write_verr, write_zonal_transfer

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
    write_omavg = self%write_omavg
    write_omega = self%write_omega
    write_parity = self%write_parity
    write_pflux_sym = self%write_pflux_sym
    write_pflux_tormom = self%write_pflux_tormom
    write_phi_over_time = self%write_phi_over_time
    write_symmetry = self%write_symmetry
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
    self%write_omavg = write_omavg
    self%write_omega = write_omega
    self%write_parity = write_parity
    self%write_pflux_sym = write_pflux_sym
    self%write_pflux_tormom = write_pflux_tormom
    self%write_phi_over_time = write_phi_over_time
    self%write_symmetry = write_symmetry
    self%write_upar_over_time = write_upar_over_time
    self%write_verr = write_verr
    self%write_zonal_transfer = write_zonal_transfer

    self%exist = exist
  end subroutine read_diagnostics_base_config

  !> Writes out a namelist representing the current state of the config object
  subroutine write_diagnostics_base_config(self, unit)
    implicit none
    class(diagnostics_base_config_type), intent(in) :: self
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
    call self%write_key_val("write_omavg", self%write_omavg, unit_internal)
    call self%write_key_val("write_omega", self%write_omega, unit_internal)
    call self%write_key_val("write_parity", self%write_parity, unit_internal)
    call self%write_key_val("write_pflux_sym", self%write_pflux_sym, unit_internal)
    call self%write_key_val("write_pflux_tormom", self%write_pflux_tormom, unit_internal)
    call self%write_key_val("write_phi_over_time", self%write_phi_over_time, unit_internal)
    call self%write_key_val("write_symmetry", self%write_symmetry, unit_internal)
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
    call broadcast(self%write_omavg)
    call broadcast(self%write_omega)
    call broadcast(self%write_parity)
    call broadcast(self%write_pflux_sym)
    call broadcast(self%write_pflux_tormom)
    call broadcast(self%write_phi_over_time)
    call broadcast(self%write_symmetry)
    call broadcast(self%write_upar_over_time)
    call broadcast(self%write_verr)
    call broadcast(self%write_zonal_transfer)

    call broadcast(self%exist)
  end subroutine broadcast_diagnostics_base_config

  !> Gets the default name for this namelist
  function get_default_name_diagnostics_base_config()
    use abstract_config, only : CONFIG_MAX_NAME_LEN
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

end module diagnostics_base_config
