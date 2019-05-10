title: gs2_diagnostics_knobs

#gs2_diagnostics_knobs

| Name | Type | Default | Description |
|:-----|:----:|:-------:|:------------|
| append_old | logical | .false.  | If true, and a copy of the .out.nc already exists,|
| | | |open it and append to it. Obviously grid sizes must be|
| | | |unchanged. |
| conv_max_step | integer | 80000 | FIXME: Add documentation |
| conv_min_step | integer | 4000 | FIXME: Add documentation |
| conv_nstep_av | integer | 4000 | FIXME: Add documentation |
| conv_nsteps_converged | integer | 10000 | FIXME: Add documentation |
| conv_test_multiplier | real | 4e-1 | FIXME: Add documentation |
| dump_check1 | logical | .false. | FIXME: Add documentation |
| dump_check2 | logical | .false. | FIXME: Add documentation |
| dump_fields_periodically | logical | .false. | FIXME: Add documentation |
| enable_parallel | logical | .false.  | If built with parallel IO capability|
| | | |enable parallel IO. Currently disabled by|
| | | |default because of problems on Helios:|
| | | |parallel IO seems to work on 1 node but not|
| | | |on >1. This problem does not affect parallel|
| | | |IO for the single restart file because it only|
| | | |affects output with infinite dimensions |
| exit_when_converged | logical | .true. | FIXME: Add documentation |
| file_safety_check | logical | .true. | FIXME: Add documentation |
| igomega | integer | 0 | FIXME: Add documentation |
| make_movie | logical | .false. | FIXME: Add documentation |
| navg | integer | 10 | FIXME: Add documentation |
| ncheck | integer | 10  | FIXME: Add documentation |
| nmovie | integer | -1 | FIXME: Add documentation |
| nsave | integer | -1 | FIXME: Add documentation |
| nwrite | integer | 10 | FIXME: Add documentation |
| nwrite_mult | integer | 10 | FIXME: Add documentation |
| ob_midplane | logical | .false. | FIXME: Add documentation |
| omegatinst | real | 1.0e6 | FIXME: Add documentation |
| omegatol | real | 1e-3 | FIXME: Add documentation |
| print_flux_line | logical | .false. | FIXME: Add documentation |
| print_line | logical | .false. | FIXME: Add documentation |
| save_distfn | logical | .false. | FIXME: Add documentation |
| save_for_restart | logical | .false. | FIXME: Add documentation |
| save_glo_info_and_grids | logical | .false.  | FIXME: Add documentation |
| save_many | logical | .false. | FIXME: Add documentation |
| save_velocities | logical | .false.  | FIXME: Add documentation |
| serial_netcdf4 | logical | .false.  | FIXME: Add documentation |
| use_nonlin_convergence | logical | .false. | FIXME: Add documentation |
| write_any | logical | .true.  | FIXME: Add documentation |
| write_apar_over_time | logical | .false. | FIXME: Add documentation |
| write_ascii | logical | .true. | FIXME: Add documentation |
| write_avg_moments | logical | .false. | FIXME: Add documentation |
| write_bpar_over_time | logical | .false. | FIXME: Add documentation |
| write_cerr | logical | .false. | FIXME: Add documentation |
| write_collisional | logical | .false.  | FIXME: Add documentation |
| write_correlation | logical | .true. | FIXME: Add documentation |
| write_correlation_extend | logical | .false. | FIXME: Add documentation |
| write_cross_phase | logical | .false. | FIXME: Add documentation |
| write_density_over_time | logical | .false.  | FIXME: Add documentation |
| write_eigenfunc | logical | .false. | FIXME: Add documentation |
| write_fields | logical | .true. | FIXME: Add documentation |
| write_final_antot | logical | .false. | FIXME: Add documentation |
| write_final_db | logical | .false. | FIXME: Add documentation |
| write_final_epar | logical | .false. | FIXME: Add documentation |
| write_final_fields | logical | .false. | FIXME: Add documentation |
| write_final_moments | logical | .false. | FIXME: Add documentation |
| write_flux_line | logical | .true. | FIXME: Add documentation |
| write_fluxes | logical | .false.  | FIXME: Add documentation |
| write_fluxes_by_mode | logical | .false.  | FIXME: Add documentation |
| write_full_moments_notgc | logical | .false. | FIXME: Add documentation |
| write_g | logical | .false. | FIXME: Add documentation |
| write_gg | logical | .false. | FIXME: Add documentation |
| write_gs | logical | .false. | FIXME: Add documentation |
| write_gyx | logical | .false. | FIXME: Add documentation |
| write_heating | logical | .false.  | FIXME: Add documentation |
| write_hrate | logical | .false. | FIXME: Add documentation |
| write_jext | logical | .false.  | FIXME: Add documentation |
| write_kpar | logical | .false. | FIXME: Add documentation |
| write_line | logical | .true. | FIXME: Add documentation |
| write_lorentzian | logical | .false. | FIXME: Add documentation |
| write_lpoly | logical | .false. | FIXME: Add documentation |
| write_max_verr | logical | .false. | FIXME: Add documentation |
| write_moments | logical | .true. | FIXME: Add documentation |
| write_nl_flux | logical | .true. | FIXME: Add documentation |
| write_nl_flux_dist | logical | .false. | FIXME: Add documentation |
| write_ntot_over_time | logical | .false.  | FIXME: Add documentation |
| write_omavg | logical | .true. | FIXME: Add documentation |
| write_omega | logical | .true. | FIXME: Add documentation |
| write_parity | logical | .false. | FIXME: Add documentation |
| write_pflux_sym | logical | .false.  | FIXME: Add documentation |
| write_pflux_tormom | logical | .false.  | FIXME: Add documentation |
| write_phi_over_time | logical | .false. | FIXME: Add documentation |
| write_symmetry | logical | .false. | FIXME: Add documentation |
| write_tperp_over_time | logical | .false.  | FIXME: Add documentation |
| write_upar_over_time | logical | .false.  | FIXME: Add documentation |
| write_verr | logical | .false. | FIXME: Add documentation |
| write_zonal_transfer | logical | .false.  | FIXME: Add documentation |