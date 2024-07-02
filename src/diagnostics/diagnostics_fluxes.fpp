#include "unused_dummy.inc"
!> Module which contains functions for calculating
!> and writing out the fluxes of heat and momentum etc.
module diagnostics_fluxes

  implicit none

  private

  public :: init_diagnostics_fluxes, finish_diagnostics_fluxes, reset_averages_and_counters
  public :: calculate_fluxes, common_calculate_fluxes
  public :: flux_vs_theta_vs_vpa, pflux_vs_theta_vs_vpa, flux_dist

  public :: qheat, qmheat, qbheat, pflux, vflux, vflux_par
  public :: vflux_perp, pflux_tormom, vflux0, vflux1, pmflux
  public :: vmflux, pbflux, vbflux, exchange

  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable ::  pflux_tormom
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:), allocatable :: esflux_vs_e, apflux_vs_e, bpflux_vs_e
  real, dimension (:,:,:), allocatable :: exchange
  real, dimension (:,:,:), allocatable :: exchange1

#ifdef NETCDF_PARALLEL
  logical, parameter :: moment_to_allprocs = .true.
#else
  logical, parameter :: moment_to_allprocs = .false.
#endif

contains
  !> Allocate and zero module-level arrays
  subroutine init_diagnostics_fluxes(gnostics)
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use le_grids, only: negrid
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout), optional :: gnostics

    if (present(gnostics)) then
      gnostics%current_results%species_heat_flux_avg = 0.0
      gnostics%current_results%species_momentum_flux_avg = 0.0
      gnostics%current_results%species_particle_flux_avg = 0.0
    end if

    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (pflux_tormom (ntheta0,naky,nspec)) ; pflux_tormom = 0.
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.

    allocate (exchange1 (ntheta0,naky,nspec)) ; exchange1 = 0.
    allocate (exchange (ntheta0,naky,nspec)) ; exchange = 0.

    allocate (vflux_par (ntheta0,naky,nspec)) ; vflux_par = 0.
    allocate (vflux_perp (ntheta0,naky,nspec)) ; vflux_perp = 0.

    allocate (vflux0 (ntheta0,naky,nspec)) ; vflux0 = 0.
    allocate (vflux1 (ntheta0,naky,nspec)) ; vflux1 = 0.

    allocate (pmflux(ntheta0,naky,nspec)) ; pmflux = 0.
    allocate (qmheat(ntheta0,naky,nspec,3)) ; qmheat = 0.
    allocate (vmflux(ntheta0,naky,nspec)) ; vmflux = 0.

    allocate (pbflux(ntheta0,naky,nspec)) ; pbflux = 0.
    allocate (qbheat(ntheta0,naky,nspec,3)) ; qbheat = 0.
    allocate (vbflux(ntheta0,naky,nspec)) ; vbflux = 0.

    allocate (esflux_vs_e(negrid,nspec)) ; esflux_vs_e = 0.
    allocate (apflux_vs_e(negrid,nspec)) ; apflux_vs_e = 0.
    allocate (bpflux_vs_e(negrid,nspec)) ; bpflux_vs_e = 0.
  end subroutine init_diagnostics_fluxes

  !> Clean up module, deallocate module-level arrays
  subroutine finish_diagnostics_fluxes
    implicit none
    if (.not. allocated(pflux)) return

    deallocate (pflux, pflux_tormom, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange, exchange1, &
         esflux_vs_e, apflux_vs_e, bpflux_vs_e)
  end subroutine finish_diagnostics_fluxes

  !> Set averages in gnostics\%current_results to 0.
  subroutine reset_averages_and_counters(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    gnostics%current_results%species_heat_flux_avg = 0.0
    gnostics%current_results%species_particle_flux_avg = 0.0
    gnostics%current_results%species_momentum_flux_avg = 0.0
  end subroutine reset_averages_and_counters

  !> Calculate and possibly write fluxes.
  !>
  !> The fluxes are calculated as a function of x, y and species. This
  !> function writes the whole array, and also various averages of them
  subroutine calculate_fluxes(gnostics)
    use run_parameters, only: has_phi, has_apar, has_bpar
    use nonlinear_terms, only: nonlin
    use diagnostics_config, only: diagnostics_type
    use unit_tests, only: debug_message
    use gs2_io, only: starts, egrid_dim, species_dim, time_dim
    use build_config, only: gs2_has_lowflow
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics

    call debug_message(4, 'diagnostics_fluxes::calculate_fluxes starting')
    gnostics%current_results%total_heat_flux = 0.0
    gnostics%current_results%total_momentum_flux = 0.0
    gnostics%current_results%total_particle_flux = 0.0
    gnostics%current_results%species_heat_flux = 0.0
    gnostics%current_results%species_momentum_flux = 0.0
    gnostics%current_results%species_particle_flux = 0.0

    if (.not. nonlin) call write_diffusive_estimates(gnostics)

    call common_calculate_fluxes()

    if (has_phi) then
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux',  'Turbulent flux of heat', 'Q_gB = ', qheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux_par',  'Turbulent flux of parallel heat', 'Q_gB = ', qheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_heat_flux_perp',  'Turbulent flux of perpendicular heat', 'Q_gB = ', qheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_part_flux',  'Turbulent flux of particles', 'n_r? ', pflux, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_part_tormom_flux',  'Ask Jung-Pyo Lee...', '? ', pflux_tormom, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux',  'Flux of toroidal angular momentum', 'Pi_gB =  ', vflux, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux_par', 'Flux of the parallel component of the toroidal angular momentum', 'Pi_gB =  ', &
            vflux_par, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_energy_exchange', '??', 'Pi_gB =  ', &
            exchange, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom_flux_perp', 'Flux of the perpendicular component of the toroidal angular momentum', 'Pi_gB =  ', &
            vflux_perp, gnostics%distributed)

#ifdef NETCDF
       if (gnostics%writing) call neasyf_write(gnostics%file_id, "es_flux_vs_e", esflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Electrostatic flux", start=starts(3, gnostics%nout))
#endif
       if (gs2_has_lowflow) then
          call calculate_standard_flux_properties(gnostics, &
               'es_mom0', 'Low-flow momentum flux 0 (Ask Michael Barnes)', 'Pi_gB =  ', &
               vflux0, gnostics%distributed)
          call calculate_standard_flux_properties(gnostics, &
               'es_mom1', 'Low-flow momentum flux 1 (Ask Michael Barnes)', 'Pi_gB =  ', &
               vflux0, gnostics%distributed)
       end if
    end if

    if (has_apar) then
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux',  'Turbulent flux of heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux_par',  'Turbulent flux of parallel heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_heat_flux_perp',  'Turbulent flux of perpendicular heat resulting from &
            & perpendicular magnetic fluctuations', 'Q_gB = ', qmheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_part_flux',  'Turbulent flux of particles resulting from &
            & perpendicular magnetic fluctuations', 'TBC ', pmflux(:,:,:), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'apar_mom_flux',  'Turbulent flux of momentum resulting from &
            & perpendicular magnetic fluctuations', 'Pi_gB = ', vmflux(:,:,:), gnostics%distributed)
#ifdef NETCDF
       if (gnostics%writing) call neasyf_write(gnostics%file_id, "apar_flux_vs_e", apflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Perpendicular magnetic flux, as function of energy", &
            start=starts(3, gnostics%nout))
#endif
    end if
    if (has_bpar) then
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux',  'Turbulent flux of heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,1), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux_par',  'Turbulent flux of parallel heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,2), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_heat_flux_perp',  'Turbulent flux of perpendicular heat resulting from &
            & parallel magnetic fluctuations', 'Q_gB = ', qbheat(:,:,:,3), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_part_flux',  'Turbulent flux of particles resulting from &
            & parallel magnetic fluctuations', 'TBC ', pbflux(:,:,:), gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'bpar_mom_flux',  'Turbulent flux of momentum resulting from &
            & parallel magnetic fluctuations', 'Pi_gB = ', vbflux(:,:,:), gnostics%distributed)
#ifdef NETCDF
       if (gnostics%writing) call neasyf_write(gnostics%file_id, "bpar_flux_vs_e", bpflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Parallel magnetic flux, as function of energy", &
            start=starts(3, gnostics%nout))
#endif
    end if

    ! Update averages
    gnostics%current_results%species_heat_flux_avg = gnostics%current_results%species_heat_flux_avg + &
         gnostics%current_results%species_heat_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_particle_flux_avg = gnostics%current_results%species_particle_flux_avg + &
         gnostics%current_results%species_particle_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_momentum_flux_avg = gnostics%current_results%species_momentum_flux_avg + &
         gnostics%current_results%species_momentum_flux * (gnostics%user_time-gnostics%user_time_old)

#ifdef NETCDF
    ! Write totals
    if (gnostics%write_fluxes .and. gnostics%writing) then
       call neasyf_write(gnostics%file_id, "heat_flux_tot", &
            gnostics%current_results%total_heat_flux, dim_names=[time_dim], start=[gnostics%nout], units="Q_gB", &
            long_name="Total heat flux")
       call neasyf_write(gnostics%file_id, "hflux_tot", &
            gnostics%current_results%total_heat_flux, dim_names=[time_dim], start=[gnostics%nout], units="Q_gB", &
            long_name="Total heat flux, same as heat_flux_tot, included for backwards compatiblity")
       call neasyf_write(gnostics%file_id, "mom_flux_tot", &
            gnostics%current_results%total_momentum_flux, dim_names=[time_dim], start=[gnostics%nout], units="Pi_gB", &
            long_name="Total momentum flux")
       call neasyf_write(gnostics%file_id, "vflux_tot", &
            gnostics%current_results%total_momentum_flux, dim_names=[time_dim], start=[gnostics%nout], units="Pi_gB", &
            long_name="Total momentum flux, same as mom_flux_tot, included for backwards compatiblity")
       call neasyf_write(gnostics%file_id, "part_flux_tot", &
            gnostics%current_results%total_particle_flux, dim_names=[time_dim], start=[gnostics%nout], &
            long_name="Total particle flux")
       call neasyf_write(gnostics%file_id, "zflux_tot", &
            gnostics%current_results%total_particle_flux, dim_names=[time_dim], start=[gnostics%nout], &
            long_name="Total particle flux, same as part_flux_tot, included for backwards compatiblity")
    end if
#endif
  end subroutine calculate_fluxes

  !> Calculate heat, particle, and momentum fluxes
  !>
  !> Shared between old and new diagnostics
  subroutine common_calculate_fluxes()
    use build_config, only: gs2_has_lowflow
    use mp, only: proc0
    use species, only: spec
    use dist_fn_arrays, only: g_adjust, gnew, to_g_gs2, from_g_gs2
    use species, only: nspec, spec
    use fields_arrays, only: phinew, bparnew, aparnew, phi
    use run_parameters, only: has_phi, has_apar, has_bpar
    use diagnostics_config, only: diagnostics_type
    use unit_tests, only: debug_message

    integer :: is

    call g_adjust (gnew, phinew, bparnew, direction = from_g_gs2)
    call flux (phinew, aparnew, bparnew, &
         pflux,  qheat,  vflux, vflux_par, vflux_perp, &
         pmflux, qmheat, vmflux, pbflux, qbheat, vbflux, pflux_tormom)
    ! Only used in new diagnostics, can't turn off separately to fluxes in general
    call flux_vs_e (phinew, aparnew, bparnew, esflux_vs_e, apflux_vs_e, bpflux_vs_e)

    ! lowflow terms only implemented in electrostatic limit at present
    if (gs2_has_lowflow) call lf_flux (phinew, vflux0, vflux1)

    call g_adjust (gnew, phinew, bparnew, direction = to_g_gs2)

    ! Note we don't use / write exchange1 anywhere, we also call with gnew unadjusted
    if (has_phi) call eexchange (phinew, phi, exchange1, exchange)

    ! Note all procs have the main results so could also do the below such that
    ! we have the right result on all processors.
    if (.not. proc0) return

    do is = 1, nspec
      associate(species => spec(is), &
           density_temperature => spec(is)%dens * spec(is)%temp, &
           momentum => spec(is)%dens * sqrt(spec(is)%mass * spec(is)%temp), &
           momentum_perp => spec(is)%dens * spec(is)%mass * spec(is)%stm &
        )

        if (has_phi) then
          qheat(:,:,is,1) = qheat(:,:,is,1) * density_temperature
          qheat(:,:,is,2) = qheat(:,:,is,2) * density_temperature
          qheat(:,:,is,3) = qheat(:,:,is,3) * density_temperature

          pflux(:,:,is) = pflux(:,:,is) * species%dens
          pflux_tormom(:,:,is) = pflux_tormom(:,:,is) * species%dens

          esflux_vs_e(:,is) = esflux_vs_e(:,is) * spec(is)%dens

          vflux(:,:,is) = vflux(:,:,is) * momentum
          vflux_par(:,:,is) = vflux_par(:,:,is) * momentum
          vflux_perp(:,:,is) = vflux_perp(:,:,is) * momentum_perp

          exchange(:,:,is) = exchange(:,:,is) * species%dens * species%z

          if (gs2_has_lowflow) then
            vflux0(:,:,is) = vflux0(:,:,is) * momentum
            vflux1(:,:,is) = vflux1(:,:,is) * species%dens * species%mass * species%tz
          end if
        end if
        if (has_apar) then
          qmheat(:,:,is,1) = qmheat(:,:,is,1) * density_temperature
          qmheat(:,:,is,2) = qmheat(:,:,is,2) * density_temperature
          qmheat(:,:,is,3) = qmheat(:,:,is,3) * density_temperature
          pmflux(:,:,is) = pmflux(:,:,is) * species%dens
          vmflux(:,:,is) = vmflux(:,:,is) * momentum_perp
          apflux_vs_e(:,is) = apflux_vs_e(:,is) * spec(is)%dens
        end if
        if (has_bpar) then
          qbheat(:,:,is,1) = qbheat(:,:,is,1) * density_temperature
          qbheat(:,:,is,2) = qbheat(:,:,is,2) * density_temperature
          qbheat(:,:,is,3) = qbheat(:,:,is,3) * density_temperature
          pbflux(:,:,is) = pbflux(:,:,is) * species%dens
          vbflux(:,:,is) = vbflux(:,:,is) * momentum_perp
          bpflux_vs_e(:,is) = bpflux_vs_e(:,is) * spec(is)%dens
        end if
      end associate
    end do
  end subroutine common_calculate_fluxes

  !> Calculate estimates of the heat and particles fluxes using
  !! gamma / k^2 estimate of the diffusivity
  subroutine write_diffusive_estimates(gnostics)
    use diagnostics_omega, only: omega_average
    use fields_parallelization, only: field_k_local
    use species, only: spec, nspec
    use kt_grids, only: kperp2, ntheta0, naky
    use theta_grid, only: grho
    use diagnostics_config, only: diagnostics_type
    use mp, only: sum_allreduce
    use gs2_io, only: flux_dims, starts, time_dim
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    real, dimension(ntheta0, naky) :: diffusivity_by_k
    real :: heat_flux_max
    real, dimension(nspec) :: heat_flux_max_by_spec
    real :: particle_flux_max
    real, dimension(nspec) :: particle_flux_max_by_spec
    real, dimension(ntheta0, naky, nspec) :: heat_flux
    real, dimension(ntheta0, naky, nspec) :: particle_flux
    real, dimension(naky, nspec) :: heat_flux_by_ky
    real, dimension(naky, nspec) :: particle_flux_by_ky
    real, dimension(ntheta0, naky) :: momentum_flux
    integer :: it, ik, is

    diffusivity_by_k = 0.0
    heat_flux = 0.0
    particle_flux = 0.0
    momentum_flux = 0.0
    heat_flux_max = 0.0

    do ik = 1,naky
       do it = 1,ntheta0
          if (.not. gnostics%distributed .or. field_k_local(it,ik)) then
             if (kperp2(gnostics%igomega,it,ik).eq.0.0) cycle
             diffusivity_by_k(it,ik) = &
                  max(aimag(omega_average(it,ik)),0.0)/kperp2(gnostics%igomega, it, ik)*2.0
             do is = 1,nspec
                ! Q = n chi grad T = n (gamma / k^2) dT / dr
                ! = dens  n_r (gamma_N v_thr / k_N**2 rho_r a) dT / drho drho/dr
                ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) T a / L_T drho/dr
                ! = dens  n_r (gamma_N v_thr rho_r **2 / k_N**2 a) temp T_r tprim drho/dr_N/a
                ! Q / (n_r  v_r T_r rho_r**2/a**2)
                ! = dens (gamma_N / k_N**2) temp tprim grho
                !
                heat_flux(it,ik,is) = diffusivity_by_k(it,ik) * &
                     spec(is)%dens * spec(is)%temp * spec(is)%tprim *  grho(gnostics%igomega)
                particle_flux(it,ik,is) = diffusivity_by_k(it,ik) * &
                     spec(is)%dens **2.0 * spec(is)%fprim * grho(gnostics%igomega)
             end do
          end if
       end do
    end do

    if (gnostics%distributed) call sum_allreduce(diffusivity_by_k)

    gnostics%current_results%diffusivity = &
      maxval(diffusivity_by_k) * grho(gnostics%igomega)

    call calculate_standard_flux_properties(gnostics, &
         'heat_flux_diff',  'Diffusive estimate of turbulent flux of heat', &
         'Q_gB = ', heat_flux, gnostics%distributed)

    call calculate_standard_flux_properties(gnostics, &
         'part_flux_diff',  'Diffusive estimate of turbulent flux of particles', &
         'n_r? ', particle_flux, gnostics%distributed)

    heat_flux_by_ky = maxval(heat_flux, 1)
    particle_flux_by_ky = maxval(particle_flux, 1)

    heat_flux_max_by_spec = maxval(heat_flux_by_ky, 1)
    particle_flux_max_by_spec = maxval(particle_flux_by_ky, 1)

    heat_flux_max = sum(heat_flux_max_by_spec)
    particle_flux_max = sum(particle_flux_max_by_spec)

#ifdef NETCDF
    if (gnostics%write_fluxes .and. gnostics%writing) then
       call neasyf_write(gnostics%file_id, "es_heat_flux_diff_max", heat_flux_max_by_spec, &
            dim_names=flux_dims, long_name="Linear estimate of the heat flux", units="Q_gB", &
            start=starts(2, gnostics%nout))
       call neasyf_write(gnostics%file_id, "heat_flux_diff_max", heat_flux_max, &
            dim_names=[time_dim], long_name="Linear estimate of the heat flux", units="Q_gB", &
            start=[gnostics%nout])
       call neasyf_write(gnostics%file_id, "es_particle_flux_diff_max", particle_flux_max_by_spec, &
            dim_names=flux_dims, long_name="Linear estimate of the particle flux", units="Q_gB", &
            start=starts(2, gnostics%nout))
       call neasyf_write(gnostics%file_id, "particle_flux_diff_max", particle_flux_max, &
            dim_names=[time_dim], long_name="Linear estimate of the particle flux", units="Q_gB", &
            start=[gnostics%nout])
       call neasyf_write(gnostics%file_id, "diffusivity", gnostics%current_results%diffusivity, &
            dim_names=[time_dim], long_name="Linear estimate of the diffusivity", &
            start=[gnostics%nout])
    end if
#endif
  end subroutine write_diffusive_estimates

  !> Writes a range of different summed and averaged properties of the given
  !! flux... e.g. the flux summed over kx as a function of ky, species and time
  subroutine calculate_standard_flux_properties(gnostics, flux_name, flux_description, &
       flux_units, flux_value, distributed)
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use fields_parallelization, only: field_k_local
    use mp, only: broadcast, sum_allreduce
    use volume_averages, only: average_all, average_kx, average_ky
    use diagnostics_config, only: diagnostics_type
    use gs2_io, only: flux_dims, fluxk_dims, mode_dims, starts, kx_dim, ky_dim, time_dim
#ifdef NETCDF
    use neasyf, only: neasyf_write
#endif
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    character(*), intent(in) :: flux_name, flux_description, flux_units
    real, dimension(ntheta0,naky,nspec), intent(inout) :: flux_value
    logical, intent(in) :: distributed
    real, dimension(ntheta0, naky) :: total_flux_by_mode
    real, dimension(naky) :: total_flux_by_ky
    real, dimension(ntheta0) :: total_flux_by_kx
    real, dimension(nspec) :: flux_by_species
    integer :: is, ik, it

    call average_all(flux_value, flux_by_species, distributed)

    total_flux_by_mode = 0.
    do ik = 1,naky
      do it = 1,ntheta0
        if (.not. distributed .or. field_k_local(it, ik)) then
          do is = 1,nspec
            total_flux_by_mode(it, ik) = &
                 total_flux_by_mode(it, ik) + flux_value(it, ik, is)
          end do
        end if
      end do
    end do

    call average_kx(total_flux_by_mode, total_flux_by_ky, distributed)
    call average_ky(total_flux_by_mode, total_flux_by_kx, distributed)
#ifdef NETCDF
    if (gnostics%write_fluxes .and. gnostics%writing) then
      call neasyf_write(gnostics%file_id, flux_name, flux_by_species, &
           dim_names=flux_dims, start=starts(2, gnostics%nout), &
           long_name=flux_description//" averaged over kx and ky", units=flux_units)

      call neasyf_write(gnostics%file_id, "total_"//flux_name//"_by_ky", total_flux_by_ky, &
           dim_names=[ky_dim, time_dim], &
           long_name=flux_description//" summed over species and averaged over kx", &
           units=flux_units, start=starts(2, gnostics%nout))

      call neasyf_write(gnostics%file_id, "total_"//flux_name//"_by_kx", total_flux_by_kx, &
           dim_names=[kx_dim, time_dim], &
           long_name=flux_description//" summed over species and averaged over ky", &
           units=flux_units, start=starts(2, gnostics%nout))

      call neasyf_write(gnostics%file_id, "total_"//flux_name, sum(total_flux_by_kx), &
           dim_names=[time_dim], &
           long_name=flux_description//" summed over species and averaged over kx and ky", &
           units=flux_units, start=[gnostics%nout])

      if (gnostics%write_fluxes_by_mode) then
        call neasyf_write(gnostics%file_id, flux_name//"_by_mode", flux_value, &
             dim_names=fluxk_dims, &
             long_name=flux_description//" as a function of species, kx and ky", &
             units=flux_units, start=starts(4, gnostics%nout))
        call neasyf_write(gnostics%file_id, "total_"//flux_name//"_by_mode", total_flux_by_mode, &
             dim_names=mode_dims, &
             long_name=flux_description//" summed over species", &
             units=flux_units, start=starts(3, gnostics%nout))
      end if
    end if
#endif

    if (flux_name .eq. 'es_heat_flux') gnostics%current_results%species_es_heat_flux = flux_by_species
    if (flux_name .eq. 'es_energy_exchange') gnostics%current_results%species_energy_exchange = flux_by_species
    if (flux_name .eq. 'apar_heat_flux') gnostics%current_results%species_apar_heat_flux = flux_by_species
    if (flux_name .eq. 'bpar_heat_flux') gnostics%current_results%species_bpar_heat_flux = flux_by_species

    if (flux_name .eq. 'es_heat_flux' &
         .or. flux_name .eq. 'apar_heat_flux' &
         .or. flux_name .eq. 'bpar_heat_flux') then
      gnostics%current_results%total_heat_flux = gnostics%current_results%total_heat_flux + sum(flux_by_species)
      gnostics%current_results%species_heat_flux = gnostics%current_results%species_heat_flux + flux_by_species
    else if (flux_name .eq. 'es_mom_flux' &
                                !.or.flux_name .eq. 'es_mom0' & ! Low flow fluxes, currently disabled
         .or.flux_name .eq. 'apar_mom_flux' &
         .or.flux_name .eq. 'bpar_mom_flux') then
      gnostics%current_results%total_momentum_flux = gnostics%current_results%total_momentum_flux + sum(flux_by_species)
      gnostics%current_results%species_momentum_flux = gnostics%current_results%species_momentum_flux + flux_by_species
    else if (flux_name .eq. 'es_part_flux' &
         .or.flux_name .eq. 'apar_part_flux' &
         .or.flux_name .eq. 'bpar_part_flux') then
      gnostics%current_results%total_particle_flux = gnostics%current_results%total_particle_flux + sum(flux_by_species)
      gnostics%current_results%species_particle_flux = gnostics%current_results%species_particle_flux + flux_by_species
    end if
  end subroutine calculate_standard_flux_properties

  !> Calculate various fluxes
  !>
  !> FIXME: add documentation of exactly what fluxes
  !> @note: Requires g_adjust to be called on gnew (from_g_gs2)
  !> before calling this routine. We _could_ move the call here
  !> to improve safety at the expense of more calls to g_adjust.
  subroutine flux (phi, apar, bpar, &
       pflux,  qflux,  vflux, vflux_par, vflux_perp, &
       pmflux, qmflux, vmflux, &
       pbflux, qbflux, vbflux, pflux_tormom)

    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, vpa, aj1, vperp2, g_work
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use run_parameters, only: has_phi, has_apar, has_bpar
    use constants, only: zi
    use theta_grid, only: Rplot, Bpol, rhoc
    use build_config, only: gs2_has_lowflow
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:,:), intent (out) :: pflux, pmflux, pbflux, pflux_tormom
    real, dimension (:,:,:), intent (out) :: vflux, vmflux, vbflux, vflux_par, vflux_perp
    real, dimension (:,:,:,:), intent (out) :: qflux, qmflux, qbflux
    integer :: it, ik, is, isgn, ig
    integer :: iglo

    !CMR, 15/1/08:
    !  Implemented Clemente Angioni's fix for fluxes by replacing g with gnew
    !  so fields and distribution function are evaluated self-consistently in time.
    !  This fixed unphysical oscillations in non-ambipolar particle fluxes
    !

    pflux = 0.0;   qflux = 0.0;   vflux = 0.0 ; vflux_par = 0.0 ; vflux_perp = 0.0
    pmflux = 0.0;  qmflux = 0.0;  vmflux = 0.0
    pbflux = 0.0;  qbflux = 0.0;  vbflux = 0.0
    pflux_tormom = 0.0

    if (has_phi) then
       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do

       call get_flux (g_work, phi, pflux)
       if (gs2_has_lowflow) call get_flux_tormom (g_work, phi, pflux_tormom)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do

       call get_flux (g_work, phi, qflux(:,:,:,1))

       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*2.*vpa(:,isgn,:)**2*aj0
       end do

       call get_flux (g_work, phi, qflux(:,:,:,2))

       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*vperp2*aj0
       end do

       call get_flux (g_work, phi, qflux(:,:,:,3))

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g_work(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do

       call get_flux (g_work, phi, vflux_par)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) = -zi*aky(ik)*gnew(:,isgn,iglo)*aj1(:,iglo) &
                  *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,iglo)*spec(is)%smz/(qval*shat*bmag**2)
          end do
       end do

       call get_flux (g_work, phi, vflux_perp)

       vflux = vflux_par + vflux_perp

    else
       pflux = 0.
       qflux = 0.
       vflux = 0.
    end if

    if (has_apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do

       call get_flux (g_work, apar, pmflux)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do

       call get_flux (g_work, apar, qmflux(:,:,:,1))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *2.*vpa(:,isgn,iglo)**2
          end do
       end do

       call get_flux (g_work, apar, qmflux(:,:,:,2))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo) &
                  *vperp2(:,iglo)
          end do
       end do

       call get_flux (g_work, apar, qmflux(:,:,:,3))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm &
                  *vpa(:,isgn,iglo)*vpac(:,isgn,iglo)
          end do
       end do

       call get_flux (g_work, apar, vmflux)

    else
       pmflux = 0.
       qmflux = 0.
       vmflux = 0.
    end if

    if (has_bpar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do

       call get_flux (g_work, bpar, pbflux)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do

       call get_flux (g_work, bpar, qbflux(:,:,:,1))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *2.*vpa(:,isgn,iglo)**2
          end do
       end do

       call get_flux (g_work, bpar, qbflux(:,:,:,2))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz &
                    *vperp2(:,iglo)
          end do
       end do

       call get_flux (g_work, bpar, qbflux(:,:,:,3))

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo) &
                  *spec(is)%tz*vpac(:,isgn,iglo)
          end do
       end do

       call get_flux (g_work, bpar, vbflux)

    else
       pbflux = 0.
       qbflux = 0.
       vbflux = 0.
    end if
  end subroutine flux

  !> FIXME : Add documentation
  subroutine flux_vs_e (phi, apar, bpar, pflux, pmflux, pbflux)
    use species, only: spec
    use theta_grid, only: ntgrid
    use dist_fn_arrays, only: gnew, aj0, vpa, aj1, vperp2, g_work
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only: has_phi, has_apar, has_bpar
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, apar, bpar
    real, dimension (:,:), intent (inout) :: pflux, pmflux, pbflux
    integer :: is, isgn
    integer :: iglo

    if (proc0) then
        pflux = 0.0;   pmflux = 0.0; pbflux = 0.0 ;
    end if

    if (has_phi) then
       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do

       call get_flux_vs_e (g_work, phi, pflux)
    else
       pflux = 0.
    end if

    if (has_apar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = -gnew(:,isgn,iglo)*aj0(:,iglo)*spec(is)%stm*vpa(:,isgn,iglo)
          end do
       end do

       call get_flux_vs_e (g_work, apar, pmflux)

    else
       pmflux = 0.
    end if

    if (has_bpar) then
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) &
                  = gnew(:,isgn,iglo)*aj1(:,iglo)*2.0*vperp2(:,iglo)*spec(is)%tz
          end do
       end do

       call get_flux_vs_e (g_work, bpar, pbflux)
    else
       pbflux = 0.
    end if

  end subroutine flux_vs_e

  !> Calculate the flux of a field
  subroutine get_flux (g_in, fld, flx)
    use theta_grid, only: ntgrid, grho, theta, jacob
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use integration, only: trapezoidal_integration
    use gs2_layouts, only: g_lo
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    !> Input weighted distribution
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    !> Input field
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    !> Output flux
    real, dimension (:,:,:), intent (out) :: flx
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))

    call integrate_moment (g_in, total, moment_to_allprocs, full_arr)

    ! This is essentially 2 pi / surfarea
    wgt = 1.0/trapezoidal_integration(theta, grho * jacob)

    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             flx(it,ik,is) = trapezoidal_integration(theta, aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) * jacob)*wgt*aky(ik)
          end do
       end do
    end do
    flx = flx*0.5

    deallocate (total)

  end subroutine get_flux

  !> FIXME : Add documentation
  subroutine get_flux_vs_e (g_in, fld, flx)
    use theta_grid, only: ntgrid, grho, theta, jacob
    use kt_grids, only: aky
    use le_grids, only: negrid, wl
    use species, only: nspec
    use mp, only: sum_reduce, proc0
    use gs2_layouts, only: g_lo,ie_idx,il_idx,is_idx,it_idx,ik_idx,isign_idx
    use integration, only: trapezoidal_integration
    use gs2_layouts, only: g_lo
    implicit none
    !> Input weighted distribution
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:), intent (inout) :: flx
    real, dimension (:,:), allocatable :: total
    real:: wgt,fac
    integer :: ik, it, is, iglo, ie, il

    allocate(total(negrid,nspec))
    total = 0.

    ! This is essentially 2 pi / surfarea
    wgt = 1.0/trapezoidal_integration(theta, grho * jacob)

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       fac =0.5
       ie = ie_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if (aky(ik) == 0.) fac = 1.0

       total(ie,is) = total(ie,is) + fac*&
            trapezoidal_integration(theta, &
            aimag(g_in(:,1,iglo)*conjg(fld(:,it,ik))) &
            * wl(:,il) * jacob) * aky(ik) * wgt
       total(ie,is) = total(ie,is) + fac*&
            trapezoidal_integration(theta, &
            aimag(g_in(:,2,iglo)*conjg(fld(:,it,ik))) &
            * wl(:,il) * jacob) * aky(ik) * wgt
    end do

    call sum_reduce(total,0)

    ! We've already weighted by 0.5 in the above, do we need a second 0.5 here?
    if (proc0) flx = 0.5*total

  end subroutine get_flux_vs_e

  !> Identical to get_flux except don't integrate over poloidal angle ! JRB
  !>
  !> @todo: It would be nice if get_flux could calculate
  !> the theta average _and_ the theta distribution in order to replace the
  !> duplication associated with this routine.
  subroutine get_flux_dist (g_in, fld, flx_dist)
    use theta_grid, only: ntgrid, grho, theta, bmag, gradpar
    use kt_grids, only: ntheta0, aky, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use integration, only: trapezoidal_integration
    use gs2_time, only: woutunits
    use gs2_layouts, only: g_lo
    implicit none
    logical, parameter :: full_arr=moment_to_allprocs
    !> Input weighted distribution
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (-ntgrid:,:,:,:), intent (in out) :: flx_dist
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    ! EGH added final parameter 'all'
    ! to the call below for new parallel output
    ! This is temporary until distributed fields
    ! are implemented 1/2014
    call integrate_moment (g_in, total, moment_to_allprocs, full_arr)

    !Note that jacob = 1/(bmag*gradpar*drhodpsi). This is similar to the factor
    !we have below without the drhodpsi.
    wgt = 1.0/trapezoidal_integration(theta, grho / (bmag * gradpar))

    ! Note: This doesn't include the field line jacobian (jacob) and as such
    ! it's not clear that the above wgt is correct (as it uses 1/(bmag*gradpar)
    ! which is jacob * drhodpsi, rather than jacob as one might expect).
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             flx_dist(:,it,ik,is) = aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                  *aky(ik) * wgt / woutunits(ik)
          end do
       end do
    end do
    flx_dist = flx_dist*0.5

    deallocate (total)

  end subroutine get_flux_dist

  !> Diagnose the poloidal distribution of the particle, angular momentum,
  !> and energy fluxes
  !>
  !> JRB
  !>
  !> @note Only considers phi related fluxes. Some degree of duplication with flux.
  !> Ideally the flux and get_flux methods could return both the theta dependent and
  !> theta averaged fluxes.
  subroutine flux_dist (phi, bpar, pflux_dist, vflux_par_dist, vflux_perp_dist, qflux_dist)
    !CMR, 15/1/08:
    !  Implemented Clemente Angioni's fix for fluxes by replacing g with gnew
    !  so fields and distribution function are evaluated self-consistently in time.
    !  This fixed unphysical oscillations in non-ambipolar particle fluxes
    use species, only: spec
    use theta_grid, only: ntgrid, bmag
    use theta_grid, only: qval, shat, gds21, gds22
    use kt_grids, only: theta0, aky
    use le_grids, only: energy
    use dist_fn_arrays, only: gnew, aj0, vpac, aj1, vperp2, g_work
    use dist_fn_arrays, only: g_adjust, from_g_gs2, to_g_gs2
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use run_parameters, only:has_phi
    use constants, only: zi
    use theta_grid, only: Rplot, Bpol, rhoc
    use array_utils, only: zero_array
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi, bpar
    real, dimension (-ntgrid:,:,:,:), intent (out) :: pflux_dist, vflux_par_dist, vflux_perp_dist, qflux_dist
    integer :: it, ik, is, isgn, ig, iglo

    call zero_array(pflux_dist) ; call zero_array(vflux_par_dist)
    call zero_array(vflux_perp_dist) ; call zero_array(qflux_dist)
    call g_adjust (gnew, phi, bpar, direction = from_g_gs2) ! convert from g to h
    if (has_phi) then
       do isgn = 1, 2
          g_work(:,isgn,:) = gnew(:,isgn,:)*aj0
       end do
       call get_flux_dist (g_work, phi, pflux_dist)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          g_work(:,:,iglo) = g_work(:,:,iglo)*energy(ie_idx(g_lo,iglo),is_idx(g_lo,iglo))
       end do
       call get_flux_dist (g_work, phi, qflux_dist)

       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             ! Not clear that this should use the grid-centre v|| (vpac)
             ! as this is not defined at ig = ntgrid and all other quantities
             ! are evaluated on the grid points.
             g_work(ig,isgn,:) = gnew(ig,isgn,:)*aj0(ig,:)*vpac(ig,isgn,:)*Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2)
          end do
       end do
       call get_flux_dist (g_work, phi, vflux_par_dist)

       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          it = it_idx(g_lo,iglo)
          ik = ik_idx(g_lo,iglo)
          is = is_idx(g_lo,iglo)
          do isgn = 1, 2
             g_work(:,isgn,iglo) = -zi*aky(ik)*gnew(:,isgn,iglo)*aj1(:,iglo) &
                  *rhoc*(gds21+theta0(it,ik)*gds22)*vperp2(:,iglo)*spec(is)%smz/(qval*shat*bmag**2)
          end do
       end do
       call get_flux_dist (g_work, phi, vflux_perp_dist)
    end if
    call g_adjust (gnew, phi, bpar, direction = to_g_gs2) ! convert back from h to g
  end subroutine flux_dist

  !> Calculate toroidal angular momentum flux due to field
  !>
  !> FIXME: confirm & expand
  subroutine get_flux_tormom (g_in, fld, flx)
    use theta_grid, only: ntgrid
#ifdef LOWFLOW
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use kt_grids, only: aky
    use theta_grid, only: theta, grho, bmag, jacob
    use lowflow, only: mach_lab
    use integration, only: trapezoidal_integration
#endif
    use gs2_layouts, only: g_lo
    implicit none
    !> Input weighted distribution
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    !> Input field
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    !> Output flux
    real, dimension (:,:,:), intent (out) :: flx
#ifdef LOWFLOW
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is
    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    call integrate_moment (g_in, total)

    ! This is essentially 2 pi / surfarea
    wgt = 1.0/trapezoidal_integration(theta, grho * jacob)

    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             flx(it,ik,is) = trapezoidal_integration(theta, aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                  *Rplot(:)**2 * jacob)*aky(ik)*wgt*mach_lab(is)
          end do
       end do
    end do
    flx = flx*0.5

    deallocate (total)
#else
    UNUSED_DUMMY(g_in); UNUSED_DUMMY(fld)
    flx = 0.0
#endif
  end subroutine get_flux_tormom !JPL

  !> Calculate energy exchange diagnostic that numerically guarantees that the total
  !> energy exchange (summed over species) is zero
  !>
  !> @note This diagnostic can be relatively expensive to calculate.
  !>
  !> @note We currently call this with an unadjusted gnew and g, is that correct?
  subroutine eexchange (phinew, phi, exchange1, exchange2)
    use constants, only: zi
    use gs2_layouts, only: g_lo, il_idx, ie_idx, it_idx, ik_idx, is_idx
    use gs2_time, only: code_dt
    use dist_fn, only: wdrift_func
    use dist_fn_arrays, only: gnew, aj0, vpac, g, g_work
    use theta_grid, only: ntgrid, gradpar, delthet, field_line_average
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment, forbid
    use species, only: spec, nspec
    use nonlinear_terms, only: nonlin
    use hyper, only: hypervisc_filter
    use array_utils, only: zero_array
    implicit none

    complex, dimension (-ntgrid:,:,:), intent (in) :: phinew, phi
    real, dimension (:,:,:), intent (out) :: exchange1, exchange2
    integer :: ig, il, ie, it, ik, is, iglo, isgn
    complex :: dgdt_hypervisc
    complex, dimension (:,:,:,:), allocatable :: total, total2

    allocate (total(-ntgrid:ntgrid, ntheta0, naky, nspec))
    allocate (total2(-ntgrid:ntgrid, ntheta0, naky, nspec))
    ! Currently zero_array only distributes over the last dimension (nspec here)
    ! so there's limited parallelisation available here.
    call zero_array(total) ; call zero_array(total2)

    ! Zero elements not set in main loops
    if (nonlin) then
       exchange1(1, 1, :) = 0.0
       exchange2(1, 1, :) = 0.0
    end if

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, it, ik, is, il, ie, isgn, ig, dgdt_hypervisc) &
    !$OMP SHARED(g_lo, nonlin, g_work, ntgrid, forbid, aj0, code_dt, &
    !$OMP gnew, spec, hypervisc_filter, vpac, gradpar, delthet) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       if (nonlin .and. it==1 .and. ik==1) then
          g_work(:, :, iglo) = 0.0
          cycle
       end if
       is = is_idx(g_lo,iglo)
       il = il_idx(g_lo,iglo)
       ie = ie_idx(g_lo,iglo)
       do isgn = 1, 2
          ! get v_magnetic piece of g_work at grid points instead of cell centers
          do ig = -ntgrid, ntgrid
             if (forbid(ig, il)) then
                g_work(ig, isgn, iglo) = 0.0
             else
                ! It would be nice to cache the aj0*zi*wdrift_func*spec%tz
                ! constant here as wdrift_func can be relatively expensive
                g_work(ig,isgn,iglo) = aj0(ig,iglo)*(zi*wdrift_func(ig, iglo)/code_dt &
                     * gnew(ig,isgn,iglo)*spec(is)%tz)
             end if
          end do

          ! add contribution to g_work from hyperviscosity at grid points
          ! this is -(dg/dt)_hypervisc, equivalent to collisions term in eqn. 5 of PRL 109, 185003
          do ig = -ntgrid, ntgrid
             if (abs(hypervisc_filter(ig,it,ik)-1.0) > epsilon(0.0)) then
                dgdt_hypervisc = (1.0-1./hypervisc_filter(ig,it,ik))*gnew(ig,isgn,iglo)/code_dt
                ! should gnew be (gnew+gold)/2?
                g_work(ig,isgn,iglo) = g_work(ig,isgn,iglo) - dgdt_hypervisc
             end if
          end do

          ! get v_magnetic piece of g_work at cell centers and add in vpar piece at cell centers
          do ig = -ntgrid, ntgrid-1
             g_work(ig,isgn,iglo) = 0.5*(g_work(ig,isgn,iglo)+g_work(ig+1,isgn,iglo)) &
                  + 0.5*vpac(ig,isgn,iglo)*(gradpar(ig)+gradpar(ig+1))/delthet(ig) &
                  * (gnew(ig+1,isgn,iglo)-gnew(ig,isgn,iglo))*spec(is)%stm
          end do

       end do
    end do
    !$OMP END PARALLEL DO

    call integrate_moment (g_work, total)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(is, ik, it) &
    !$OMP SHARED(nspec, naky, ntheta0, nonlin, exchange1, total, phinew) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             if (nonlin .and. it==1 .and. ik==1) cycle
             exchange1(it,ik,is) = field_line_average( &
                  real(total(:,it,ik,is)*conjg(phinew(:,it,ik))))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, isgn) &
    !$OMP SHARED(g_lo, g_work, aj0, gnew, g) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*0.25*(gnew(:,isgn,iglo)+g(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (g_work, total)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(iglo, isgn) &
    !$OMP SHARED(g_lo, g_work, aj0, gnew, g) &
    !$OMP SCHEDULE(static)
    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       do isgn = 1, 2
          g_work(:,isgn,iglo) = aj0(:,iglo)*0.25*(gnew(:,isgn,iglo)-g(:,isgn,iglo))
       end do
    end do
    !$OMP END PARALLEL DO
    call integrate_moment (g_work, total2)

    ! exchange2 is a symmetrized form of energy exchange,
    ! which guarantees species-summed energy exchange is zero
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(is, ik, it) &
    !$OMP SHARED(nspec, naky, ntheta0, nonlin, exchange2, total, &
    !$OMP phinew, total2, phi, code_dt) &
    !$OMP COLLAPSE(3) &
    !$OMP SCHEDULE(static)
    do is = 1, nspec
       do ik = 1, naky
          do it = 1, ntheta0
             if (nonlin .and. it==1 .and. ik==1) cycle
             exchange2(it,ik,is) = field_line_average( &
                  real(total(:,it,ik,is) &
                  *conjg(phinew(:,it,ik)-phi(:,it,ik)) &
                  - (phinew(:,it,ik)+phi(:,it,ik))*conjg(total2(:,it,ik,is))) &
                  )/code_dt
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate (total, total2)

  end subroutine eexchange

  !> FIXME : Add documentation
  !!
  !! @Warning This routine may only make sense when compiling with LOWFLOW on?
  subroutine lf_flux (phi, vflx0, vflx1)
    use species, only: nspec
    use theta_grid, only: ntgrid, gradpar, delthet
    use theta_grid, only: drhodpsi, IoB
    use kt_grids, only: naky, ntheta0
    use dist_fn_arrays, only: gnew, aj0, vpa, g_work
    use gs2_layouts, only: g_lo, ie_idx, is_idx, it_idx, ik_idx
    use mp, only: proc0
    use run_parameters, only:has_phi, rhostar
    use constants, only: zi
    implicit none
    complex, dimension (-ntgrid:,:,:), intent (in) :: phi
    real, dimension (:,:,:), intent (out) :: vflx0, vflx1
    real, dimension (:,:), allocatable :: dum
    complex, dimension (:,:,:), allocatable :: dphi
    integer :: isgn, ig
    integer :: iglo

    allocate (dum (-ntgrid:ntgrid,nspec))
    allocate (dphi (-ntgrid:ntgrid,ntheta0,naky))

    if (proc0) then
       vflx0 = 0.0 ; vflx1 = 0.0 ; dum = 0.0
    end if

    do ig = -ntgrid, ntgrid-1
       dphi(ig,:,:) = (phi(ig+1,:,:)-phi(ig,:,:))/delthet(ig)
    end do
    ! not sure if this is the right way to handle ntgrid point -- MAB
    dphi(ntgrid,:,:) = (phi(ntgrid,:,:)-phi(ntgrid-1,:,:))/delthet(-ntgrid)

    if (has_phi) then
       ! this is the second term in Pi_0^{tb} in toroidal_flux.pdf notes
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g_work(:,isgn,iglo) = -zi*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo) &
                  *drhodpsi*IoB**2*gradpar*rhostar
          end do
       end do
       call get_lfflux (g_work, dphi, vflx0)

       ! this is the bracketed part of the first term in Pi_0^{tb} in toroidal_flux.pdf notes
       do iglo = g_lo%llim_proc, g_lo%ulim_proc
          do isgn = 1, 2
             g_work(:,isgn,iglo) = 0.5*gnew(:,isgn,iglo)*aj0(:,iglo)*vpa(:,isgn,iglo)**2 &
                  *drhodpsi*IoB**2*rhostar
          end do
       end do
       call get_flux (g_work, phi, vflx1)
    else
       vflx0 = 0. ; vflx1 = 0.
    end if

    deallocate (dum,dphi)
  end subroutine lf_flux

  !> FIXME : Add documentation
  subroutine get_lfflux (g_in, fld, flx)
    use theta_grid, only: ntgrid, grho, theta, bmag, gradpar
    use kt_grids, only: ntheta0, naky
    use le_grids, only: integrate_moment
    use species, only: nspec
    use mp, only: proc0
    use integration, only: trapezoidal_integration
    use gs2_layouts, only: g_lo
    implicit none
    complex, dimension (-ntgrid:,:,g_lo%llim_proc:), intent (in) :: g_in
    complex, dimension (-ntgrid:,:,:), intent (in) :: fld
    real, dimension (:,:,:), intent (out) :: flx
    complex, dimension (:,:,:,:), allocatable :: total
    real :: wgt
    integer :: ik, it, is

    allocate (total(-ntgrid:ntgrid,ntheta0,naky,nspec))
    total = 0.0
    call integrate_moment (g_in, total)

    wgt = 1.0/trapezoidal_integration(theta, grho/(bmag*gradpar))
    if (proc0) then
       do is = 1, nspec
          do ik = 1, naky
             do it = 1, ntheta0
                flx(it,ik,is) = trapezoidal_integration(theta, &
                     aimag(total(:,it,ik,is)*conjg(fld(:,it,ik))) &
                     /(bmag*gradpar))*wgt !No ky factor?
             end do
          end do
       end do

       flx = flx*0.5

    end if

    deallocate (total)

  end subroutine get_lfflux

  !> Calculate the momentum flux as a function of \((v_\parallel, \theta, t)\)
  subroutine flux_vs_theta_vs_vpa (phinew,vflx)
    use constants, only: zi
    use dist_fn_arrays, only: gnew, vperp2, aj1, aj0, vpac, g_work
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use theta_grid, only: ntgrid, bmag, gds21, gds22, qval, shat
    use theta_grid, only: Rplot, Bpol, rhoc
    use kt_grids, only: aky, theta0
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only: spec, nspec

    implicit none
    complex, dimension(-ntgrid:,:,:), intent(in) :: phinew
    real, dimension (-ntgrid:,:,:), intent (out) :: vflx

    integer :: all = 1
    integer :: iglo, isgn, ig, it, ik, is

    real, dimension (:,:,:), allocatable :: g0r
    real, dimension (:,:,:,:,:), allocatable :: gavg

    allocate (g0r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gavg(-ntgrid:ntgrid,nlambda,negrid,2,nspec))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             ! Not clear that this should use the grid-centre v|| (vpac)
             ! as this is not defined at ig = ntgrid and all other quantities
             ! are evaluated on the grid points.
             g_work(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)*vpac(ig,isgn,iglo) &
                  *Rplot(ig)*sqrt(1.0-Bpol(ig)**2/bmag(ig)**2) &
                  -zi*aky(ik)*gnew(ig,isgn,iglo)*aj1(ig,iglo) &
                  *rhoc*(gds21(ig)+theta0(it,ik)*gds22(ig))*vperp2(ig,iglo)*spec(is)%smz/(qval*shat*bmag(ig)**2)
             g0r(ig,isgn,iglo) = aimag(g_work(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)
          end do
       end do
    end do

    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)

  end subroutine flux_vs_theta_vs_vpa

  !> Diagnose particle flux contribution to toroidal momentum flux
  !> in the lab frame in terms of vpar and theta.
  subroutine pflux_vs_theta_vs_vpa (vflx)
#ifdef LOWFLOW
    use dist_fn_arrays, only: gnew, aj0, g_work
    use gs2_layouts, only: g_lo
    use gs2_layouts, only: it_idx, ik_idx, is_idx
    use theta_grid, only: Rplot !JPL
    use fields_arrays, only: phinew
    use kt_grids, only: aky
    use le_grids, only: integrate_volume, nlambda, negrid
    use le_grids, only: get_flux_vs_theta_vs_vpa
    use species, only:  nspec
    use lowflow, only: mach_lab
#endif
    use theta_grid, only: ntgrid
    implicit none

    real, dimension (-ntgrid:,:,:), intent (out) :: vflx
#ifdef LOWFLOW
    integer :: all = 1
    integer :: iglo, isgn, ig, it, ik, is

    real, dimension (:,:,:), allocatable :: g0r
    real, dimension (:,:,:,:,:), allocatable :: gavg
    allocate (g0r(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc))
    allocate (gavg(-ntgrid:ntgrid,nlambda,negrid,2,nspec))

    do iglo = g_lo%llim_proc, g_lo%ulim_proc
       it = it_idx(g_lo,iglo)
       ik = ik_idx(g_lo,iglo)
       is = is_idx(g_lo,iglo)
       do isgn = 1, 2
          do ig = -ntgrid, ntgrid
             g_work(ig,isgn,iglo) = gnew(ig,isgn,iglo)*aj0(ig,iglo)

             g0r(ig,isgn,iglo) = aimag(g_work(ig,isgn,iglo)*conjg(phinew(ig,it,ik)))*aky(ik)*Rplot(ig)**2*mach_lab(is)

          end do
       end do
    end do

    call integrate_volume (g0r, gavg, all)
    call get_flux_vs_theta_vs_vpa (gavg, vflx)

    deallocate (gavg)
    deallocate (g0r)
#else
    vflx=0.
#endif
  end subroutine pflux_vs_theta_vs_vpa

end module diagnostics_fluxes
