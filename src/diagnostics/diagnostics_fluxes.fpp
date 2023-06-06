!> Module which contains functions for calculating
!! and writing out the fluxes of heat and momentum etc.
module diagnostics_fluxes
 
  implicit none

  private

  !> Allocate arrays
  public :: init_diagnostics_fluxes

  !> Deallocate arrays
  public :: finish_diagnostics_fluxes

  !> Set averages in gnostics\%current_results to 0.
  public :: reset_averages_and_counters

  !> Calculate and possibly write fluxes.  The fluxes are actually calculated
  !!  in the module dist_fn. They are returned as a function of
  !!  x, y and species. This function writes the whole array,
  !!  and also various averages of them.
  public :: calculate_fluxes

  real, dimension (:,:,:,:), allocatable ::  qheat, qmheat, qbheat
  real, dimension (:,:,:), allocatable ::  pflux,  vflux, vflux_par, vflux_perp
  real, dimension (:,:,:), allocatable ::  pflux_tormom
  real, dimension (:,:,:), allocatable :: vflux0, vflux1  ! low flow correction to turbulent momentum flux
  real, dimension (:,:,:), allocatable :: pmflux, vmflux
  real, dimension (:,:,:), allocatable :: pbflux, vbflux
  real, dimension (:,:), allocatable :: esflux_vs_e, apflux_vs_e, bpflux_vs_e
  real, dimension (:,:,:), allocatable :: exchange
  real, dimension (:,:,:), allocatable :: exchange1

contains
  !> FIXME : Add documentation  
  subroutine init_diagnostics_fluxes(gnostics)
    use kt_grids, only: naky, ntheta0
    use species, only: nspec
    use le_grids, only: negrid
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    
    gnostics%current_results%species_heat_flux_avg = 0.0
    gnostics%current_results%species_momentum_flux_avg = 0.0
    gnostics%current_results%species_particle_flux_avg = 0.0
  
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

  !> FIXME : Add documentation  
  subroutine finish_diagnostics_fluxes
    implicit none
    deallocate (pflux, pflux_tormom, qheat, vflux, vflux_par, vflux_perp, pmflux, qmheat, vmflux, &
         pbflux, qbheat, vbflux, vflux0, vflux1, exchange, exchange1, &
         esflux_vs_e, apflux_vs_e, bpflux_vs_e)
  end subroutine finish_diagnostics_fluxes


  subroutine reset_averages_and_counters(gnostics)
    use diagnostics_config, only: diagnostics_type
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    gnostics%current_results%species_heat_flux_avg = 0.0
    gnostics%current_results%species_particle_flux_avg = 0.0
    gnostics%current_results%species_momentum_flux_avg = 0.0
  end subroutine reset_averages_and_counters

  !> FIXME : Add documentation  
  subroutine calculate_fluxes(gnostics)
    use dist_fn, only: flux, lf_flux, eexchange, flux_vs_e
    use dist_fn_arrays, only: g_adjust, gnew, to_g_gs2, from_g_gs2
    use species, only: nspec, spec
    use fields_arrays, only: phinew, bparnew, aparnew, phi
    use run_parameters, only: has_phi, has_apar, has_bpar
    use nonlinear_terms, only: nonlinear_mode_switch, nonlinear_mode_none
    use diagnostics_config, only: diagnostics_type
    use unit_tests, only: debug_message
    use gs2_io, only: starts, egrid_dim, species_dim, time_dim
    use neasyf, only: neasyf_write
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    integer :: is

    call debug_message(4, 'diagnostics_fluxes::calculate_fluxes starting')
    gnostics%current_results%total_heat_flux = 0.0
    gnostics%current_results%total_momentum_flux = 0.0
    gnostics%current_results%total_particle_flux = 0.0
    gnostics%current_results%species_heat_flux = 0.0
    gnostics%current_results%species_momentum_flux = 0.0
    gnostics%current_results%species_particle_flux = 0.0

    if (nonlinear_mode_switch .eq. nonlinear_mode_none) &
         call write_diffusive_estimates(gnostics)

    call g_adjust (gnew, phinew, bparnew, direction = from_g_gs2)
    call flux (phinew, aparnew, bparnew, &
         pflux,  qheat,  vflux, vflux_par, vflux_perp, &
         pmflux, qmheat, vmflux, pbflux, qbheat, vbflux, pflux_tormom)
    call flux_vs_e (phinew, aparnew, bparnew, esflux_vs_e, apflux_vs_e, bpflux_vs_e)

#ifdef LOWFLOW
    !lowflow terms only implemented in electrostatic limit at present
    call lf_flux (phinew, vflux0, vflux1)
#endif

    call g_adjust (gnew, phinew, bparnew, direction = to_g_gs2)
    call eexchange (phinew, phi, exchange1, exchange)
      
      if (has_phi) then
         do is = 1, nspec
            qheat(:,:,is,1) = qheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
            qheat(:,:,is,2) = qheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
            qheat(:,:,is,3) = qheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
            
            pflux(:,:,is) = pflux(:,:,is) * spec(is)%dens
            pflux_tormom(:,:,is) = pflux_tormom(:,:,is) * spec(is)%dens  

            esflux_vs_e(:,is) = esflux_vs_e(:,is) * spec(is)%dens
            
            vflux(:,:,is) = vflux(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
            vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
            vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
            exchange(:,:,is) = exchange(:,:,is) * spec(is)%dens*spec(is)%z
            
#ifdef LOWFLOW
            vflux0(:,:,is) = vflux0(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
            vflux1(:,:,is) = vflux1(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%temp/spec(is)%z
#endif
         end do

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

       if (gnostics%writing) call neasyf_write(gnostics%file_id, "es_flux_vs_e", esflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Electrostatic flux", start=starts(3, gnostics%nout))

#ifdef LOWFLOW
       call calculate_standard_flux_properties(gnostics, &
            'es_mom0', 'Low-flow momentum flux 0 (Ask Michael Barnes)', 'Pi_gB =  ', &
            vflux0, gnostics%distributed)
       call calculate_standard_flux_properties(gnostics, &
            'es_mom1', 'Low-flow momentum flux 1 (Ask Michael Barnes)', 'Pi_gB =  ', &
            vflux0, gnostics%distributed)
#endif
    end if

    if (has_apar) then
       do is = 1, nspec
          qmheat(:,:,is,1)=qmheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          qmheat(:,:,is,2)=qmheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          qmheat(:,:,is,3)=qmheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          pmflux(:,:,is)=pmflux(:,:,is) * spec(is)%dens
          vmflux(:,:,is)=vmflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm

          apflux_vs_e(:,is) = apflux_vs_e(:,is) * spec(is)%dens
       end do
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

       if (gnostics%writing) call neasyf_write(gnostics%file_id, "apar_flux_vs_e", apflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Perpendicular magnetic flux, as function of energy", &
            start=starts(3, gnostics%nout))
    end if
    if (has_bpar) then
       do is = 1, nspec
          qbheat(:,:,is,1)=qbheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
          qbheat(:,:,is,2)=qbheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
          qbheat(:,:,is,3)=qbheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
          pbflux(:,:,is)=pbflux(:,:,is) * spec(is)%dens
          vbflux(:,:,is)=vbflux(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
          bpflux_vs_e(:,is) = bpflux_vs_e(:,is) * spec(is)%dens
       end do
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
       if (gnostics%writing) call neasyf_write(gnostics%file_id, "bpar_flux_vs_e", bpflux_vs_e, &
            dim_names=[egrid_dim, species_dim, time_dim], &
            long_name="Parallel magnetic flux, as function of energy", &
            start=starts(3, gnostics%nout))
    end if
    
    ! Update averages
    gnostics%current_results%species_heat_flux_avg = gnostics%current_results%species_heat_flux_avg + &
         gnostics%current_results%species_heat_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_particle_flux_avg = gnostics%current_results%species_particle_flux_avg + &
         gnostics%current_results%species_particle_flux * (gnostics%user_time-gnostics%user_time_old)
    gnostics%current_results%species_momentum_flux_avg = gnostics%current_results%species_momentum_flux_avg + &
         gnostics%current_results%species_momentum_flux * (gnostics%user_time-gnostics%user_time_old)
    
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
  end subroutine calculate_fluxes

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
    use neasyf, only: neasyf_write
    implicit none
    type(diagnostics_type), intent(inout) :: gnostics
    real, dimension(ntheta0, naky) :: diffusivity_by_k
    !real :: diffusivity
    real :: heat_flux_max
    real, dimension(nspec) :: heat_flux_max_by_spec
    real :: particle_flux_max
    real, dimension(nspec) :: particle_flux_max_by_spec
    !real :: grho
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
    use neasyf, only: neasyf_write
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
end module diagnostics_fluxes
