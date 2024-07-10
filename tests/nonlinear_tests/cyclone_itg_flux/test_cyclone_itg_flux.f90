program test_cyclone_itg_flux
  use gs2_main, only: run_gs2, gs2_program_state_type
  use gs2_main, only: finalize_diagnostics, finalize_equations, finalize_gs2
  use unit_tests, only: announce_test, process_test
  use mp, only: init_mp, mp_comm, finish_mp
  implicit none
  character(len=*), parameter :: test_name = 'Cyclone itg flux'
  real, parameter :: absolute_tolerance = 10.0
  type(gs2_program_state_type) :: state

  call init_mp()
  state%mp_comm = mp_comm

  call run_gs2(state, finalize=.false.)
  call finalize_diagnostics(state)

  call announce_test(test_name)
  call process_test(check_fluxes(), test_name)

  call finalize_equations(state)
  call finalize_gs2(state, quiet=.true.)

  call finish_mp()

contains

  logical function check_fluxes() result(check)
    use species, only: nspec
    use diagnostics_fluxes, only: init_diagnostics_fluxes, finish_diagnostics_fluxes, &
         common_calculate_fluxes, qheat
    use gs2_diagnostics, only: get_volume_average
    use run_parameters, only: has_phi
    use mp, only: proc0
    implicit none
    real, dimension (nspec) ::  heat_fluxes

    real, parameter :: expected_es_hflux_tot = 100.0
    integer :: is

    call init_diagnostics_fluxes()
    call common_calculate_fluxes()

    heat_fluxes = 0.0

    if (proc0) then
       if (has_phi) then
          do is = 1, nspec
             call get_volume_average (qheat(:,:,is,1), heat_fluxes(is))

             !Commented out for now as we don't test these yet
             ! call get_volume_average (qheat(:,:,is,2), heat_par(is))
             ! call get_volume_average (qheat(:,:,is,3), heat_perp(is))
             ! call get_volume_average (pflux(:,:,is), part_fluxes(is))
             ! call get_volume_average (pflux_tormom(:,:,is), part_tormom_fluxes(is))
             ! call get_volume_average (vflux(:,:,is), mom_fluxes(is))
             ! call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))
             ! call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))
             ! call get_volume_average (exchange(:,:,is), energy_exchange(is))
          end do
       end if

       check = abs(sum(heat_fluxes) - expected_es_hflux_tot) <= absolute_tolerance
       if(.not. check) print*,'Fail with',sum(heat_fluxes),'but expected',expected_es_hflux_tot
    else
       check = .true. ! Just pass for iproc/=0
    end if

    call finish_diagnostics_fluxes()

  end function check_fluxes
end program test_cyclone_itg_flux
