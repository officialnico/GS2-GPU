program test_cyclone_itg_flux
  use gs2_main, only: run_gs2, gs2_program_state_type
  use gs2_main, only: finalize_diagnostics, finalize_equations, finalize_gs2
  use unit_tests, only: announce_test, process_test, functional_test_flag
  use mp, only: test_driver_flag, init_mp, mp_comm, finish_mp
  implicit none
  character(len=*), parameter :: test_name = 'Cyclone itg flux'
  real, parameter :: absolute_tolerance = 5.0
  type(gs2_program_state_type) :: state

  test_driver_flag = .true.
  functional_test_flag = .true.

  call init_mp()
  state%mp_comm = mp_comm

  call run_gs2(state, finalize=.false.)
  call finalize_diagnostics(state)

  call announce_test('results')
  call process_test(check_fluxes(), 'results')

  call finalize_equations(state)
  call finalize_gs2(state, quiet=.true.)

  call finish_mp()

contains

  logical function check_fluxes() result(check)
    use kt_grids, only: ntheta0, naky
    use species, only: nspec, spec
    use dist_fn, only: flux
    use dist_fn_arrays, only: g_adjust, gnew
    use gs2_diagnostics, only: get_volume_average
    use run_parameters, only: fphi, fbpar
    use fields_arrays, only: phinew, aparnew, bparnew
    use mp, only: proc0
    implicit none
    real, dimension(:,:,:,:), allocatable :: qheat, qmheat, qbheat
    real, dimension(:,:,:), allocatable   :: pflux, vflux, vflux_par, vflux_perp
    real, dimension(:,:,:), allocatable   :: pmflux, vmflux, pbflux, vbflux, pflux_tormom
    real, dimension(:), allocatable :: heat_fluxes
    real, parameter :: expected_es_hflux_tot = 108.0
    integer :: is

    allocate (pflux (ntheta0,naky,nspec)) ; pflux = 0.
    allocate (pflux_tormom (ntheta0,naky,nspec)) ; pflux_tormom = 0. 
    allocate (qheat (ntheta0,naky,nspec,3)) ; qheat = 0.
    allocate (vflux (ntheta0,naky,nspec)) ; vflux = 0.

    allocate (vflux_par (ntheta0,naky,nspec)) ; vflux_par = 0.
    allocate (vflux_perp (ntheta0,naky,nspec)) ; vflux_perp = 0.

    allocate (pmflux(ntheta0,naky,nspec)) ; pmflux = 0.
    allocate (qmheat(ntheta0,naky,nspec,3)) ; qmheat = 0.
    allocate (vmflux(ntheta0,naky,nspec)) ; vmflux = 0.

    allocate (pbflux(ntheta0,naky,nspec)) ; pbflux = 0.
    allocate (qbheat(ntheta0,naky,nspec,3)) ; qbheat = 0.
    allocate (vbflux(ntheta0,naky,nspec)) ; vbflux = 0.

    allocate (heat_fluxes(nspec)) ; heat_fluxes = 0.
    call g_adjust (gnew, phinew, bparnew, fphi, fbpar)
    call flux (phinew, aparnew, bparnew, &
         pflux,  qheat,  vflux, vflux_par, vflux_perp, &
         pmflux, qmheat, vmflux, pbflux, qbheat, vbflux, pflux_tormom)

    if (proc0) then
       if (fphi > epsilon(0.0)) then
          do is = 1, nspec
             qheat(:,:,is,1) = qheat(:,:,is,1) * spec(is)%dens*spec(is)%temp
             call get_volume_average (qheat(:,:,is,1), heat_fluxes(is))

             !Commented out for now as we don't test these yet
             ! qheat(:,:,is,2) = qheat(:,:,is,2) * spec(is)%dens*spec(is)%temp
             ! call get_volume_average (qheat(:,:,is,2), heat_par(is))

             ! qheat(:,:,is,3) = qheat(:,:,is,3) * spec(is)%dens*spec(is)%temp
             ! call get_volume_average (qheat(:,:,is,3), heat_perp(is))

             ! pflux(:,:,is) = pflux(:,:,is) * spec(is)%dens
             ! call get_volume_average (pflux(:,:,is), part_fluxes(is))

             ! pflux_tormom(:,:,is) = pflux_tormom(:,:,is) * spec(is)%dens  
             ! call get_volume_average (pflux_tormom(:,:,is), part_tormom_fluxes(is))

             ! vflux(:,:,is) = vflux(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
             ! call get_volume_average (vflux(:,:,is), mom_fluxes(is))

             ! vflux_par(:,:,is) = vflux_par(:,:,is) * spec(is)%dens*sqrt(spec(is)%mass*spec(is)%temp)
             ! call get_volume_average (vflux_par(:,:,is), parmom_fluxes(is))

             ! vflux_perp(:,:,is) = vflux_perp(:,:,is) * spec(is)%dens*spec(is)%mass*spec(is)%stm
             ! call get_volume_average (vflux_perp(:,:,is), perpmom_fluxes(is))
          end do
       end if

       check = abs(sum(heat_fluxes) - expected_es_hflux_tot) <= absolute_tolerance
       if(.not. check) print*,'Fail with',sum(heat_fluxes),'but expected',expected_es_hflux_tot
    else
       check = .true. ! Just pass for iproc/=0
    end if
    

  end function check_fluxes
end program test_cyclone_itg_flux
