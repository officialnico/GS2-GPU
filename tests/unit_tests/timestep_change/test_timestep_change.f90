program test_timestep_change
  use gs2_main
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  use gs2_time, only: code_dt, save_dt_cfl
  use dist_fn_arrays, only: g, gnew, gexp_1, gexp_2, gexp_3
  use fields_arrays, only: phi, apar, bpar  
  use fields_arrays, only: phinew, aparnew, bparnew
  use gs2_reinit, only: reset_time_step
  use fields, only: advance
  implicit none
  real, parameter :: tolerance = 0.0
  type(gs2_program_state_type) :: state
  complex, dimension(:, :, :) , allocatable :: g_bak, gexp_1_bak, gexp_2_bak, gexp_3_bak
  complex, dimension(:, :, :) , allocatable :: phi_bak, apar_bak, bpar_bak  
  integer :: istep
  real :: code_dt_bak
  call init_mp
  
  call announce_module_test("timestep_change")

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_wall_clock_timer 
  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)

  ! Run for five steps
  call evolve_equations(state, 5)
  istep = 5
  
  ! Make a copy of the current state
  g_bak = gnew
  gexp_1_bak = gexp_1
  gexp_2_bak = gexp_2
  gexp_3_bak = gexp_3  
  phi_bak = phinew
  apar_bak = aparnew
  bpar_bak = bparnew
  code_dt_bak = code_dt
  
  ! Change the timestep
  call save_dt_cfl(code_dt / 1.5)
  
  call prepare_initial_values_overrides(state)
  call set_initval_overrides_to_current_vals(state%init%initval_ov)
  state%init%initval_ov%override = .true.
  if (state%is_external_job) then
     call reset_time_step (state%init, istep, state%exit, state%external_job_id)
  else
     call reset_time_step (state%init, istep, state%exit)
  end if

  ! First check the restored state
  call process_test(code_dt <= code_dt_bak / 1.5, 'timestep_change : Code dt')
  call process_test(maxval(abs(gnew-g_bak)) <= tolerance, 'timestep_change : gnew')
  call process_test(maxval(abs(gexp_1-gexp_1_bak)) <= tolerance, 'timestep_change : gexp_1')
  call process_test(maxval(abs(gexp_2-gexp_2_bak)) <= tolerance, 'timestep_change : gexp_2')
  call process_test(maxval(abs(gexp_3-gexp_3_bak)) <= tolerance, 'timestep_change : gexp_3')
  call process_test(maxval(abs(phinew-phi_bak)) <= tolerance, 'timestep_change : phinew')
  call process_test(maxval(abs(aparnew-apar_bak)) <= tolerance, 'timestep_change : aparnew')
  call process_test(maxval(abs(bparnew-bpar_bak)) <= tolerance, 'timestep_change : bparnew')

  ! For now the following will fail due to a bug in the handling of istep_last in nonlinear_terms

  ! ! Now try repeating the step and checking initial values still match
  ! call advance(istep)

  ! call process_test(code_dt <= code_dt_bak / 1.5, 'timestep_change : Code dt')
  ! call process_test(maxval(abs(g-g_bak)) <= tolerance, 'timestep_change : g')
  ! call process_test(maxval(abs(gexp_1-gexp_1_bak)) <= tolerance, 'timestep_change : gexp_1 2')
  ! call process_test(maxval(abs(gexp_2-gexp_2_bak)) <= tolerance, 'timestep_change : gexp_2 2')
  ! call process_test(maxval(abs(gexp_3-gexp_3_bak)) <= tolerance, 'timestep_change : gexp_3 2')
  ! call process_test(maxval(abs(phi-phi_bak)) <= tolerance, 'timestep_change : phi')
  ! call process_test(maxval(abs(apar-apar_bak)) <= tolerance, 'timestep_change : apar')
  ! call process_test(maxval(abs(bpar-bpar_bak)) <= tolerance, 'timestep_change : bpar')
  
  call close_module_test("timestep_change")
  call finish_mp
end program test_timestep_change

