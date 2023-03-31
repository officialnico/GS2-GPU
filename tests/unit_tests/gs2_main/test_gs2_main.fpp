
!> This unit tests the gs2 main interface,
!! and by repeatedly initializing and finalizing
!! gs2, tests that gs2 is being properly tidied up,
!! variables deallocated etc.
!
program test_gs2_main
  use gs2_main
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  implicit none
  type(gs2_program_state_type) :: state

  call init_mp()
  state%mp_comm = mp_comm

  call announce_module_test("gs2_main")

  call initialize_gs2(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call initialize_equations(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call initialize_equations(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)
  call evolve_equations(state, state%nstep/2)
  call evolve_equations(state, state%nstep/2)
  !! This call should do nothing and print a warning
  call evolve_equations(state, state%nstep/2)
  call evolve_equations(state, state%nstep/2)

  call calculate_outputs(state)

  call finalize_diagnostics(state)
  call finalize_equations(state)
  call finalize_gs2(state)

  call run_gs2(state)
  call run_gs2(state)

  call run_gs2(state, finalize=.false.)
  call finish_gs2(state)

  call finalize_overrides(state)
  call close_module_test("gs2_main")

  call finish_mp
end program test_gs2_main
