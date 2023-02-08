
!> A program to test the replay feature
!! (see the documentation for replay in 
!! gs2_main).
!!
!! This is free software release under the 
!! MIT license. 
!! Written by:
!!     Edmund Highcock (edmundhighcock@users.sourceforge.net)

program cyclone_itg_replay
  use mp, only: init_mp, mp_comm, finish_mp
  use gs2_main
  implicit none
  type(gs2_program_state_type) :: state

  call init_mp
  state%mp_comm = mp_comm
  ! This is what we are testing
  state%replay = .true.

  call initialize_gs2(state)
  call initialize_equations(state)
  call initialize_diagnostics(state)
  state%print_times = .false.
  call evolve_equations(state, state%nstep)
  call finalize_diagnostics(state)
  call finalize_equations(state)
  state%print_times = .true.
  state%print_full_timers = .true.
  call finalize_gs2(state)

  call finish_mp

end program cyclone_itg_replay
