
!> A program to test the replay feature
!! (see the documentation for replay in 
!! gs2_main).
!!
!! This is free software release under the 
!! MIT license. 
!! Written by:
!!     Edmund Highcock (edmundhighcock@users.sourceforge.net)

program cyclone_itg_replay
  use gs2_main
  implicit none
  type(gs2_program_state_type) :: state

  ! This is what we are testing
  state%replay = .true.

  write(*,*) "initialize_gs2"
  call initialize_gs2(state)
  write(*,*) "initialize_equations"
  call initialize_equations(state)
  write(*,*) "initialize_diagnostics"
  call initialize_diagnostics(state)
  state%print_times = .false.
  write(*,*) "evolve_equations"
  call evolve_equations(state, state%nstep)
  write(*,*) "finalize_diagnostics"
  call finalize_diagnostics(state)
  write(*,*) "finalize_equations"
  call finalize_equations(state)
  state%print_times = .true.
  state%print_full_timers = .true.
  write(*,*) "finalize_gs2"
  call finalize_gs2(state)



end program cyclone_itg_replay
