!> This unit runs a tiny grid with nstep=large and max_sim_time small
!> Checks that the run ends when the simulation time is reached, rather than number of steps.
program test_max_sim_time
  use gs2_main, only: run_gs2, gs2_program_state_type
  use unit_tests
  use mp, only: proc0, init_mp, mp_comm, finish_mp
  use job_manage, only: timer_local
  use exit_codes, only: EXIT_MAX_SIM_TIME
  use gs2_time, only: user_time, user_dt
  use run_parameters, only: max_sim_time
  implicit none
  type(gs2_program_state_type) :: state

  ! Tell make system to run this test
  functional_test_flag = .true.

  call announce_module_test("max_sim_time")

  call init_mp()
  state%mp_comm = mp_comm

  ! Run test:
  call announce_test("that gs2 halts when max_sim_time reached")
  call run_gs2(state)

  ! The simulation time should match the expected value. Note we allow
  ! for cases where the maximum time isn't hit exactly. Because of how
  ! GS2 updates the time, if we hit the max_sim_time exactly then the
  ! time on exit will read as max_sim_time + dt (essentially we update
  ! the time to the next value before we reach the next iteration)
  print*,"USER TIME",user_time, user_dt, max_sim_time
  if(proc0) call process_test((user_time <= max_sim_time + user_dt) .and. (user_time >= max_sim_time) , &
       "gs2 ends when max_sim_time reached")
  
  ! The reason for stopping should be given in "test_max_sim_time.exit_reason".
  if(proc0) call process_test(EXIT_MAX_SIM_TIME%code_matches(get_exit_code()) , &
       "gs2 ends when max_sim_time reached and generates the correct exit code")

  call close_module_test("max_sim_time")
  call finish_mp()

contains

  !> Read error code from the file "run_name.exit_reason"
  function get_exit_code()
    use file_utils, only: get_unused_unit
    implicit none
    integer :: get_exit_code
    integer :: unit

    call get_unused_unit(unit)

    open (unit=unit, file="test_max_sim_time.exit_reason")
    read (unit, '(I3.1)') get_exit_code
    close (unit)
  end function get_exit_code

end program test_max_sim_time
