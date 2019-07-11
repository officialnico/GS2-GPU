
!> A program that repeatedly calls a timestep for benchmarking
!!
program time_init
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use fields, only: advance
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use run_parameters, only: nstep
  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  integer :: i
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  if (proc0) call announce_module_test('time_timestep')  

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm


  call initialize_gs2(state)
  if (proc0) call time_message(.false., time_taken, "advance time")

  call init(state%init, init_level_list%full)

  if (proc0) then
    call time_message(.false., time_taken, "advance time")
    write(*, '(" Time for initialization on ",I6," procs: ",F5.1," s")') nproc, time_taken(1)
    call append_output_file(timing_unit, &
      benchmark_identifier())
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    call close_output_file(timing_unit)
  end if

  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  if (proc0) call close_module_test('time_timestep')

  call finish_mp

end program time_init
