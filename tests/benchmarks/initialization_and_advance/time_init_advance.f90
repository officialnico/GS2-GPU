
!> A program that separately times gs2 initialziation and the gs2 timestep.
!!
program time_init_advance
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use fields, only: advance
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm, barrier
  use run_parameters, only: nstep
  implicit none
  type(gs2_program_state_type) :: state
  real :: time_init(2) = 0.0
  real :: time_advance(2) = 0.0
  integer :: i
  integer :: timing_unit
  character(1) :: benchmark_level

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  ! Set benchmark_level variable for use in output file names
  if(proc0) then
    call get_environment_variable("BENCHMARK_LEVEL", benchmark_level)
    if (benchmark_level == "") then
      benchmark_level = "1"
    end if
  end if

  if (proc0) call announce_module_test('time_timestep')  

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  ! This initializes the state variable to the lowest initialization level
  call initialize_gs2(state)

  ! Synchronize processors for start of test
  call barrier

  if (proc0) call time_message(.false., time_init, "init time")
  ! This fully initializes gs2
  call init(state%init, init_level_list%full)

  if (proc0) then
    call time_message(.false., time_init, "init time")
    write(*, '(" Time for initialization on ",I6," procs: ",F5.1," s")') nproc, time_init(1)
    call append_output_file(timing_unit, benchmark_identifier(), "time_initialization_level_"//benchmark_level)
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_init(1)
    call close_output_file(timing_unit)
  end if

  ! Ensure tests are separate
  call barrier

  if (proc0) call time_message(.false., time_advance, "advance time")
  do i = 1,nstep
    call advance(i)
  end do

  if (proc0) then
    call time_message(.false., time_advance, "advance time")
    write(*, '(" Time for ",I6," advance steps on ",I6," procs: ",F5.1," s")') nstep, nproc, time_advance(1)
    call append_output_file(timing_unit, benchmark_identifier(), "time_advance_level_"//benchmark_level)
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_advance(1)
    call close_output_file(timing_unit)
  end if

  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  if (proc0) call close_module_test('time_timestep')

  call finish_mp

end program time_init_advance
