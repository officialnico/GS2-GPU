
!> A program that repeatedly calls g_lo to le_lo redistributes for benchmarking
!!
!! This is free software released under the MIT licence
!!   Written by: Joseph Parker (joseph.parker@stfc.ac.uk)
program time_redists
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: g_lo, le_lo, layout
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use le_grids, only: nxi, negrid, g2le, check_g2le
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use redistribute, only: gather, scatter
  use run_parameters, only: nstep
  use theta_grid, only: ntgrid
  implicit none
  type(gs2_program_state_type) :: state
  complex, dimension (:,:,:), allocatable :: gle, g
  real :: time_taken(2) = 0.0
  real :: time_init(2) = 0.0
  integer :: i
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  call announce_module_test('time_redists')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  call init(state%init, init_level_list%collisions)

  ! Set up test-specific arrays
  allocate (gle(nxi+1,negrid+1,le_lo%llim_proc:le_lo%ulim_alloc)) ; gle = 0.
  allocate (g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; g = 1.

  ! Begin test
  if (proc0) call time_message(.false., time_taken, "redist time")

  !write(*,*) g2le%to_low

  do i = 1,nstep
    call gather (g2le, g, gle)
    call scatter (g2le, gle, g)
  end do

  if (proc0) then
    call time_message(.false., time_taken, "advance time")
    write(*, '(" Time for ",I6," gather/scatters from g_lo::",A5," to le_lo on ",I6," procs: ",F5.1," s")') nstep, layout, nproc, time_taken(1)
    call append_output_file(timing_unit, &
      benchmark_identifier())
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    call close_output_file(timing_unit)
  end if

  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call close_module_test('time_redists')

  call finish_mp

end program time_redists
