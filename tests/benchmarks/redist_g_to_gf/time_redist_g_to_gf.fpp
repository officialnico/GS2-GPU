
!> A program that repeatedly calls g_lo to le_lo redistributes for benchmarking
!!
!! This is free software released under the MIT licence
!!   Written by: Joseph Parker (joseph.parker@stfc.ac.uk)
program time_redists
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: g_lo, gf_lo, layout
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use le_grids, only: nlambda, negrid, g2gf
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use redistribute, only: gather, scatter
  use run_parameters, only: nstep
  use species, only: nspec
  use theta_grid, only: ntgrid
  implicit none
  type(gs2_program_state_type) :: state
  complex, dimension (:,:,:), allocatable :: g
  complex, dimension (:,:,:,:,:,:), allocatable :: gf
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
  allocate (g(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_alloc)) ; g = 1.
  allocate (gf(-ntgrid:ntgrid, 2, nspec, negrid, nlambda, gf_lo%llim_proc:gf_lo%ulim_alloc)) ; gf = 0.

  ! Begin test
  if (proc0) call time_message(.false., time_taken, "redist time")

  do i = 1,nstep
    call gather (g2gf, g, gf, ntgrid)
    call scatter (g2gf, gf, g, ntgrid)
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
