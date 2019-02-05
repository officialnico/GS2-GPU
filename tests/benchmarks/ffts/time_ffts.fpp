
!> A program that repeatedly calls add_nl for benchmarking the ffts and
!! transposes
!!
!! This is free software released under MIT
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
program time_ffts
  use unit_tests
#ifdef SHMEM
  use shm_mpi3
#endif
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use file_utils, only: append_output_file, close_output_file
  use kt_grids, only: naky, ntheta0
  use theta_grid, only: ntgrid
  use gs2_layouts, only: g_lo
  use nonlinear_terms, only: nonlinear_terms_unit_test_time_add_nl
  use job_manage, only: time_message
  use benchmarks, only: benchmark_identifier
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use gs2_init, only: init, init_level_list
  use run_parameters, only: nstep
  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  integer :: i
  integer :: timing_unit

  !complex, dimension (:,:,:), allocatable :: integrate_species_results
#ifdef SHMEM
  complex, dimension (:,:,:), pointer, contiguous :: g1
#else
  complex, dimension (:,:,:), allocatable :: g1
#endif
  complex, dimension (:,:,:), allocatable :: phi, apar, bpar



  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  call announce_module_test('time_ffts')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  call init(state%init, init_level_list%dist_fn_level_2)


#ifdef SHMEM
  call shm_alloc(g1,(/-ntgrid, ntgrid, 1, 2, g_lo%llim_proc, g_lo%ulim_proc/))
#else
  allocate(g1(-ntgrid:ntgrid,2,g_lo%llim_proc:g_lo%ulim_proc))
#endif
  allocate(phi(-ntgrid:ntgrid,ntheta0,naky))
  allocate(apar(-ntgrid:ntgrid,ntheta0,naky))
  allocate(bpar(-ntgrid:ntgrid,ntheta0,naky))

  if (proc0) call time_message(.false., time_taken, "FFT time")

  do i = 1,nstep
    call nonlinear_terms_unit_test_time_add_nl(g1, phi, apar, bpar)
  end do

  if (proc0) then
    call time_message(.false., time_taken, "FFT time")
    write(*, '(" Time for ",I6," nonlinear_terms on ",I6," procs: ",F5.1," s")') nstep, nproc, time_taken(1)
    call append_output_file(timing_unit, &
      benchmark_identifier())
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    call close_output_file(timing_unit)
  end if

  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call close_module_test('time_ffts')

  call finish_mp

end program time_ffts
