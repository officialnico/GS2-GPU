
!> A program that repeatedly calls xxf_lo to yxf_lo redistributes for benchmarking
!!
!! This is free software released under the MIT licence
!!   Written by: Joseph Parker (joseph.parker@stfc.ac.uk)
program time_redists
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: yxf_lo, xxf_lo, layout
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use gs2_transforms, only: x2y
  use job_manage, only: time_message
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm, proc0
  use redistribute, only: gather, scatter, write_connection_matrix, report_map_property
  use redistribute, only: get_redistname
  use run_parameters, only: nstep
  use species, only: nspec
  use theta_grid, only: ntgrid
  implicit none
  type(gs2_program_state_type) :: state
  complex, dimension (:,:), allocatable :: xxf, yxf
  real :: time_taken(2) = 0.0
  real :: time_init(2) = 0.0
  integer :: i
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  if (proc0) call announce_module_test('time_redists')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  call init(state%init, init_level_list%collisions)

  ! Set up test-specific arrays
  allocate (xxf(xxf_lo%nx,xxf_lo%llim_proc:xxf_lo%ulim_alloc)) ; xxf = 0.
  allocate (yxf(yxf_lo%ny/2+1, yxf_lo%llim_proc:yxf_lo%ulim_alloc)) ; yxf = 1.

  ! Begin test
  if (proc0) call time_message(.false., time_taken, "redist time")

  if( get_redistname(x2y)=='' ) then
    if(proc0) then
      write(*,*) "x2y has not been initialized. This probably means that the"
      write(*,*) "input file for this test uses accelx_lo rather than xxf_lo."
      write(*,*) "Now exiting cleaning, producing no results."
    end if

  else

    do i = 1,nstep
      call gather (x2y, xxf, yxf)
      call scatter (x2y, yxf, xxf)
    end do

    if (proc0) then
      call time_message(.false., time_taken, "advance time")
      write(*, '(" Time for ",I10," gather/scatters from xxf_lo to yxf_lo on ",I6," procs: ",F5.1," s")') nstep, nproc, time_taken(1)
      call append_output_file(timing_unit, &
        benchmark_identifier())
      write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
      call close_output_file(timing_unit)
    end if

    ! now timing is finishes, write info on the redist matrix
    call write_connection_matrix(x2y)

  end if

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  if (proc0) call close_module_test('time_redists')

  call finish_mp

end program time_redists
