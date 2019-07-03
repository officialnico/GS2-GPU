!> A program that repeatedly calls a local iglo loop for benchmarking
!!
!! This is free software released under the MIT licence
!!   Written by: Joseph Parker (joseph.parker@stfc.ac.uk)
program time_indexing
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: g_lo, layout, is_idx, it_idx, ik_idx, il_idx, ie_idx, idx
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use run_parameters, only: nstep

  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  real :: time_init(2) = 0.0
  integer :: i, j, k, l, it, ik, il, ie, is, iglo, work
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  call announce_module_test('time_indexing')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  ! initialize up to nonlinear terms level so our loop
  ! can copy the load_kx_bpar loop
  call init(state%init, init_level_list%le_grids)

  ! Begin test
  if (proc0) call time_message(.false., time_taken, "indexing time")

  do j = 1,nstep
    do iglo=g_lo%llim_world, g_lo%ulim_world
      it=it_idx(g_lo,iglo)
      ik=ik_idx(g_lo,iglo)
      il=il_idx(g_lo,iglo)
      ie=ie_idx(g_lo,iglo)
      is=is_idx(g_lo,iglo)
      i=idx(g_lo,ik,it,il,ie,is)
    end do
  end do

  if (proc0) then
    call time_message(.false., time_taken, "indexing time")
    write(*, '(" Time for ",I6," loops in layout ",A5," on ",I6," procs: ",F5.1," s")') nstep, layout, nproc, time_taken(1)
    call append_output_file(timing_unit, &
      benchmark_identifier())
    write(timing_unit, '(I6,"   ",F9.3)') nproc, time_taken(1)
    call close_output_file(timing_unit)
  end if

  ! Without this write statement, the compiler will optimize away all the work in the loop...
  write(*,*) it, ik, il, ie, is, i

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call close_module_test('time_indexing')

  call finish_mp

end program time_indexing
