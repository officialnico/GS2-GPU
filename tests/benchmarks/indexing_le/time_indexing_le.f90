!> A program that repeatedly calls an le_lo loop for benchmarking

! Copyright Joseph Parker 2019 (joseph.parker@stfc.ac.uk)

! This file is part of gs2.

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
! (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify,
! merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

program time_indexing_le
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: le_lo, layout, is_idx, it_idx, ik_idx, il_idx, ie_idx, idx
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use run_parameters, only: nstep

  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  real :: time_init(2) = 0.0
  integer :: i, j, ig, it, ik, is, ile
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  call announce_module_test('time_indexing')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  ! initialize up to le_grids terms level so indexing is initialized
  call init(state%init, init_level_list%collisions)

  ! Begin test
  if (proc0) call time_message(.false., time_taken, "indexing time")

  do j = 1,nstep
    do ile=le_lo%llim_world, le_lo%ulim_world
      ig=it_idx(le_lo,ile)
      it=it_idx(le_lo,ile)
      ik=ik_idx(le_lo,ile)
      is=is_idx(le_lo,ile)
      i=idx(le_lo,ig,ik,it,is)
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
  write(*,*) ig, it, ik, is, i

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call close_module_test('time_indexing')

  call finish_mp

end program time_indexing_le
